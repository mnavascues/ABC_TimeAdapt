#!/usr/bin/python3.6

# python3 msprimeNstats.py -i data/SampleInfoTest.txt -s 1 -b 1 -p test
# -t 0 22 46 71 78 85 119 146 208 290 305 384
# -z 4 1 1 2 2 1 1 1 1 1 1 1
# -d 123456789 -n 200

import allel
import msprime
import numpy as np
import pyslim
import myfun


def main():
    # Get options from command line arguments and info from input file
    options = myfun.get_arguments()
    # options = myfun.get_arguments(interactive=True)
    sample_id, coverage, is_ancient, is_modern, is_dr, total_ancient, \
        sample_size, group_levels, \
        groups = myfun.read_sample_info(sample_info_file=options.info_file)
    # print(groups)

    # initial settings and verifications
    np.random.seed(options.seed)
    na = len(options.ss) - 1
    if sum(options.ss) != sample_size:
        msg = "Number of samples from command line (sum of ss=" + str(sum(options.ss)) + \
              ") and number of samples from file (sample_size=" + str(sample_size) + \
              ") do not match"
        raise ValueError(msg)

    chrono_order_coverage = [coverage[i] for i in options.sample_order]
    chrono_order_is_dr = [is_dr[i] for i in options.sample_order]

    chrono_order_groups = np.zeros([group_levels, sample_size], dtype='int')
    groups_in_level = {}
    num_of_pair_comparisons = 0
    number_of_groups = np.zeros(group_levels, dtype='int')
    total_number_of_groups = 0
    for lev in range(0, group_levels):
        chrono_order_groups[lev] = [groups[lev][i] for i in options.sample_order]
        number_of_groups[lev] = len(np.unique(groups[lev]))
        total_number_of_groups += number_of_groups[lev]
        num_of_pair_comparisons += int((number_of_groups[lev] * (number_of_groups[lev] - 1)) / 2)
        for g in range(0, number_of_groups[lev]):
            groups_in_level['level' + str(lev) + 'group' + str(g)] = np.where(chrono_order_groups[lev] == g)[0]

    # TODO mark redundant groups between levels so summary statistics are not calculated twice
    # unique_groups = {}
    # for g in range(0,total_number_of_groups):
        


    # add demography here
    # demogr_event = [msprime.PopulationParametersChange(time=1000, initial_size=300, population_id=0)]

    # read tree sequence from SLiM output file
    treesq = pyslim.load(
        "results/" + options.project + "/" + str(options.batch_id) + "/slim_" + str(options.sim_i) + ".tree")
    # tree = treesq.first()
    # print(tree.draw(format="unicode"))

    # simplify tree sequence keeping nodes for the sampled individuals and their roots
    sample_individuals = np.random.choice(treesq.individuals_alive_at(options.ts[0]),
                                          options.ss[0], replace=False)
    sample_individuals.sort()
    for x in range(1, na + 1):
        sample_individuals = np.concatenate([sample_individuals, treesq.individuals_alive_at(options.ts[x])])
    keep_nodes = []
    for samp_i in sample_individuals:
        keep_nodes.extend(treesq.individual(samp_i).nodes)
    treesq = treesq.simplify(keep_nodes, keep_input_roots=True)

    # read genome intervals (e.g. chromosome arms start and end)
    # and recombination rates from file
    start_chr_arm, end_chr_arm, \
        rec_rate_chr_arm = myfun.read_genome_intervals(genome_info_file=options.genome_file)
    num_of_genome_intervals = len(start_chr_arm)
    # num_of_genome_intervals = 1 # KEEP THIS LINE ONLY FOR TEST

    # loop over genome intervals: recapitate, mutate, calculate sumstats
    for gi in range(0, num_of_genome_intervals):
        length_interval = end_chr_arm[gi] - start_chr_arm[gi]
        genome_interval = np.array([[start_chr_arm[gi], end_chr_arm[gi]]])
        # genome_interval = np.array([[0, 100000]]) # KEEP THIS LINE ONLY FOR TEST
        gi_treesq = treesq.keep_intervals(genome_interval, simplify=False)
        gi_treesq = gi_treesq.ltrim()
        gi_treesq = pyslim.SlimTreeSequence(gi_treesq.rtrim())
        gi_treesq = gi_treesq.recapitate(recombination_rate=rec_rate_chr_arm[0],
                                         Ne=options.ne,
                                         # demographic_events=demogr_event,
                                         model="dtwf",
                                         random_seed=np.random.randint(1, 2 ^ 32 - 1))
        gi_treesq = pyslim.SlimTreeSequence(msprime.mutate(gi_treesq,
                                                           rate=options.mu,
                                                           random_seed=np.random.randint(1, 2 ^ 32 - 1)))
        #
        # sample_individuals = np.empty(0,dtype=int)
        # for x in range(0, na + 1):
        #    sample_individuals = np.concatenate([sample_individuals, gi_treesq.individuals_alive_at(options.ts[x])])
        # for ind in sample_individuals:
        #    print(gi_treesq.individual(ind))

        print("Number of mutations " + str(gi_treesq.num_mutations))
        if gi_treesq.num_mutations == 0:
            print("No mutations")
            # TODO: Create empty sumstats
        else:
            geno_data, positions = myfun.sequencing(ts=gi_treesq,
                                                    ssize=sample_size,
                                                    ttr=options.ttratio,
                                                    seq_error=options.seq_error,
                                                    dr=chrono_order_is_dr,
                                                    cov=chrono_order_coverage)

            # print("Genotype data matrix:")
            # print(geno_data)
            allele_counts = geno_data.count_alleles()
            # print("Allele count=\n"+str(ac))
            # TODO : make a function to calculate single sample summary stats, call it here and below
            pi_per_lg = allel.sequence_diversity(positions, allele_counts,
                                                 start=1,
                                                 stop=length_interval)
            he_per_site = allel.mean_pairwise_difference(allele_counts)
            pi_per_window, windows, n_bases, \
            n_sites = allel.windowed_diversity(positions, allele_counts,
                                               size=50000,
                                               start=1,
                                               stop=length_interval)
            allele_counts_per_group = {}
            for lev in range(0, group_levels):
                for g in range(0, number_of_groups[lev]):
                    allele_counts_per_group['level' + str(lev) + 'group' + str(g)] = \
                        geno_data.count_alleles(subpop=groups_in_level['level' + str(lev) + 'group' + str(g)])
                    # TODO : call here function for single sample sumstats

            for lev in range(0, group_levels):
                for g in range(0, number_of_groups[lev]):
                    for h in range(g + 1, number_of_groups[lev]):
                        gr1 = allele_counts_per_group['level' + str(lev) + 'group' + str(g)]
                        gr2 = allele_counts_per_group['level' + str(lev) + 'group' + str(h)]
                        pairwise_diff = allel.mean_pairwise_difference_between(gr1, gr2)
                        # print("Level: " + str(lev) + ". Groups: " + str(g) + " " + str(h) +
                        #      ". Pairwise difference: " + str(pairwise_diff))

    outfile = open("results/" + options.project + "/" + str(options.batch_id) + "/stats_" + str(options.sim_i) + ".txt",
                   "w")


# n_obs, minmax, mean, var, skew, kurt = st.describe(pi_w)

############################################################################################################
############################################################################################################
if __name__ == "__main__":
    main()
