#!/usr/bin/python3.6

# python3 src/add.mutations.py -i data/SampleInfoTest.txt -s 1 -b 1 -p test
# -t 0 22 46 71 78 85 119 146 208 290 305 384
# -z 4 1 1 2 2 1 1 1 1 1 1 1
# -d 123456789 -n 200

import allel
import msprime
import numpy as np
import pyslim
import scipy.stats as st
import argparse


def read_sample_info(sample_info_file="data/SampleInfoTest.txt"):
    '''
    Read csv file with information on sample. Format of the file:
    --------------------------------------------------------------------------
    sampleID           age14C  age14Cerror  year  coverage  damageRepair  group
    B_Ju_hoan_North-4  NA      NA           2010  40.57     TRUE          00
    S_Ju_hoan_North-1  NA      NA           2010  46.49     TRUE          00
    BallitoBayA        1980    20           NA    12.94     FALSE         11
    BallitoBayB        2110    30           NA    1.25      TRUE          12
    --------------------------------------------------------------------------

    :param sample_info_file: path of file
    :return:
    '''
    info_file = open(sample_info_file, "r")
    next(info_file) # header_line = next(info_file) # if header needed
    sample_id = []
    # age14C = []
    # age14Cerror = []
    # year = []
    coverage = []
    is_dr = []
    is_ancient = []
    is_modern = []
    total_ancient = 0
    sample_size = 0
    # TODO: check for blank lines at the end of the file
    first_line = True
    for line in info_file:
        sample_size = sample_size + 1
        v1, v2, v3, v4, v5, v6, v7 = line.split()
        if (first_line):
            group_levels = len(v7)
        sample_id.append(v1)
        if v2 == "NA":
            # age14C.append(float('nan'))
            is_modern.append(True)
        else:
            # age14C.append(float(v2))
            is_modern.append(False)
        # if v3=="NA":
        # age14Cerror.append(float('nan'))
        # else:
        # age14Cerror.append(float(v3))
        if v4 == "NA":
            # year.append(float('nan'))
            is_ancient.append(True)
            total_ancient = total_ancient + 1
        else:
            # year.append(int(v4))
            is_ancient.append(False)
        coverage.append(float(v5))
        if (v6 == "FALSE") or (v6 == "F"):
            is_dr.append(False)
        elif (v6 == "TRUE") or (v6 == "T"):
            is_dr.append(True)
        else:
            is_dr.append(False)
        first_line = False
    info_file.seek(0)
    next(info_file)
    groups = np.full((group_levels,sample_size),0,"int")
    sample_counter = 0
    for line in info_file:
        v1, v2, v3, v4, v5, v6, v7 = line.split()
        for level in range(0,group_levels):
            groups[level,sample_counter] = v7[level]
        sample_counter = sample_counter + 1

    # t0 = max(year)

    return sample_id, coverage, is_ancient, is_modern, is_dr, total_ancient,\
           sample_size, group_levels, groups


def read_recombination_map(recombination_map_file="data/recombination_map_msprime.txt"):
    '''
    Read file with positions for start and end of chromosomes and
    recombination rates for each chromosome. Format:
    -------------------------------------
    0 1.14856e-08
    249218992 0.5
    249218993 1.10543e-08
    492309989 0.5
    492309990 1.12796e-08
    690184517 0.5
    -------------------------------------

    :param recombination_map_file:
    :return:
    '''
    file_recomb_map = open(recombination_map_file, "r")
    positions = []
    rates = []
    for line in file_recomb_map:
        p, r = line.split()
        positions.append(int(p))
        rates.append(float(r))
    return positions, rates


def read_genome_intervals(genome_info_file="data/genome_test.txt"):
    '''
    read file with information on the starting position and end position of
    each chromosome arm (centromeres removed)

    :param genome_info_file:
    :return:
    '''
    genome_file = open(genome_info_file, "r")
    next(genome_file) # header_line = next(info_file) # if header needed
    start = []
    end = []
    rates = []
    for line in genome_file:
        #print(line.split())
        v1, v2, v3, v4, v5, v6 = line.split()
        start.append(int(v2))
        start.append(int(v5))
        end.append(int(v4))
        end.append(int(v3))
        rates.append(float(v6))
        rates.append(float(v6))
        # TODO: check for blank lines at the end of the file
    return start, end, rates


def snp_calling(true_genotype, f_num_reads, error_rate=0.005, reads_th=1, score_th=10,
                ratio_th=3, dr=True, transversion=True):
    '''
    snp_calling function takes perfect simulated data from one locus of one diploid individual and
    adds missing data and error according to the number of reads of the site, error rate of the
    sequencing technology and, for ancient DNA not sequenced from damage repair (dr) libraries,
    creates missing data for transition SNPs (since they cannot be distinguished from aDNA damage)

    :param true_genotype:
    :param f_num_reads:
    :param error_rate:
    :param reads_th:
    :param score_th:
    :param ratio_th:
    :param dr:
    :param transversion:
    :return:
    '''
    if dr is False and transversion is False:
        genotype_call = [-1, -1]
    elif f_num_reads >= reads_th:
        derived_count = sum(true_genotype)
        p_derived = derived_count / 2. * (1 - error_rate) + (1 - derived_count / 2.) * error_rate
        derived_reads = st.binom.rvs(f_num_reads, p_derived)
        ancestral_reads = f_num_reads - derived_reads
        if f_num_reads >= score_th:
            if derived_reads == 0:
                genotype_call = [0, 0]
            elif ancestral_reads == 0:
                genotype_call = [1, 1]
            else:
                ratio_of_scores = derived_reads / ancestral_reads
                if (ratio_of_scores >= 1 / ratio_th) & (ratio_of_scores <= ratio_th):
                    genotype_call = [0, 1]
                elif derived_reads > ancestral_reads:
                    genotype_call = [1, 1]
                else:
                    genotype_call = [0, 0]
        else:
            random_allele = st.binom.rvs(1, derived_reads / f_num_reads)
            genotype_call = np.full(2, random_allele)
    else:
        genotype_call = [-1, -1]
    return genotype_call


def empty_genotype_array(n_loci,n_samples,ploidy = 2,allele = -1):
    '''
    Creates a genotype array with all values as missing (-1) for a given number
    of samples, loci and ploidy

    :return: empty_ga
    '''
    empty_ga = allel.GenotypeArray(np.full((n_loci, n_samples, ploidy), allele), dtype='i1')
    return empty_ga

def sequencing(ts,ssize,ttr,seq_error,dr,cov):
    if len(cov) != ssize:
        msg="Number of coverage values (length="+ str(len(cov)) +\
            ") and number of samples (ssize="+ str(ssize) +\
            ") do not match"
        raise ValueError(msg)

    geno_data = empty_genotype_array(n_loci=ts.num_mutations,
                                     n_samples=ssize,
                                     ploidy=2)
    positions = []
    locus = 0
    for variant in ts.variants():
        positions.append(round(variant.position))
        # print(positions)
        var_genotypes = variant.genotypes
        # print("--------------------------------")
        # print(var_genotypes)
        num_reads = np.random.poisson(lam=cov, size=ssize)
        transversion_SNP = True
        if np.random.random() < ttr / (ttr + 1):
            transversion_SNP = False
        # print(num_reads)
        for i in range(0, 2 * ssize, 2):
            geno_data[locus, int(i / 2)] = snp_calling(true_genotype=var_genotypes[i:(i + 2)],
                                                       f_num_reads=num_reads[int(i / 2)],
                                                       error_rate=seq_error,
                                                       dr=dr[int(i / 2)],
                                                       transversion=transversion_SNP)
        locus = locus + 1
    return geno_data, positions


def get_arguments(interactive=False):
    parser = argparse.ArgumentParser(description='Gets SLiM tree sequences, recapitates, add mutations '
                                                 'and calculates summary statistics')
    parser.add_argument('-b', '--batch_number',
                        dest='batch_id',
                        required=True,
                        type=int,
                        help='[type: %(type)s] Number used to identify the batch of simulations. '
                             'It is used as part of the input files (coming from SLiM)')
    parser.add_argument('-d', '--seed',
                        dest='seed',
                        required=True,
                        type=int,
                        help='[type: %(type)s] Seed for random number generator')
    parser.add_argument('-e', '--sequencing_error',
                        dest='seq_error',
                        default=0.005,
                        type=float,
                        help='[type: %(type)s] Sequencing error rate. '
                             '[default: %(default)s]')
    parser.add_argument("-g", "--genome_info_file",
                        dest='genome_file',
                        type=str,
                        required=True,
                        help='[type: %(type)s] Text file with genome information organised as in the '
                             'example below:\n'
                             '--------------------------------------------------------------------\n'
                             'ID chromosome_start chromosome_end centromere_start centromere_end recombination_rate\n'
                             'chr1	0	249218991	121700000	125099999	1.15E-08\n'
                             'chr2	249218992	492309988	341018992	345218991	1.11E-08\n'
                             'chr3	492309989	690184516	580109989	586309988	1.13E-08\n'
                             'chr4	690184517	881213161	738384517	741984516	1.12E-08\n'
                             'chr5	881213162	1061928971	927313162	932613161	1.13E-08\n'
                             '--------------------------------------------------------------------')
    parser.add_argument("-i", "--sample_info_file",
                        dest='info_file',
                        type=str,
                        required=True,
                        help='[type: %(type)s] Text file with sample information organised as in the '
                             'example below:\n'
                             '--------------------------------------------------------------------\n'
                             'sampleID           age14C  age14Cerror  year  coverage  damageRepair\n'
                             'B_Ju_hoan_North-4  NA      NA           2010  40.57     TRUE\n'
                             'S_Ju_hoan_North-1  NA      NA           2010  46.49     TRUE\n'
                             'BallitoBayA        1980    20           NA    12.94     FALSE\n'
                             'BallitoBayB        2110    30           NA    1.25      TRUE\n'
                             '--------------------------------------------------------------------')
    parser.add_argument('-k', '--trans_transv_ratio',
                        dest='ttratio',
                        default=2.0,
                        type=float,
                        help='[type: %(type)s] Transition/transversion ratio in the studied species. '
                             'This value is used to simulate the '
                             'filtering of transition polymorphisms from ancient DNA as they are not '
                             'distinguishable from DNA damage (except for damage repair libraries). '
                             '[default: %(default)s]')
    parser.add_argument('-n', '--effective_population_size',
                        dest='ne',
                        required=True,
                        type=int,
                        help='[type: %(type)s] Effective population size used for recapitation.')
    parser.add_argument('-o', '--sample_order',
                        dest='sample_order',
                        required=True,
                        type=int,
                        nargs='*',
                        help='[type: %(type)s] Sample chronological order in the simulation '
                             '(tree sequence). It will be different from the order in the sample '
                             'information file.')
    parser.add_argument('-p', '--project_name',
                        dest='project',
                        required=True,
                        type=str,
                        help='[type: %(type)s] Name of the project analysis. '
                             'It is used as directory in the path to the input files (coming from SLiM)')
    parser.add_argument('-s', '--simulation_number',
                        dest='sim_i',
                        required=True,
                        type=int,
                        help='[type: %(type)s] Number used to identify the simulation within a batch. '
                             'It is used as part of the input files (coming from SLiM)')
    parser.add_argument('-t', '--time_samples',
                        dest='ts',
                        required=True,
                        type=int,
                        nargs='*',
                        help='[type: %(type)s] Time (backwards) of samples, in generations.')
    parser.add_argument('-u', '--mutation_rate',
                        dest='mu',
                        default=1.25e-08,
                        type=float,
                        help='[type: %(type)s] Mutation rate per bp. '
                             '[default: %(default)s]')
    parser.add_argument('-z', '--sample_size',
                        dest='ss',
                        required=True,
                        type=int,
                        nargs='*',
                        help='[type: %(type)s] Sample size, in number of diploid individuals.')
    if interactive:
        options = parser.parse_args(['-i', 'data/SampleInfoTest.txt',
                                     '-g', 'data/genome_test.txt',
                                     '-s', '1',
                                     '-b', '1',
                                     '-p', 'test',
                                     '-t', '0','22','46','71','78','85','119','146','208','290','305','384',
                                     '-z', '4','1','1','2','2','1','1','1','1','1','1','1',
                                     '-o', '0','1','2','3','10','9','4','8','5','6','7','12','13','16','11','15','14',
                                     '-d', '1762431206',
                                     '-n', '124',
                                     '-u', '5e-08'])
    else:
        options = parser.parse_args()
    if len(options.ts) != len(options.ss):
        msg="Number of samples (length of ss="+ str(len(options.ss)) +\
            ") and number of sampling times (length of st="+ str(len(options.ts)) +\
            ") do not match"
        raise ValueError(msg)
    return options

############################################################################################################
############################################################################################################
def main():
    # Get options from command line arguments and info from input file
    options = get_arguments()
    options = get_arguments(interactive=True)
    sample_id, coverage, is_ancient, is_modern, is_dr, total_ancient, \
    sample_size, group_levels, groups = read_sample_info(sample_info_file=options.info_file)
    # print(groups)

    # initial settings and verifications
    np.random.seed(options.seed)
    na = len(options.ss) - 1
    if sum(options.ss) != sample_size:
        msg = "Number of samples from command line (sum of ss=" + str(sum(options.ss)) + \
              ") and number of samples from file (sample_size=" + str(sample_size) + \
              ") do not match"
        raise ValueError(msg)

    # add demography here
    # demogr_event = [msprime.PopulationParametersChange(time=1000, initial_size=300, population_id=0)]

    # read tree sequence from SLiM output file
    treesq = pyslim.load(
        "results/" + options.project + "/" + str(options.batch_id) + "/slim" + str(options.sim_i) + ".tree")
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
    start_chr_arm, end_chr_arm, rec_rate_chr_arm = read_genome_intervals()
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
        #sample_individuals = np.empty(0,dtype=int)
        #for x in range(0, na + 1):
        #    sample_individuals = np.concatenate([sample_individuals, gi_treesq.individuals_alive_at(options.ts[x])])
        #for ind in sample_individuals:
        #    print(gi_treesq.individual(ind))

        chrono_order_coverage = [coverage[i] for i in options.sample_order]
        chrono_order_is_dr = [is_dr[i] for i in options.sample_order]

        print("Number of mutations " + str(gi_treesq.num_mutations))
        if gi_treesq.num_mutations == 0:
            print("No mutations")
            # TODO: Create empty sumstats
        else:
            geno_data, positions = sequencing(ts=gi_treesq,
                                              ssize=sample_size,
                                              ttr=options.ttratio,
                                              seq_error=options.seq_error,
                                              dr=chrono_order_is_dr,
                                              cov=chrono_order_coverage)

            # print("Genotype data matrix:")
            # print(geno_data)
            ac = geno_data.count_alleles()
            # print("Allele count=\n"+str(ac))
            pi = allel.mean_pairwise_difference(ac) / length_interval
            mean_pi = sum(pi) / length_interval

            print("pi = " + str(mean_pi) + " for genome interval " + str(gi))

    # ac = geno_data.count_alleles()
    # pi = allel.sequence_diversity(positions, ac, start=1, stop=1000000)
    # print(pi)

    # pi_w, windows, n_bases, counts = allel.windowed_diversity(positions, ac, size=100000, start=1, stop=1000000)
    # n_obs, minmax, mean, var, skew, kurt = st.describe(pi_w)
    # print(n_obs)
    # print(minmax)
    # print(mean)
    # print(var)
    # print(skew)
    # print(kurt)


############################################################################################################
############################################################################################################
if __name__ == "__main__":
    main()


