import numpy as np
import scipy.stats as st
import allel
import argparse
import pandas as pd

def read_sample_info(sample_info_file="data/sample_info_test.txt"):
    """
    Read text file with information on sample. Format of the file:
    --------------------------------------------------------------------------
    sampleID           age14C  age14Cerror  year  coverage  damageRepair  groups
    B_Ju_hoan_North-4  NA      NA           2010  40.57     TRUE          00
    S_Ju_hoan_North-1  NA      NA           2010  46.49     TRUE          00
    BallitoBayA        1980    20           NA    12.94     FALSE         11
    BallitoBayB        2110    30           NA    1.25      TRUE          12
    --------------------------------------------------------------------------

    :param sample_info_file: path of file
    :return:
    """
    # TODO : change damageRepair for ancientDamage (which should be more intuitive)
    # TODO : make it work with only ancient data (i.e. no year column)
    
    info = pd.read_table(filepath_or_buffer=sample_info_file, sep="\s+",
                         converters={'groups': lambda x: str(x)})

    sample_id = info["sampleID"]
    coverage = info["coverage"]
    is_dr = info["damageRepair"]

    sample_size = len(info)
    group_levels = len(info["groups"][1])

    groups = np.full((group_levels, sample_size), 0, "int")
    is_ancient = []
    is_modern = []
    total_ancient = 0
    for i, row in info.iterrows():
        if len(row["groups"]) != group_levels:
            print("Error: verify that all individuals are assigned to groups")
            # TODO : interrupt program here
        if np.isnan(row["year"]):
            is_ancient.append(True)
            is_modern.append(False)
            total_ancient = total_ancient + 1
        else:
            is_ancient.append(False)
            is_modern.append(True)
        for level in range(0, group_levels):
            groups[level, i] = row["groups"][level]


    return sample_id, coverage, is_ancient, is_modern, is_dr, total_ancient, \
        sample_size, group_levels, groups


# TODO : use pandas for reading table


def read_genome_intervals(genome_info_file="data/genome_info_test.txt"):
    """
    read file with information on the starting position and end position of
    each chromosome arm (centromeres removed)

    :param genome_info_file:
    :return:
    """
    genome_file = open(genome_info_file, "r")
    next(genome_file)  # header_line = next(info_file) # if header needed
    start = []
    end = []
    rates = []
    for line in genome_file:
        # print(line.split())
        v1, v2, v3, v4, v5, v6 = line.split()
        start.append(int(v2))
        start.append(int(v5))
        end.append(int(v4))
        end.append(int(v3))
        rates.append(float(v6))
        rates.append(float(v6))
        # TODO: check for blank lines at the end of the file
    return start, end, rates


# TODO : use pandas for reading table


def snp_calling(true_genotype, f_num_reads, error_rate=0.005, reads_th=1,
                score_th=10, ratio_th=3, dr=True, transversion=True):
    """
    snp_calling function takes perfect simulated data from one locus of one 
    diploid individual and adds missing data and error according to the number 
    of reads of the site, error rate of the sequencing technology and, for 
    ancient DNA not sequenced from damage repair (dr) libraries, creates 
    missing data for transition SNPs (since they cannot be distinguished from
    aDNA damage)

    :param true_genotype:
    :param f_num_reads:
    :param error_rate:
    :param reads_th:
    :param score_th:
    :param ratio_th:
    :param dr:
    :param transversion:
    :return:
    """
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


def empty_genotype_array(n_loci, n_samples, ploidy=2, allele=-1):
    """
    Creates a genotype array with all values as missing (-1) for a given number
    of samples, loci and ploidy

    :return: empty_ga
    """
    empty_ga = allel.GenotypeArray(np.full((n_loci, n_samples, ploidy), allele), dtype='i1')
    return empty_ga


def sequencing(ts, ssize, ttr, seq_error, dr, cov):
    if len(cov) != ssize:
        msg = "Number of coverage values (length=" + str(len(cov)) + \
              ") and number of samples (ssize=" + str(ssize) + \
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
        transversion_snp = True
        if np.random.random() < ttr / (ttr + 1):
            transversion_snp = False
        # print(num_reads)
        for i in range(0, 2 * ssize, 2):
            geno_data[locus, int(i / 2)] = snp_calling(true_genotype=var_genotypes[i:(i + 2)],
                                                       f_num_reads=num_reads[int(i / 2)],
                                                       error_rate=seq_error,
                                                       dr=dr[int(i / 2)],
                                                       transversion=transversion_snp)
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
        options = parser.parse_args(['-i', '../data/sample_info_test.txt',
                                     '-g', 'data/genome_test.txt',
                                     '-s', '1',
                                     '-b', '1',
                                     '-p', 'test',
                                     '-t', '0', '22', '46', '71', '78', '85', '119', '146', '208', '290', '305', '384',
                                     '-z', '4', '1', '1', '2', '2', '1', '1', '1', '1', '1', '1', '1',
                                     '-o', '0', '1', '2', '3', '10', '9', '4', '8', '5', '6', '7', '12', '13', '16',
                                     '11', '15', '14',
                                     '-d', '1762431206',
                                     '-n', '124',
                                     '-u', '5e-08'])
    else:
        options = parser.parse_args()
    if len(options.ts) != len(options.ss):
        msg = "Number of samples (length of ss=" + str(len(options.ss)) + \
              ") and number of sampling times (length of st=" + str(len(options.ts)) + \
              ") do not match"
        raise ValueError(msg)
    return options
