#!/usr/bin/python3.6

import allel
import msprime
import numpy as np
import pyslim
import scipy.stats as st
import pandas as pd

#
#
def read_sample_info(sample_info_file="data/SampleInfoTest.csv"):
    info = pd.read_csv(sample_info_file,skipinitialspace=True)
    sampleID = info.sampleID
    age14C = info.age14C
    age14Cerror = info.age14Cerror
    ageBCAD = info.year
    coverage = info.coverage
    is_ancient = np.isnan(info.year)
    is_modern = np.isnan(info.age14C)
    is_dr = info.damageRepair
    total_ancient = sum(is_ancient)
    size = len(info)
    t0 = max(info.year)
    return sampleID, age14C, age14Cerror, ageBCAD, coverage, is_ancient, is_modern, is_dr, total_ancient, sample_size, t0


sampleID, age14C, age14Cerror, ageBCAD, coverage, is_ancient, is_modern, is_dr, total_ancient, sample_size, t0 = read_sample_info()


# snp_calling function takes perfect simulated data from one locus of one diploid individual and
# adds missing data and error according to the number of reads of the site, error rate of the
# sequencing technology and, for ancient DNA not sequenced from damage repair (dr) libraries,
# creates missing data for transition SNPs (since they cannot be distinguished from aDNA damage)
def snp_calling(true_genotype, f_num_reads, error_rate=0.005, reads_th=1, score_th=10,
                ratio_th=3, dr=True, transversion=True):
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


ttratio = 2.0/1.0
# i     = int(sys.argv[1])
i = 1
batch_ID = 1
# N     = int(sys.argv[2])
N = 200
# mu    = float(sys.argv[3])
mu = 5e-9  # 1.25e-08
# seed  = int(sys.argv[4])
seed = 123456789
# na    = int(sys.argv[6])
na = 11
# ts = ()
# for x in range(0, na+1):
#  ts=ts+( int(sys.argv[7+x]), )
ts = (0,  22,  46,  71,  78,  85, 119, 146, 208, 290, 305, 384)
# ss = ()
# for x in range(0, na+1):
#  ss=ss+( int(sys.argv[7+na+x]), )
ss = (4,  1,  1,  2,  2,  1,  1,  1,  1,  1,  1,  1)

coverage = (40.570000, 46.490000, 34.220000, 37.240000, 12.940000,
            1.250000, 0.010000, 1.100000, 2.300000, 0.800000, 1.000000,
            11.096470, 0.227619, 11.581280, 4.998657, 0.130293, 22.193640)


sequencing_error = 0.005

np.random.seed(seed)

# read recombination map and load it for msprime model
file_recomb_map = open("data/recombination_map_msprime.txt", "r")
positions = []
rates = []
for line in file_recomb_map:
    p, r = line.split()
    positions.append(int(p))
    rates.append(float(r))

recomb_map = msprime.RecombinationMap(positions=positions, rates=rates)  # num_loci=positions[-1])

# add demography here
demogr_event = [msprime.PopulationParametersChange(time=1000, initial_size=300, population_id=0)]

# read tree, recapitate & add mutations
treesq = pyslim.load("results/" + str(batch_ID) + "/slim" + str(i) + ".tree")
treesq = treesq.recapitate(Ne=N, recombination_map=recomb_map, demographic_events=demogr_event, model="dtwf")

sample_ind = np.random.choice(treesq.individuals_alive_at(ts[0]), ss[0], replace=False)
for x in range(1, na + 1):
    sample_ind = np.concatenate([sample_ind, treesq.individuals_alive_at(ts[x])])

keep_nodes = []
for i in sample_ind:
    keep_nodes.extend(treesq.individual(i).nodes)

treesq = treesq.simplify(keep_nodes)

# check age of sample:  treesq.individual(0).time

sample_ind = treesq.individuals_alive_at(ts[0])
for x in range(1, na + 1):
    sample_ind = np.concatenate([sample_ind, treesq.individuals_alive_at(ts[x])])


treesq = pyslim.SlimTreeSequence(msprime.mutate(treesq, rate=mu))

ploidy = 2
n_samples = sample_ind.shape[0]
n_loci = treesq.num_mutations

geno_data = allel.GenotypeArray(np.full((n_loci, n_samples, ploidy), -1), dtype='i1')

# check age of sample:  treesq.individual(0).time
positions = []
locus = 0
for variant in treesq.variants():
    positions.append(round(variant.position))
    # print(positions)
    var_genotypes = variant.genotypes
    # print(var_genotypes)
    num_reads = st.poisson.rvs(mu=coverage, size=n_samples)
    transversion_SNP = True
    if st.uniform.rvs() < ttratio/(ttratio+1):
        transversion_SNP = False
    # print(num_reads)
    for i in range(0, 2 * n_samples, 2):
        geno_data[locus, int(i / 2)] = snp_calling(true_genotype=var_genotypes[i:(i + 2)],
                                                   f_num_reads=num_reads[int(i/2)],
                                                   transversion=transversion_SNP)
    locus = locus + 1


ac = geno_data.count_alleles()
pi = allel.sequence_diversity(positions, ac, start=1, stop=1000000)
print(pi)

pi_w, windows, n_bases, counts = allel.windowed_diversity(positions, ac, size=100000, start=1, stop=1000000)
n_obs, minmax, mean, var, skew, kurt = st.describe(pi_w)
print(n_obs)
print(minmax)
print(mean)
print(var)
print(skew)
print(kurt)
