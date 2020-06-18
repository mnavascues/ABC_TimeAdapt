import allel
import msprime
import numpy as np
import pyslim
import scipy.stats as st


def snp_calling(true_genotype, f_num_reads, error_rate=0.005, reads_th=1, score_th=30, ratio_th=3):
    if f_num_reads >= reads_th:
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


# i     = int(sys.argv[1])
i = 1
# N     = int(sys.argv[2])
N = 200
# mu    = float(sys.argv[3])
mu = 5e-9
# seed  = int(sys.argv[4])
seed = 1234
# na    = int(sys.argv[6])
na = 11
# ts = ()
# for x in range(0, na+1):
#  ts=ts+( int(sys.argv[7+x]), )
ts = (0,  22,  47,  71,  78,  86, 119, 146, 208, 291, 305, 384)
# ss = ()
# for x in range(0, na+1):
#  ss=ss+( int(sys.argv[7+na+x]), )
ss = (12,  1,  1,  2,  2,  1,  1,  1,  1,  1,  1,  1)

coverage = (33., 34., 33., 34., 33., 34., 33., 34., 33., 34.,
            33., 34., 12.94,  1.25,  0.01,  1.1,  2.3,  0.8,  1., 11.096470,
            0.227619, 11.58128, 4.998657,  0.130293, 22.19364031)

sequencing_error = 0.005

np.random.seed(seed)

recomb_map = msprime.RecombinationMap.read_hapmap("data/RecombinationMap2010/genetic_map_GRCh37_chr22.txt")
demogr_event = [msprime.PopulationParametersChange(time=1000, initial_size=300, population_id=0)]

# read tree, recapitate & add mutations
treesq = pyslim.load("results/slim" + str(i) + ".tree")
treesq = treesq.recapitate(Ne=N, recombination_map=recomb_map, demographic_events=demogr_event)

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
    #print(positions)
    var_genotypes = variant.genotypes
    #print(var_genotypes)
    num_reads = st.poisson.rvs(mu=coverage, size=n_samples)
    #print(num_reads)
    for i in range(0, 2 * n_samples, 2):
        geno_data[locus, int(i / 2)] = snp_calling(true_genotype=var_genotypes[i:(i + 2)], f_num_reads=num_reads[int(i/2)])
    locus = locus + 1




ac = geno_data.count_alleles()
pi = allel.sequence_diversity(positions, ac, start=1, stop=1000000)

pi_w, windows, n_bases, counts = allel.windowed_diversity(positions, ac, size=100000, start=1, stop=1000000)
n_obs, minmax, mean, var, skew, kurt = st.describe(pi_w)
