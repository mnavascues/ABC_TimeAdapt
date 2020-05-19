import allel
import msprime
import numpy as np
import pyslim
import scipy.stats as st

# i     = int(sys.argv[1])
i = 1
# N     = int(sys.argv[2])
N = 200
# mu    = float(sys.argv[3])
mu = 5e-9
# seed  = int(sys.argv[4])
seed = 1234
# n0    = int(sys.argv[5])
n0 = 10
# na    = int(sys.argv[6])
na = 8
# ts = ()
# for x in range(0, na):
#  ts=ts+( int(sys.argv[7+x]), )
ts = (105, 139, 166, 259, 292, 337, 356, 412)

np.random.seed(seed)

recomb_map = msprime.RecombinationMap.read_hapmap("data/RecombinationMap2010/genetic_map_GRCh37_chr22.txt")
demogr_event = [msprime.PopulationParametersChange(time=1000, initial_size=300, population_id=0)]

# read tree, recapitate & add mutations
treesq = pyslim.load("results/slim" + str(i) + ".tree")
treesq = treesq.recapitate(Ne=N, recombination_map=recomb_map, demographic_events=demogr_event)

# sample_ind = treesq.individuals_alive_at(0)
sample_ind = np.random.choice(treesq.individuals_alive_at(0), n0, replace=False)
for x in range(0, na):
    sample_ind = np.concatenate([sample_ind, treesq.individuals_alive_at(ts[x])])

keep_nodes = []
for i in sample_ind:
    keep_nodes.extend(treesq.individual(i).nodes)

treesq = treesq.simplify(keep_nodes)

sample_ind = np.random.choice(treesq.individuals_alive_at(0), n0, replace=False)
for x in range(0, na):
    sample_ind = np.concatenate([sample_ind, treesq.individuals_alive_at(ts[x])])

# check age of sample:  treesq.individual(0).time


treesq = pyslim.SlimTreeSequence(msprime.mutate(treesq, rate=mu))

ploidy = 2
n_samples = sample_ind.shape[0]
n_loci = treesq.num_mutations

geno_data = allel.GenotypeArray(np.full((n_loci, n_samples, ploidy), -1), dtype='i1')

coverage = (10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 1, 2, 3, 4, 0.1, 0.2, 0.3, 0.4)
sequencing_error = 0.005

# check age of sample:  treesq.individual(0).time
positions = []
locus = 0
for variant in treesq.variants():
    positions.append(round(variant.position))
    var_genotypes = variant.genotypes
    print(var_genotypes)
    num_reads = st.poisson.rvs(mu=coverage, size=n_samples)
    print(num_reads)
    for i in range(0, 2 * n_samples, 2):
        if num_reads[int(i / 2)] >= 1:
            current_genotype = var_genotypes[i:(i + 2)]
            derived_count = sum(current_genotype)
            p_derived = derived_count / 2. * (1 - sequencing_error) + (1 - derived_count / 2.) * sequencing_error
            derived_reads = st.binom.rvs(num_reads[int(i / 2)], p_derived)
            ancestral_reads = num_reads[int(i / 2)] - derived_reads
            if num_reads[int(i / 2)] > 10:
                if derived_reads==0:
                    geno_data[locus, int(i / 2)] = [0, 0]
                elif ancestral_reads==0:
                    geno_data[locus, int(i / 2)] = [1, 1]
                else:
                    ratio_of_scores = derived_reads/ancestral_reads
                    if (ratio_of_scores >= 1/3) & (ratio_of_scores <= 3):
                        geno_data[locus, int(i / 2)] = [0, 1]
                    elif derived_reads > ancestral_reads:
                        geno_data[locus, int(i / 2)] = [1, 1]
                    else:
                        geno_data[locus, int(i / 2)] = [0, 0]
            else:
                random_allele = st.binom.rvs(1, derived_reads / num_reads[int(i / 2)])
                geno_data[locus, int(i / 2)] = np.full(2, random_allele)
    locus = locus + 1


ac = geno_data.count_alleles()
pi = allel.sequence_diversity(positions, ac, start=1, stop=1000000)

pi_w, windows, n_bases, counts = allel.windowed_diversity(positions, ac, size=100000, start=1, stop=1000000)
n_obs, minmax, mean, var, skew, kurt = st.describe(pi_w)

