import msprime
import pyslim

# read recombination map and load it for msprime model
file_recomb_map = open("data/recombination_map_msprime.txt", "r")
positions = []
rates = []
for line in file_recomb_map:
    p, r = line.split()
    positions.append(int(p))
    rates.append(float(r))

num_of_chr = 7

recomb_map = msprime.RecombinationMap(positions=positions[0:(num_of_chr*2)],
                                      rates = rates[0:(num_of_chr*2)],
                                      num_loci = positions[num_of_chr*2-1])

# read tree, recapitate & add mutations
treesq = pyslim.load("results/slim_test.tree")
treesq = treesq.recapitate(Ne=N, recombination_map=recomb_map, model="dtwf")
