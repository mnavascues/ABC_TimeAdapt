# run
# python3 src/test_recombination_map.py 7

import msprime
import pyslim
import sys

# read recombination map and load it for msprime model
file_recomb_map = open("data/recombination_map_msprime.txt", "r")
positions = []
rates = []
for line in file_recomb_map:
    p, r = line.split()
    positions.append(int(p))
    rates.append(float(r))

num_of_chr = int(sys.argv[1])

recomb_map = msprime.RecombinationMap(positions=positions[0:(num_of_chr*2)],
                                      rates = rates[0:(num_of_chr*2)])

# read tree, recapitate & add mutations
treesq = pyslim.load("results/test_recombination_map.tree")
treesq = treesq.recapitate(Ne=200, recombination_map=recomb_map, model="dtwf")

