import numpy as np
import msprime

samples=[ (0,  0),(0,  0),
          (0, 10),(0, 10),
          (0,100),(0,100) ]

tree_seq = msprime.simulate(samples=samples,
                            Ne=1000,
                            length=5e4,
                            recombination_rate=2e-8,
                            mutation_rate=1e-8,
                            random_seed=30)


with open("results/neutral_genotypes.vcf", "w") as vcffile:
  tree_seq.write_vcf(vcffile)
