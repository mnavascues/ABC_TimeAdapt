import sys
import msprime
import pyslim
import numpy as np
import timeadapt

## how to run:
#  python scripts/coalsim.py tests/input/config_project.ini tests/input/sim_1.ini

def main():
  # get options for project and simulation:
  project, batch, sim, genome_file, sample_file, N, seed_coal, seed_mut = \
           timeadapt.get_options(proj_options_file = sys.argv[1], sim_options_file = sys.argv[2])
  
  # set random seed:
  np.random.seed(seed_coal)
  # get recombination map:
  num_of_chr, chrom_rates, chrom_ends = timeadapt.get_genome_map(gf = genome_file)
  # make_recombination map:
  positions, rates = timeadapt.make_rec_map(nchr = num_of_chr,
                                            c_rates = chrom_rates,
                                            c_ends = chrom_ends)
  rate_map = msprime.RateMap(position = positions, rate = rates) 
  # simulate with msprime
  # https://tskit.dev/msprime/docs/latest/intro.html
  msp_ts = msprime.sim_ancestry(samples            = N[0],
                                population_size    = N[0],
                                model              = "dtwf",
                                recombination_rate = rate_map,
                                random_seed        = np.random.randint(1, 2 ^ 32 - 1))
  # make tree sequence a SLiM tree sequence
  slim_ts = pyslim.annotate_defaults(msp_ts, model_type="WF", slim_generation=1)
  
  # save tree
  slim_ts.dump("results/"+project+"/"+batch+"/coalsim_"+sim+".trees")


############################################################################################################
############################################################################################################
if __name__ == "__main__":
    main()
