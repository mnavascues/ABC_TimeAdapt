import sys
import configparser
import msprime
import pyslim
import math


def main():
    
  # Get project and sim parameters from ini files
  proj_options_file = sys.argv[1]
  sim_options_file = sys.argv[2]
  proj_options = configparser.ConfigParser()
  proj_options.read(proj_options_file)
  sim_options = configparser.ConfigParser()
  sim_options.read(sim_options_file)

  project     = proj_options.get('Settings','project')
  print(project)
  batch       = proj_options.get('Settings','batch')
  print(batch)
  N           = int(sim_options.get('Demography','N'))
  print(N)




  num_of_chr  = 6

  chrom_ends = [249218992,492309989,690184517,881213162,1061928972,
              1232980288,1392107062,1538408218,1679518197,1815021890,
              1949967412,2083746208,2198853204,2306141581,2408647689,
              2498811277,2579856546,2657871726,2716968886,2779918331,
              2828018486,2879248291]
  chrom_rates = [1.14856e-08,1.10543e-08,1.12796e-08,1.12312e-08,1.12809e-08,
               1.12229e-08,1.17646e-08,1.14785e-08,1.17807e-08,1.33651e-08,
               1.17193e-08,1.30502e-08,1.09149e-08,1.11973e-08,1.38358e-08,
               1.48346e-08,1.58249e-08,1.5076e-08,1.82201e-08,1.71783e-08,
               1.30452e-08,1.4445e-08]             
             
  positions = [0]
  rates = []
  for i in range(num_of_chr-1):
    print(i)
    positions.append(chrom_ends[i])
    positions.append(chrom_ends[i]+1)
    rates.append(chrom_rates[i])
    rates.append(math.log(2))

  positions.append(chrom_ends[num_of_chr-1])
  rates.append(chrom_rates[num_of_chr-1])

  print(positions)
  print(rates)

  rate_map = msprime.RateMap(
    position = positions,
    rate     = rates
  ) 

  msp_ts = msprime.sim_ancestry(
    samples            = N,
    population_size    = N,
    model              = "dtwf",
    recombination_rate = rate_map,
    random_seed        = 1234
  )
  print(msp_ts)

  slim_ts = pyslim.annotate_defaults(msp_ts, model_type="WF", slim_generation=1)
  print(slim_ts)

  slim_ts.dump("results/"+project+"/"+batch+"/coalsim_1.trees")
  #msp_ts.draw_svg("results/"+project+"/"+batch+"/coalsim_1.svg")












############################################################################################################
############################################################################################################
if __name__ == "__main__":
    main()
