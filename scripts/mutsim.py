#    TimeAdapt: joint inference of demography and selection
#    Copyright (C) 2021  Miguel de Navascués, Uppsala universitet, INRAE
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <https://www.gnu.org/licenses/>.
#
import sys
import msprime
import pyslim
import numpy as np
import allel
import timeadapt

def main():
  # get options for project and simulation:
  project, batch, sim, genome_file, sample_file, ss, chrono_order, N, mu, ttratio, seq_error, seed_coal, seed_mut = \
           timeadapt.get_options(proj_options_file = sys.argv[1], sim_options_file = sys.argv[2])

  # set random seed:
  np.random.seed(seed_mut)

  # read sample info file
  sample_id, coverage, is_ancient, is_modern, is_dr, total_ancient, \
    sample_size, group_levels, \
    groups = timeadapt.read_sample_info(sample_info_file=sample_file)

  _, _, map_positions = timeadapt.get_recombination_map(gf = genome_file)


  # check sample size from sample file equal to simulated sampled size in config file
  if sum(ss) != sample_size:
    msg = "Number of samples from *.ini file (sum of ss=" + str(sum(ss)) + \
          ") and number of samples from sample info file (sample_size=" + str(sample_size) + \
          ") do not match"
    raise ValueError(msg)

  # set sample features in chronological order (which can be different for each simulation)
  chrono_order_coverage = [coverage[i] for i in chrono_order]
  chrono_order_is_dr = [is_dr[i] for i in chrono_order]
  chrono_order_groups = np.zeros([group_levels, sample_size], dtype='int')
  groups_in_level = {}
  num_of_pair_comparisons = 0
  number_of_groups = np.zeros(group_levels, dtype='int')
  total_number_of_groups = 0
  for lev in range(0, group_levels):
    chrono_order_groups[lev] = [groups[lev][i] for i in chrono_order]
    number_of_groups[lev] = len(np.unique(groups[lev]))
    total_number_of_groups += number_of_groups[lev]
    num_of_pair_comparisons += int((number_of_groups[lev] * (number_of_groups[lev] - 1)) / 2)
    for g in range(0, number_of_groups[lev]):
      groups_in_level['level' + str(lev) + 'group' + str(g)] = np.where(chrono_order_groups[lev] == g)[0]


  # read tree sequence from SLiM output file:
  treesq = pyslim.load("results/"+project+"/"+batch+"/forwsim_"+sim+".trees")

  mut_treesq = msprime.sim_mutations(treesq,
                                     rate = mu,
                                     random_seed = np.random.randint(1, 2 ^ 32 - 1))

  print("Number of mutations " + str(mut_treesq.num_mutations))
  print("Number of sites " + str(mut_treesq.num_sites))
  # TODO : what do you do when num of mutations > num of sites ?????
  if mut_treesq.num_sites == 0:
    print("No mutations")
    # TODO: Create empty sumstats
  else:
    geno_data, positions = timeadapt.sequencing(ts=mut_treesq,
                                                ssize=sample_size,
                                                ttr=ttratio,
                                                seq_error=seq_error,
                                                dr=chrono_order_is_dr,
                                                cov=chrono_order_coverage)
                                                
    print("length geno_data "+str(len(geno_data)))
    print("length positions "+str(len(positions)))
    # print("geno_data"+str(geno_data))
    print("map_positions[-1] "+str(map_positions[-1]))
    segsites, pi, min_pi, max_pi, mean_pi, variance_pi, skewness_pi,\
         kurtosis_pi, min_W_theta, max_W_theta, mean_W_theta,\
         variance_W_theta, skewness_W_theta, kurtosis_W_theta, Taj_D,\
         min_Taj_D, max_Taj_D, mean_Taj_D, variance_Taj_D, skewness_Taj_D,\
         kurtosis_Taj_D, roh_distribution\
         = timeadapt.single_sample_sumstats(geno_data, positions, map_positions[-1], 50000)
    print("segsites: "+ str(segsites))
    print("pi: "+ str(pi))
    print("Taj_D: "+ str(Taj_D))
    print("roh_distribution: "+ str(roh_distribution))
    print(np.all(positions[1:] > positions[:-1]))
 
    print("-------------------------------------------------------------------")
    #print("Genotype data matrix:")
    #print(geno_data)
    #allele_counts = geno_data.count_alleles()
    #print("Allele count=\n"+str(allele_counts))
    # TODO : make a function to calculate single sample summary stats, call it here and below
    #pi_per_lg = allel.sequence_diversity(positions, allele_counts,
    #                                     start=1,
    #                                     stop=40000000)
    #he_per_site = allel.mean_pairwise_difference(allele_counts)
    #pi_per_window, windows, n_bases, \
    #  n_sites = allel.windowed_diversity(positions, allele_counts,
    #                                     size=50000,
    #                                     start=1,
    #                                     stop=40000000)
    #allele_counts_per_group = {}
    #for lev in range(0, group_levels):
    #  for g in range(0, number_of_groups[lev]):
    #    allele_counts_per_group['level' + str(lev) + 'group' + str(g)] = \
    #      geno_data.count_alleles(subpop=groups_in_level['level' + str(lev) + 'group' + str(g)])
        # TODO : call here function for single sample sumstats

    #for lev in range(0, group_levels):
    #  for g in range(0, number_of_groups[lev]):
    #    for h in range(g + 1, number_of_groups[lev]):
    #      gr1 = allele_counts_per_group['level' + str(lev) + 'group' + str(g)]
    #      gr2 = allele_counts_per_group['level' + str(lev) + 'group' + str(h)]
    #      pairwise_diff = allel.mean_pairwise_difference_between(gr1, gr2)
    #      print("Level: " + str(lev) + ". Groups: " + str(g) + " " + str(h) +
    #            ". Pairwise difference: " + str(pairwise_diff))
                                              

  outfile = open("results/" + project + "/" + batch + "/sumstats_" + sim + ".txt",
                   "w")



############################################################################################################
if __name__ == "__main__":
    main()
