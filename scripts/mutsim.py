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
import dill
import timeadapt

def main():
  verbose=10
  
  # get options for project and simulation:
  project, batch, sim, genome_file, sample_file, ss, chrono_order, N, mu, ttratio, seq_error, seed_coal, seed_mut = \
           timeadapt.get_options(proj_options_file = sys.argv[1], sim_options_file = sys.argv[2])

  # set random seed:
  np.random.seed(seed_mut)

  # read sample info file
  sample_id, coverage, is_ancient, is_modern, is_dr, total_ancient, \
    sample_size, group_levels, \
    groups = timeadapt.read_sample_info(sample_info_file=sample_file)

  nchr, chr_ends, map_rates, map_positions = timeadapt.get_recombination_map(gf = genome_file)

  if verbose>=100 : print(nchr)
  if verbose>=100 : print(chr_ends)
  if verbose>=100 : print(map_rates)
  if verbose>=100 : print(map_positions)

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
  unique_groups = {}
  num_of_pair_comparisons = 0
  number_of_groups = np.zeros(group_levels, dtype='int')
  total_number_of_groups = 0
  for lev in range(0, group_levels):
    if verbose>=100 : print("level "+str(lev))
    chrono_order_groups[lev] = [groups[lev][i] for i in chrono_order]
    number_of_groups[lev] = len(np.unique(groups[lev]))
    # total_number_of_groups += number_of_groups[lev]
    num_of_pair_comparisons += int((number_of_groups[lev] * (number_of_groups[lev] - 1)) / 2)
    for g in range(0, number_of_groups[lev]):
      new_group = np.where(chrono_order_groups[lev] == g)[0]
      if verbose>=100 : print("  group "+str(g)+":"+str(new_group))
      unique_groups['level' + str(lev) + 'group' + str(g)] = True
      for key, value in groups_in_level.items():
        if verbose>=100 : print("key: "+str(key)+" ; value: "+str(value))
        if np.size(value)==np.size(new_group):
          if (value==new_group).all() : unique_groups['level' + str(lev) + 'group' + str(g)] = False
      groups_in_level['level' + str(lev) + 'group' + str(g)] = new_group

  if verbose>=100 : print(unique_groups)

  # read tree sequence from SLiM output file:
  treesq = pyslim.load("results/"+project+"/"+batch+"/forwsim_"+sim+".trees")

  # Simulate neutral mutation over the tree sequence
  mut_treesq = msprime.sim_mutations(treesq,
                                     rate = mu,
                                     random_seed = np.random.randint(1, 2 ^ 32 - 1))
  if verbose>=100 : print("Number of mutations " + str(mut_treesq.num_mutations))
  if verbose>=100 : print("Number of sites " + str(mut_treesq.num_sites))


  #### CALCULATE SUMMARY STATISTICS
  ####------------------------------
  
  ref_table_sumstats = {}
  
  if mut_treesq.num_sites == 0:
    print("No mutations")
    # TODO: Create empty sumstats
  else:
    geno_data, positions = timeadapt.sequencing(ts = mut_treesq,
                                                ssize = sample_size,
                                                ttr = ttratio,
                                                seq_error = seq_error,
                                                dr = chrono_order_is_dr,
                                                cov = chrono_order_coverage)
    if verbose>=100 : print(map_positions)
    # Calculate summary statistics from total sample
    timeadapt.single_sample_sumstats(ga = geno_data,
                                     pos = positions,
                                     nchr = 1,
                                     chr_ends = [map_positions[-1]],
                                     w_size = 50000,
                                     sumstats = ref_table_sumstats,
                                     sep = "")
    # Calculate summary statistics from single group
    for lev in range(0, group_levels):
      for g in range(0, number_of_groups[lev]):
        if unique_groups['level'+str(lev)+'group'+str(g)] is True:
          timeadapt.single_sample_sumstats(ga = geno_data[:, groups_in_level['level'+str(lev)+'group'+str(g)]],
                                           pos = positions,
                                           nchr = 1,
                                           chr_ends = [map_positions[-1]],
                                           w_size = 50000,
                                           sumstats = ref_table_sumstats,
                                           name = 'l'+str(lev)+'g'+str(g))

    # Calculate summary statistics from pair of groups
    for lev in range(0, group_levels):
      for g1 in range(0, number_of_groups[lev]):
        for g2 in range(g1+1, number_of_groups[lev]):
          if verbose>=10 : print("PAIR:")
          if verbose>=10 : print("  1st group: "+str(groups_in_level['level'+str(lev)+'group'+str(g1)]))
          if verbose>=10 : print("  2nd group: "+str(groups_in_level['level'+str(lev)+'group'+str(g2)]))
 
                                              
  # save binary with row of summary statistics for the reference table:
  dill.dump(ref_table_sumstats, file = open("results/" + project + "/" + batch + "/sumstats_" + sim + ".pkl", "wb"))
  # THIS IS HOW YOU LOAD THE DATA AGAIN: ref_table_sumstats = dill.load(file = open("results/" + project + "/" + batch + "/sumstats_" + sim + ".pkl", "rb"))

############################################################################################################
if __name__ == "__main__":
    main()
