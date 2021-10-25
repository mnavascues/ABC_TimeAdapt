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
  window_size = 50000

  # Ignoring some warnings when calculating summary stats with missing data,
  # typically on Tajima's D and Fst (division by zero, etc) 
  np.seterr(invalid='ignore',divide='ignore')
  
  # get options for project and simulation:
  options = timeadapt.get_options(proj_options_file = sys.argv[1], sim_options_file = sys.argv[2])

  # print program name
  timeadapt.print_info(sys.argv[0],options["verbose"],sim=options["sim"])

  # set random seed:
  np.random.seed(options["seed_mut"])

  # read sample and genome info file
  sample = timeadapt.read_sample_info(sample_info_file=options["sample_file"])
  genome_info = timeadapt.get_genome_info(options["genome_file"])

  if options["verbose"]>=100 : print(genome_info["nchr"])
  if options["verbose"]>=100 : print(genome_info["chr_ends"])
  if options["verbose"]>=100 : print(genome_info["msprime_r_map"]["rates"])
  if options["verbose"]>=100 : print(genome_info["msprime_r_map"]["positions"])

  # check sample size from sample file equal to simulated sampled size in config file
  if sum(options["ss"]) != sample["size"]:
    msg = "Number of samples from *.ini file (sum of ss=" + str(sum(options["ss"])) + \
          ") and number of samples from sample info file (sample size=" + str(sample["size"]) + \
          ") do not match"
    raise ValueError(msg)

  # set sample features in chronological order (which can be different for each simulation)
  chrono_order_coverage = [sample["coverage"][i] for i in options["chrono_order"]]
  chrono_order_is_dr = [sample["is_dr"][i] for i in options["chrono_order"]]
  chrono_order_groups = np.zeros([sample["group_levels"], sample["size"]], dtype='int')
  groups_in_level = {}
  unique_groups = {}
  num_of_pair_comparisons = 0
  number_of_groups = np.zeros(sample["group_levels"], dtype='int')
  # total_number_of_groups = 0
  for lev in range(0, sample["group_levels"]):
    if options["verbose"]>=100 : print("level "+str(lev))
    chrono_order_groups[lev] = [sample["groups"][lev][i] for i in options["chrono_order"]]
    number_of_groups[lev] = len(np.unique(sample["groups"][lev]))
    # total_number_of_groups += number_of_groups[lev]
    num_of_pair_comparisons += int((number_of_groups[lev] * (number_of_groups[lev] - 1)) / 2)
    for g in range(0, number_of_groups[lev]):
      new_group = np.where(chrono_order_groups[lev] == g)[0]
      if options["verbose"]>=100 : print("  group "+str(g)+":"+str(new_group))
      unique_groups['level' + str(lev) + 'group' + str(g)] = True
      for key, value in groups_in_level.items():
        if options["verbose"]>=100 : print("key: "+str(key)+" ; value: "+str(value))
        if np.size(value)==np.size(new_group):
          if (value==new_group).all() : unique_groups['level' + str(lev) + 'group' + str(g)] = False
      groups_in_level['level' + str(lev) + 'group' + str(g)] = new_group

  if options["verbose"]>=100 : print(unique_groups)

  # read tree sequence from SLiM output file:
  treesq = pyslim.load("results/"+options["project"]+"/"+options["batch"]+"/forwsim_"+options["sim"]+".trees")

  # Simulate neutral mutation over the tree sequence
  mut_treesq = msprime.sim_mutations(treesq,
                                     rate = options["mu"],
                                     random_seed = np.random.randint(1, 2**32-1))
  if options["verbose"]>=100 : print("Number of mutations " + str(mut_treesq.num_mutations))
  if options["verbose"]>=100 : print("Number of sites " + str(mut_treesq.num_sites))


  #### CALCULATE SUMMARY STATISTICS
  ####------------------------------
  
  ref_table_sumstats = {}
  
  if mut_treesq.num_sites == 0:
    print("No mutations")
    # TODO: Create empty sumstats
  else:
    geno_data, positions = timeadapt.sequencing(ts = mut_treesq,
                                                ssize = sample["size"],
                                                ttr = options["ttratio"],
                                                seq_error = options["seq_error"],
                                                dr = chrono_order_is_dr,
                                                cov = chrono_order_coverage)
    if options["verbose"]>=100 : print(genome_info["msprime_r_map"]["positions"])
    # Calculate summary statistics from total sample
    timeadapt.single_sample_sumstats(ga = geno_data,
                                     pos = positions,
                                     nchr = genome_info["nchr"],
                                     chr_ends = genome_info["chr_ends"],
                                     w_size = window_size,
                                     sumstats = ref_table_sumstats,
                                     sep = "")
    # Calculate summary statistics from single group
    for lev in range(0, sample["group_levels"]):
      for g in range(0, number_of_groups[lev]):
        if unique_groups['level'+str(lev)+'group'+str(g)] is True:
          timeadapt.single_sample_sumstats(ga = geno_data[:, groups_in_level['level'+str(lev)+'group'+str(g)]],
                                           pos = positions,
                                           nchr = genome_info["nchr"],
                                           chr_ends = genome_info["chr_ends"],
                                           w_size = window_size,
                                           sumstats = ref_table_sumstats,
                                           name = 'l'+str(lev)+'g'+str(g))

    # Calculate summary statistics from pair of groups
    for lev in range(0, sample["group_levels"]):
      if options["verbose"]>=100 : print("Level: " + str(lev))
      for g1 in range(0, number_of_groups[lev]):
        for g2 in range(g1+1, number_of_groups[lev]):
          pair = [groups_in_level['level'+str(lev)+'group'+str(g1)],
                  groups_in_level['level'+str(lev)+'group'+str(g2)]]
          if options["verbose"]>=100 : print(" Pair:")
          if options["verbose"]>=100 : print("  1st group: " + str(pair[0]))
          if options["verbose"]>=100 : print("  2nd group: " + str(pair[1]))
          timeadapt.two_samples_sumstats(ga = geno_data,
                                         pair_of_groups = pair,
                                         pos = positions,
                                         nchr = genome_info["nchr"],
                                         chr_ends = genome_info["chr_ends"],
                                         w_size =  window_size,
                                         sumstats = ref_table_sumstats,
                                         name = 'l'+str(lev)+'p'+str(g1)+str(g2) )

  # Print summary statistics to screen
  if options["verbose"]>=100:
    for key,value in ref_table_sumstats.items():
      print(key, ':', value)
  # save binary with row of summary statistics for the reference table:
  dill.dump(ref_table_sumstats, file = open("results/" + options["project"] + "/" + options["batch"] + "/sumstats_" + options["sim"] + ".pkl", "wb"))

############################################################################################################
if __name__ == "__main__":
    main()
