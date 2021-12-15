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
import timeadapt

# import scripts.timeadapt as timeadapt
# options = timeadapt.get_options('results/test/project_options.ini', 'results/test/1/sim_1.ini')

def main():
  # get options for project and simulation:
  options = timeadapt.get_options(proj_options_file = sys.argv[1], sim_options_file = sys.argv[2])

  # print program name
  timeadapt.print_info(sys.argv[0],options["verbose"],batch=options["batch"],sim=options["sim"])
  
  # get recombination map:
  rate_map = msprime.RateMap(position = options["msprime_r_map"]["positions"],
                             rate = options["msprime_r_map"]["rates"]) 
  # get demography
  demography = msprime.Demography()
  demography.add_population(name="focal", initial_size=options["N"][0])
  for i in range(0,len(options["times_of_change_back"])):
    demography.add_population_parameters_change(time         = options["times_of_change_back"][i],
                                                initial_size = options["N"][i+1],
                                                growth_rate  = None,
                                                population   = "focal")
  
  # simulate with msprime
  # https://tskit.dev/msprime/docs/latest/intro.html
  msp_ts = msprime.sim_ancestry(samples            = options["N"][0],
                                demography         = demography,
                                model              = "dtwf", # because sample size = pop size
                                recombination_rate = rate_map,
                                random_seed        = options["seed_coal"])

  # make tree-sequence a SLiM-tree-sequence
  slim_ts = pyslim.annotate_defaults(msp_ts, model_type = "WF", slim_generation = 1)
  
  # save tree
  slim_ts.dump("results/"+options["project"]+"/"+options["batch"]+"/coalsim_"+options["sim"]+".trees")


############################################################################################################
############################################################################################################
if __name__ == "__main__":
    main()
