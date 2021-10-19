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
import os
import numpy as np
import timeadapt

def main():
  # get options for project and simulation:
  options = timeadapt.get_project_options(proj_options_file = sys.argv[1])
  sims        = range(1,options["num_of_sims"]+1)

  # print program name
  timeadapt.print_info(sys.argv[0],options["verbose"])

  # set seed for RNG
  np.random.seed(options["seed"])
  
  # create results directory and project/batch subdirectories
  try : os.mkdir("results")
  except FileExistsError: 
    if options["verbose"] >=0 : print("results folder already exists")
  project_dir = "results/"+options["project"]  
  try : os.mkdir(project_dir)
  except FileExistsError: 
    if options["verbose"] >=0 : print("{} folder already exists".format(project_dir))
  batch_dir = project_dir+"/"+options["batch"]  
  try : os.mkdir(batch_dir)
  except FileExistsError: 
    if options["verbose"] >=0 : print("{} folder already exists".format(batch_dir))

  # read sample and genome information files
  sample_id, coverage, is_ancient, is_modern, is_dr, total_ancient, \
           sample_size, group_levels, groups = timeadapt.read_sample_info(options["sample_file"])
  nchr, chr_ends, rates, positions = timeadapt.get_recombination_map(options["genome_file"])


  for sim in sims:
    if options["verbose"] >=0 : print("simulation "+str(sim))

    

############################################################################################################
if __name__ == "__main__":
    main()

