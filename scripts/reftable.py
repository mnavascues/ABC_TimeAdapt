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
import dill
import pandas as pd
import timeadapt

def main():
  # get options for project and simulation:
  options = timeadapt.get_project_options(proj_options_file = sys.argv[1])
  sims        = range(1,options["num_of_sims"]+1)

  # print program name
  timeadapt.print_info(sys.argv[0],options["verbose"])
    
  for sim in sims:
    if options["verbose"] >=0 : print("simulation "+str(sim))
    sumstats_sim = dill.load(file = open("results/"+options["project"]+"/"+options["batch"]+"/sumstats_"+str(sim)+".pkl", "rb"))
    print(pd.DataFrame.from_dict(sumstats_sim))
    
    

############################################################################################################
if __name__ == "__main__":
    main()

