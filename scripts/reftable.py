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
  timeadapt.print_info(sys.argv[0],options["verbose"],batch=options["batch"])
  
  # read latent variables file
  latent_variables_file = "results/"+options["project"]+"/"+options["batch"]+"/latent_variables.txt"
  ref_table_latent_variables = pd.read_table(filepath_or_buffer=latent_variables_file, index_col="sim", sep=r'\s+')
  if options["verbose"] >=10 : print(ref_table_latent_variables)    
  ref_table_latent_variables = ref_table_latent_variables.sort_index(axis=0,ascending=True)
  if options["verbose"] >=0 : print(ref_table_latent_variables)    
 
  
  ref_table_sumstats  = pd.DataFrame()

  for sim in sims:
    if options["verbose"] >=0 : print("simulation "+str(sim))
    sumstats_file = "results/"+options["project"]+"/"+options["batch"]+"/sumstats_"+str(sim)+".pkl"
    sumstats_sim = dill.load(file = open(sumstats_file, "rb"))
    if options["verbose"] >=100 : print(sumstats_sim)
    if options["verbose"] >=100 : print(pd.DataFrame(sumstats_sim,index=[sim]))
    ref_table_sumstats = ref_table_sumstats.append(pd.DataFrame(sumstats_sim,index=[sim]))
  
  dill.dump(ref_table_sumstats, file = open("results/" + options["project"] + "/" + options["batch"] + "/ref_table_sumstats.pkl", "wb"))
  dill.dump(ref_table_latent_variables, file = open("results/" + options["project"] + "/" + options["batch"] + "/ref_table_latent_variables.pkl", "wb"))
  
  if options["verbose"] >=0 : print(ref_table_sumstats)    
    

############################################################################################################
if __name__ == "__main__":
    main()

