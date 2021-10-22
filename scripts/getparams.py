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
import scypy.stats as st
from datetime import date
import timeadapt

def main():
  # get options for project and simulation:
  options = timeadapt.get_project_options(proj_options_file = sys.argv[1])
  sims        = range(1,options["num_of_sims"]+1)

  # print program name
  timeadapt.print_info(sys.argv[0],options["verbose"],batch=options["batch"])

  # set seed for RNG
  np.random.seed(options["seed"])
  
  # create results directory and project/batch subdirectories
  try : os.mkdir("results")
  except FileExistsError: 
    if options["verbose"] >=1 : print("results folder already exists")
  project_dir = "results/"+options["project"]  
  try : os.mkdir(project_dir)
  except FileExistsError: 
    if options["verbose"] >=1 : print("{} folder already exists".format(project_dir))
  batch_dir = project_dir+"/"+options["batch"]  
  try : os.mkdir(batch_dir)
  except FileExistsError: 
    if options["verbose"] >=1 : print("{} folder already exists".format(batch_dir))

  # read sample and genome information files
  sample_info = timeadapt.read_sample_info(sample_info_file=options["sample_file"])
  genome_info = timeadapt.get_genome_info(options["genome_file"])

  # number of generations to simulate in forward (in SLiM)
  if options["verbose"] >=10 : print("generations_forward:"+str(options["generations_forward"]))
  if options["verbose"] >=10 : print("times_of_change_forw:"+str(options["times_of_change_forw"]))

  # calculate probability distribution curves for calibrated age of ancient samples
  # (from 14C ages)
  if any(sample_info["is_ancient"]):
    age14C_NAremoved = [x for x in sample_info["age14C"] if np.isnan(x) == False]
    age14Cerror_NAremoved = [x for x in sample_info["age14Cerror"] if np.isnan(x) == False]
    age_pdf = timeadapt.get_age_pdf(age14C_NAremoved,age14Cerror_NAremoved,'shcal13')
    if options["verbose"] >=100 : print(age_pdf)
    
  # get the oldest possible age of samples and verify that falls within forward time simulation period
  if any(sample_info["is_modern"]):
    oldest_sample_age = [min(sample_info["ageBCAD"])]
    if any(sample_info["is_ancient"]):
      for i in range(0,sum(sample_info["is_ancient"])):
        if options["verbose"] >=10 : print("oldest_sample_age " + str(oldest_sample_age))
        if options["verbose"] >=10 : print("age_pdf " + str(age_pdf))
        oldest_sample_age = [min(oldest_sample_age + age_pdf[i]["ageBCAD"])]
  else:
    oldest_sample_age = [date.today().year]
    for i in range(0,sum(sample_info["is_ancient"])):
      oldest_sample_age = [min(oldest_sample_age + age_pdf[i]["ageBCAD"])]
  oldest_sample_age = abs(oldest_sample_age-sample_info["t0"])/options["gen_len_min"]
  assert oldest_sample_age+2 < options["generations_forward"],"Larger period in forward simulation necessary to include oldest sample"
  if options["verbose"] >=10 : print("oldest_sample_age " + str(oldest_sample_age))
  if options["verbose"] >=10 : print("generations_forward " + str(options["generations_forward"]))

  ###### SAMPLE FROM PRIORS ############################
  # sample generation length (generation time)
  
  gen_len = st.beta.rvs(a=1, b=1, loc=0, scale=1, size=1, random_state=None)
  


  
  


  for sim in sims:
    if options["verbose"] >=0 : print("simulation "+str(sim))

    

############################################################################################################
if __name__ == "__main__":
    main()

