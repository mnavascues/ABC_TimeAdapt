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
import configparser  # for writing ini files
import numpy as np
import pandas as pd
import dill
import scipy.stats as st
from datetime import date
import timeadapt

def main():
  # get options for project and simulation:
  options = timeadapt.get_project_options(proj_options_file = sys.argv[1])
  sims  = range(0,options["num_of_sims"])
  batch = sys.argv[2]
  
  # print program name
  timeadapt.print_info(sys.argv[0],options["verbose"],batch=batch)

  # set seed for RNG
  np.random.seed(options["seed"])
  
  # create results directory and project/batch subdirectories
  try : os.mkdir("results")
  except FileExistsError : 
    if options["verbose"] >=1 : print("results folder already exists")
  project_dir = "results/"+options["project"]  
  try : os.mkdir(project_dir)
  except FileExistsError : 
    if options["verbose"] >=1 : print("{} folder already exists".format(project_dir))
  batch_dir = project_dir+"/"+str(batch)
  try : os.mkdir(batch_dir)
  except FileExistsError : 
    if options["verbose"] >=1 : print("{} folder already exists".format(batch_dir))

  # read sample and genome information files
  sample = timeadapt.read_sample_info(sample_info_file=options["sample_file"])
  genome_info = timeadapt.get_genome_info(options["genome_file"])

  # create Eidos file with SLiM recombination map
  with open(project_dir+'/recombination_map.eidos', 'w') as rec_map_file:
    rec_map_file.write('defineConstant("ends",c(' + str(','.join(map(str, genome_info["slim_r_map"]["positions"] ))) + '));\n')
    rec_map_file.write('defineConstant("rates",c(' + str(','.join(map(str, genome_info["slim_r_map"]["rates"] ))) + '));\n')

  # number of generations to simulate in forward (in SLiM)
  if options["verbose"] >=10 : print("generations_forward:"+str(options["generations_forward"]))
  if options["verbose"] >=10 : print("times_of_change_forw:"+str(options["times_of_change_forw"]))

  # calculate probability distribution curves for calibrated age of ancient samples
  # (from 14C ages)
  if any(sample["is_ancient"]):
    age14C_NAremoved = [x for x in sample["age14C"] if np.isnan(x) == False]
    age14Cerror_NAremoved = [x for x in sample["age14Cerror"] if np.isnan(x) == False]
    age_pdf = timeadapt.get_age_pdf(age14C_NAremoved,age14Cerror_NAremoved,'shcal13')
    if options["verbose"] >=100 : print(age_pdf)
  else:
    age_pdf=None
    
  # get the oldest possible age of samples and verify that falls within forward time simulation period
  if any(sample["is_modern"]):
    oldest_sample_age = [min(sample["ageBCAD"])]
    if any(sample["is_ancient"]):
      for i in range(0,sum(sample["is_ancient"])):
        if options["verbose"] >=10 : print("oldest_sample_age " + str(oldest_sample_age))
        if options["verbose"] >=10 : print("age_pdf " + str(age_pdf))
        oldest_sample_age = [min(oldest_sample_age + age_pdf[i]["ageBCAD"])]
  else:
    oldest_sample_age = [date.today().year]
    for i in range(0,sum(sample["is_ancient"])):
      oldest_sample_age = [min(oldest_sample_age + age_pdf[i]["ageBCAD"])]
  oldest_sample_age = abs(oldest_sample_age-sample["t0"])/options["gen_len_min"]
  assert oldest_sample_age+2 < options["generations_forward"],"Larger period in forward simulation necessary to include oldest sample"
  if options["verbose"] >=10 : print("oldest_sample_age " + str(oldest_sample_age))
  if options["verbose"] >=10 : print("generations_forward " + str(options["generations_forward"]))

  total_number_of_periods = options["periods_forward"] + options["periods_coalescence"]

  # create latent variable file with headers
  latent_variables_file = batch_dir+"/latent_variables.txt"
  header = "sim"
  for i in range(0,options["periods_forward"]):
    header = header + " Ne_" + str(i)
  header = header+'\n'
  with open(batch_dir+"/latent_variables.txt", 'w') as latent_variables_file:
    latent_variables_file.write(header)

  ###### SAMPLE FROM PRIORS ############################
  # sample generation length (generation time)
  gen_len = st.beta.rvs(a=options["gen_len_sh1"],
                        b=options["gen_len_sh2"],
                        loc=options["gen_len_min"],
                        scale=options["gen_len_max"]-options["gen_len_min"],
                        size=options["num_of_sims"])
  if options["verbose"] >=10 : print("gen_len " + str(gen_len))
  # sample mutation rate
  u = st.lognorm.rvs(s=options["mut_rate_sd"],
                     scale=options["mut_rate_mean"],
                     size=options["num_of_sims"])
  # sample ages and demography for each simulation
  N = [None]*options["num_of_sims"]
  for sim in sims:
    if options["verbose"] >=0 : print("simulation "+str(sim+1))
    # sample seeds for msprime and SLiM
    seed_coal = np.random.randint(1, 2**32-1)
    seed_forw = np.random.randint(1, 2**32-1)
    seed_mut  = np.random.randint(1, 2**32-1)
    if options["verbose"] >=10 : print("seeds: "+str(seed_coal)+" -- "+str(seed_forw)+" -- "+str(seed_mut))
    # sample ages in generations
    sample_ages = timeadapt.get_sample_ages(sample,age_pdf,gen_len[sim])
    sampling_times_slim = sorted([options["generations_forward"]-age for age in set(sample_ages)])
    sampling_times_msprime = sorted(set(sample_ages))
    if options["verbose"] >=10 : print("sample_ages: "+str(sample_ages))
    # get chronological order of samples and sample size per generation
    chrono_order = sorted(range(len(sample_ages)), key=lambda k: sample_ages[k])
    if options["verbose"] >=10 : print("chrono_order: "+str(chrono_order))
    sample_size_per_gen = [sample_ages.count(age) for age in sorted(set(sample_ages))]
    if options["verbose"] >=10 : print("sample_size_per_gen: "+str(sample_size_per_gen))
    # sample effective population size
    N[sim] = timeadapt.sample_param_trajectory(times=total_number_of_periods,
                                               minimum=options["pop_size_min"],
                                               maximum=options["pop_size_max"]).astype(int)
    if options["verbose"] >=10 : print("N at present simulation: "+str(N[sim]))
    
    backward_demography=np.flip(N[sim][0:options["periods_coalescence"]])
    
    sim_config = configparser.ConfigParser()
    sim_config['Simulation'] = {'sim'          : sim+1,
                                'batch'        : batch}
    sim_config['Sample']     = {'ss'           : ' '.join(map(str, sample_size_per_gen)),
                                'msprime_ts'   : ' '.join(map(str, sampling_times_msprime)),
                                'chrono_order' : ' '.join(map(str, chrono_order))}
    sim_config['Demography'] = {'N'            : ' '.join(map(str, backward_demography))}
    sim_config['Genome']     = {'mu'           : u[sim],
                                'ttratio'      : 2.0,
                                'seq_error'    : 0.005}
    sim_config['Seeds']      = {'seed_coal'    : seed_coal,
                                'seed_mut'     : seed_mut}
    with open(batch_dir+'/sim_'+str(sim+1)+'.ini', 'w') as msprime_config_file:
      sim_config.write(msprime_config_file)

    forward_demography = N[sim][total_number_of_periods-options["periods_forward"]:total_number_of_periods]

    with open(batch_dir+'/sim_'+str(sim+1)+'.eidos', 'w') as slim_config_file:
      slim_config_file.write('setSeed(' + str(seed_forw) + ');\n')
      slim_config_file.write('defineConstant("L",' + str(genome_info["L"]) + ');\n')
      slim_config_file.write('defineConstant("N",c(' + str(','.join(map(str, forward_demography)) ) + '));\n')
      slim_config_file.write('defineConstant("tc",c(' + str(','.join(map(str, options["times_of_change_forw"])) ) + '));\n')
      slim_config_file.write('defineConstant("ts",c(' + str(','.join(map(str, sampling_times_slim ))) + '));\n')
      slim_config_file.write('defineConstant("ss",c(' + str(','.join(map(str, reversed(sample_size_per_gen) ))) + '));\n')
      slim_config_file.write('defineConstant("np",' + str(options["periods_forward"]) + ');\n')
      slim_config_file.write('defineConstant("na",' + str(len(sampling_times_slim)-1) + ');\n')
      slim_config_file.write('defineConstant("i",' + str(sim+1) + ');\n')
      slim_config_file.write('defineConstant("batch",' + str(batch) + ');\n')
      slim_config_file.write('defineConstant("project","' + str(options["project"]) + '");\n')
      slim_config_file.write('source("'+project_dir+'/recombination_map.eidos");\n')

  # export reference table: parameters
  ref_table_params = pd.DataFrame(data    = N,
                                  index   = range(1,options["num_of_sims"]+1),
                                  columns = ["N"+str(x) for x in range(0,total_number_of_periods)])
  ref_table_params = ref_table_params.assign(u=u)
  ref_table_params = ref_table_params.assign(generation_length=gen_len)
  if options["verbose"] >=0 : print(ref_table_params)
  dill.dump(ref_table_params, file = open("results/" + options["project"] + "/" + str(batch) + "/ref_table_params.pkl", "wb"))

############################################################################################################
if __name__ == "__main__":
    main()

