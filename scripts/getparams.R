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

library(extraDistr, quietly=TRUE)
library(ini, quietly = TRUE)
source("scripts/timeadapt.R")

# read command line arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0){
  stop("Two positional arguments necessary in the command (RDS options file and batch number)")
  quit(save="no")
} else {
  options_file = args[1]             # options_file = "results/test/project_options.RDS"
  batch        = as.integer(args[2]) # batch = 1
}

opts = readRDS(options_file)
Settings = opts$Settings
Model    = opts$Model
Priors   = opts$Priors
Sample   = opts$Sample
Genome   = opts$Genome

total_periods = Model$periods_forward+Model$periods_coalescence

# print script info to screen
print_info("getparams.R", Settings$verbose, project = Settings$project, batch = batch)

# set seed for random number generator and get seeds for msprime and SLiM
set.seed(Settings$seed + as.numeric(batch))
seed_coal = as.integer(round(runif(Settings$num_of_sims, 0, 2^31)))
seed_forw = as.integer(round(runif(Settings$num_of_sims, 0, 2^31)))
seed_mut  = as.integer(round(runif(Settings$num_of_sims, 0, 2^31)))
if (Settings$verbose >= 10) write(paste("Seed (msprime-coalescence)", seq_len(Settings$num_of_sims), ":", seed_coal), stdout())
if (Settings$verbose >= 10) write(paste("Seed (SLiM)", seq_len(Settings$num_of_sims), ":", seed_forw), stdout())
if (Settings$verbose >= 10) write(paste("Seed (msprime-mutation)", seq_len(Settings$num_of_sims), ":", seed_mut), stdout())

# create results directory
batch_dir <- paste(Settings$project_dir, batch, sep="/")
dir.create(batch_dir, showWarnings = FALSE)

#### SAMPLE FROM PRIORS #### 

# Generation length = generation interval in years
gen_len = rnsbeta(Settings$num_of_sims,
                  Priors$gen_len_sh1,
                  Priors$gen_len_sh2,
                  Priors$gen_len_min,
                  Priors$gen_len_max)
if (Settings$verbose >= 10) write(paste("Generation length", seq_len(Settings$num_of_sims), ":", gen_len), stdout())
# mutation rate
u = rlnorm(Settings$num_of_sims,
           log(Priors$mut_rate_mean),
           Priors$mut_rate_sd)
if (Settings$verbose >= 10) write(paste("Mutation rate", seq_len(Settings$num_of_sims), ":", u), stdout())

# census population size (empty data frame)
N = as.data.frame(matrix(nrow = Settings$num_of_sims,
                         ncol = total_periods))
names(N) = paste0("N_",seq_len(total_periods))

# ages of samples in generations (empty data frame)
ages = as.data.frame(matrix(nrow = Settings$num_of_sims,
                            ncol = Sample$size))
names(ages) = paste0(Sample$id,"_age")

for (sim in seq_len(Settings$num_of_sims)){
  if (Settings$verbose >= 10) write(paste("Simulation", sim), stdout())
  
  # sample demographic (N) trajectory
  N_trajectory = round(sample_param_trajectory(total_periods,
                                               Priors$pop_size_min,
                                               Priors$pop_size_max))
  N[sim,] = N_trajectory
  backward_demography = rev(N_trajectory[1:Model$periods_coalescence])
  forward_demography = N_trajectory[(Model$periods_coalescence+1):total_periods]
  if (Settings$verbose >= 10){
    write("Demography (N):", stdout())  
    write(paste(N_trajectory, collapse = " "), stdout())  
    write("Demography on coalscent period:", stdout())  
    write(paste(backward_demography, collapse = " "), stdout())  
  }
  
  # get ages of samples in number of generation in the past
  # ancient sample ages are taken from the calibrated age PDF
  sample_ages = get_sample_ages(Sample, gen_len[sim])
  ages[sim,] = sample_ages
  chrono_order = order(sample_ages)-1
  sampling_times_msprime = sort(unique(sample_ages))
  sampling_times_slim = rev(Model$generations_forward - sampling_times_msprime)
  sample_size_per_gen = rep(NA,length(sampling_times_msprime))
  for (i in seq_along(sampling_times_msprime)){
    sample_size_per_gen[i] = sum(sample_ages==sampling_times_msprime[i])
  }
  if (Settings$verbose >= 10){
    write("Ages of samples (in generations):", stdout())  
    write(paste(sample_ages, collapse = " "), stdout())  
    write("Chronological order of samples:", stdout())  
    write(paste(chrono_order, collapse = " "), stdout())  
    write("Sample size at each sampling time:", stdout())  
    write(paste(sample_size_per_gen, collapse = " "), stdout())  
  }

  # write sim *.ini file
  sim_ini <- list()
  sim_ini[["Simulation"]] = list(batch        = batch,
                                 sim          = sim)
  sim_ini[["Sample"]]     = list(ss           = paste(sample_size_per_gen, collapse=" "),
                                 msprime_ts   = paste(sampling_times_msprime, collapse=" "),
                                 chrono_order = paste(chrono_order, collapse=" ")) 
  sim_ini[["Demography"]] = list(N            = paste(backward_demography, collapse=" ")) 
  sim_ini[["Genome"]]     = list(mu           = u[sim],
                                 ttratio      = 2.0,
                                 seq_error    = 0.005)
  sim_ini[["Seeds"]]      = list(seed_coal    = seed_coal[sim],
                                 seed_mut     = seed_mut[sim])
  
  sim_ini_file <- paste0(batch_dir,"/sim_",sim,".ini")
  write.ini(sim_ini, sim_ini_file)
  
  # write source file for SLiM
  source4slim <- paste0("setSeed(", seed_forw[sim], ");\n",
                        "defineConstant(\"L\",", Genome$L, ");\n",
                        "defineConstant(\"N\",c(", paste(forward_demography, collapse = ","), "));\n",
                        "defineConstant(\"tc\",c(", paste(Model$times_of_change_forw, collapse = ","), "));\n",
                        "defineConstant(\"ts\",c(", paste(sampling_times_slim, collapse = ","), "));\n",
                        "defineConstant(\"ss\",c(", paste(rev(sample_size_per_gen), collapse = ","), "));\n",
                        "defineConstant(\"np\",", Model$periods_forward, ");\n",
                        "defineConstant(\"na\",", length(sampling_times_slim) - 1, ");\n",
                        "defineConstant(\"i\",", sim, ");\n",
                        "defineConstant(\"batch\",", batch, ");\n",
                        "defineConstant(\"project\",\"", Settings$project, "\");\n",
                        "defineConstant(\"ends\",c(", paste(Genome$slim_r_map$positions, collapse = ","), "));\n",
                        "defineConstant(\"rates\",c(", paste(Genome$slim_r_map$rates, collapse = ","), "));\n")
  write(source4slim, file = paste0(batch_dir, "/sim_", sim, ".eidos"))
}

params_RDS_file = paste0(batch_dir,"/ref_table_params.RDS")
params = cbind(N,u,gen_len,ages)
saveRDS(params,file=params_RDS_file)
