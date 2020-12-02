#!/usr/bin/env Rscript

# TimeAdapt - simulation.R
# Miguel Navascués
# Uppsala universitet & INRAE
# 2020

library(argparser, quietly=TRUE)
library(extraDistr, quietly=TRUE)
source("R/myfunctions.R")

# read arguments from command line (gets default values in interactive)
argv <- get_arguments()


# write header
if (!interactive() & !argv$quiet){
  write("\n\n",stdout())
  write("###############################",stdout())
  write("# TimeAdapt - simulation.R    #",stdout())
  write("# by Miguel Navascués         #",stdout())
  write("# Uppsala universitet & INRAE #",stdout())
  write("# miguel.navascues@inrae.fr   #",stdout())
  write("###############################",stdout())
  write("\n",stdout())
}

# set seed for random number generator
set.seed(argv$seed)

# create results directory
dir.create("results", showWarnings = FALSE)
project_dir <- paste("results",argv$project_name,sep="/")
dir.create(project_dir, showWarnings = FALSE)
batch_dir <- paste("results",argv$project_name,argv$batch_ID,sep="/")
dir.create(batch_dir, showWarnings = FALSE)

# read sample and genome information from tables in text files
Sample <- read_sample_info(argv$sample_info_file)
if (!argv$quiet) {cat("Sample:\n");(Sample)}
Genome <- read_genome_info(argv$genome_info_file)
if (!argv$quiet) {cat("Genome:\n");(Genome)}

# number of generations to simulate in forward (in SLiM)
times_of_change_forw     <- seq(from = argv$num_of_gen_in_forw_sim/argv$num_of_periods_forw,
                                to   = argv$num_of_gen_in_forw_sim-1,
                                by   = argv$num_of_gen_in_forw_sim/argv$num_of_periods_forw)
if (!argv$quiet) {cat("Periods in forward:\n");(times_of_change_forw)}

# calculate probability distribution curves for calibrated age of ancient samples
# (from 14C ages)
if(any(Sample$is_ancient)){
  if (file.exists(paste(project_dir,"cal_age_PDF.RDS",sep="/"))){
    if (!argv$quiet) cat("RDS file exists: importing calibrated age PDF\n")
    cal_age_PDF <- readRDS(file = paste(project_dir,"cal_age_PDF.RDS",sep="/"))
  }else{
    if (!argv$quiet) cat("RDS file does not exists: calculating calibrated age PDF\n")
    cal_age_PDF <- get_sample_cal_age_PDF(Sample)
    saveRDS(cal_age_PDF,file = paste(project_dir,"cal_age_PDF.RDS",sep="/"))
  }  
}else{
  cal_age_PDF <- NULL
}


# verify that samples cannot be older than the number of generations simulated in forward
if(check_ts_lower_gen_in_for_sim(argv$num_of_gen_in_forw_sim,
                                 Sample,
                                 cal_age_PDF,
                                 argv$generation_length_prior_params[3])){
  if (!argv$quiet) cat("Length of simulation forward OK\n")
}else{
  quit("no",status=30)
} 

# SAMPLE FROM PRIORS  
# Generation length = generation interval in years
sim_gen_length  <- rnsbeta(argv$num_of_sims,
                           argv$generation_length_prior_params[1],
                           argv$generation_length_prior_params[2],
                           argv$generation_length_prior_params[3],
                           argv$generation_length_prior_params[4])
# census population size
sim_N <- sample_demography_from_prior(argv$num_of_sims,
                                      argv$num_of_periods_forw,
                                      argv$population_size_prior_params[1],
                                      argv$population_size_prior_params[2])
# mutation rate
sim_u <- rlnorm(argv$num_of_sims,
                log(argv$mutation_rate_prior_params[1]),
                argv$mutation_rate_prior_params[2])

for (sim in seq_len(argv$num_of_sims)){
  if (!argv$quiet) cat(paste("\n\nSimulation",sim,"\n----------------------------------\n"))
  # simulate ages of aDNA from their calibrated age distribution
  sim_sample_time <- sample_ages_from_prior(Sample,
                                            argv$num_of_gen_in_forw_sim,
                                            cal_age_PDF,
                                            gen_length=sim_gen_length[sim])
  # sample seeds for SLiM and pyslim
  seed_slim <- round(runif(1,0,2^32-1))
  seed_pyslim <- round(runif(1,0,2^32-1))
  
  # write command line for SLiM
  command_slim <- paste0("slim ",
                         " -d ", paste0("N=\"c("), paste(sim_N[sim,],collapse=","), paste0(")\""),
                         " -d ", paste0("tc=\"c("), paste(times_of_change_forw,collapse=","), paste0(")\""),
                         " -d ", paste0("ts=\"c("), paste(sim_sample_time$slim_ts,collapse=","), paste0(")\""),
                         " -d ", paste0("ss=\"c("), paste(rev(sim_sample_time$sample_sizes),collapse=","), paste0(")\""),
                         " -d ", paste0("i=", sim),
                         " -d ", paste0("batch_ID=",argv$batch_ID),
                         " -d ", paste0("project=\"'",argv$project_name,"'\""),
                         " -d ", paste0("np=", argv$num_of_periods_forw),
                         " -d ", paste0("na=", sim_sample_time$na),
                         " -d ", paste0("L=", Genome$L),
                         " -d ", paste0("ends=\"c("), paste(Genome$rec_map_SLiM[,1],collapse=","), paste0(")\""),
                         " -d ", paste0("rates=\"c("), paste(Genome$rec_map_SLiM[,2],collapse=","), paste0(")\""),
                         " -s ", seed_slim,
                         " slim/simulation.slim > /tmp/slimout.txt")
  write(command_slim, file = paste0(batch_dir,"/slim_",sim,".sh"))

  # write command line for pyslim
  command_pyslim <- paste("python3", "python/msprimeNstats.py",
                           "-i", argv$sample_info_file,
                           "-g", argv$genome_info_file,
                           "-s", sim,
                           "-b", argv$batch_ID,
                           "-p", argv$project_name,
                           "-t", paste(sim_sample_time$msprime_ts, collapse=" "),
                           "-z", paste(sim_sample_time$sample_sizes, collapse=" "),
                           "-o", paste(sim_sample_time$chrono_order, collapse=" "),
                           "-d", seed_pyslim,
                           "-n", sim_N[sim,1],
                           "-u", sim_u[sim])
  write(command_pyslim, file = paste0(batch_dir,"/pyslim_",sim,".sh"))

}

