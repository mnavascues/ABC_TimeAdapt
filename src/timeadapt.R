#!/usr/bin/env Rscript

# Approximate Bayesian Computation inference 
# from ancient and modern DNA (NGS data)
# assuming a single continuous population
# using SLiM (and msprime via pyslim)

# Miguel Navascués
# Uppsala universitet & INRAE
# 2020

# Rscript timeadapt.R -h
# Rscript timeadapt.R -p test -q T -d 1234567890 -b 1 -s 1 -i data/SampleInfoTest.txt -g data/genome_test.txt -l 2 1.4 26 30 -n 10 200 -f 400 -w 8

library(argparser, quietly=TRUE)
library(abcrf, quietly=TRUE)
library(extraDistr, quietly=TRUE)
library(rcarbon, quietly=TRUE)
source("src/fun.R")

# read arguments from command line (gets default values in interactive)
argv <- get_arguments()

# write header
if (!interactive() & !argv$quiet){
  write("\n\n",stdout())
  write("###############################",stdout())
  write("# TimeAdapt                   #",stdout())
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

# read smaple and genome information from tables in text files
Sample <- read_sample_info(argv$sample_info_file)
Genome <- read_genome_info(argv$genome_info_file)

# number of generations to simulate in forward (in SLiM)
times_of_change_forw     <- seq(from = argv$num_of_gen_in_forw_sim/argv$num_of_periods_forw,
                                to   = argv$num_of_gen_in_forw_sim-1,
                                by   = argv$num_of_gen_in_forw_sim/argv$num_of_periods_forw)

# write header of files
N_header <- c("sim",paste0("N",seq_len(argv$num_of_periods_forw)))
Ne_header <- paste("sim",paste0("Ne",seq_len(argv$num_of_periods_forw),collapse = " "))
write(Ne_header,file=paste(batch_dir,"Ne.txt",sep="/"),append=F)

# calculate probability distribution curves for calibrated age of ancient samples
# (from 14C ages)
if (file.exists(paste(project_dir,"cal_age_PDF.RDS",sep="/"))){
  cal_age_PDF <- readRDS(file = paste(project_dir,"cal_age_PDF.RDS",sep="/"))
}else{
  cal_age_PDF <- get_sample_cal_age_PDF(Sample)
  saveRDS(cal_age_PDF,file = paste(project_dir,"cal_age_PDF.RDS",sep="/"))
}

# verify that samples cannot be older than the number of generations simulated in forward
check_ts_lower_gen_in_for_sim(argv$num_of_gen_in_forw_sim,
                              Sample,
                              cal_age_PDF,
                              argv$generation_length_prior_params[3])

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

ref_table_N <- cbind(seq_len(argv$num_of_sims),sim_N)
colnames(ref_table_N) <- N_header
saveRDS(ref_table_N,file = paste(batch_dir,"ref_table_N.RDS",sep="/"))
rm(ref_table_N); gc_out<-gc(verbose=FALSE)

# mutation rate
sim_u <- rep(5e-08,argv$num_of_sims) #rep(1.25e-08,argv$num_of_sims)

#sim<-1
for (sim in seq_len(argv$num_of_sims)){
  if (!argv$quiet) cat(paste("\n\nSimulation",sim,"\n----------------------------------\n"))
  # simulate ages of aDNA from their calibrated age distribution
  sim_sample_time <- sample_ages_from_prior(Sample,
                                            argv$num_of_gen_in_forw_sim,
                                            cal_age_PDF,
                                            gen_length=sim_gen_length[sim])
  #plot(times_of_change_forw,rep(1,length(times_of_change_forw)),xlim=c(0,argv$num_of_gen_in_forw_sim))
  #points(sim_sample_time$slim_ts,rep(1,length(sim_sample_time$slim_ts)),col="red")
  
  # SLiM (forward simulation of last generations) 
  seed_slim <- round(runif(1,0,2^32-1))
  args_slim <- c("-d", paste0("N=\"c("), paste(sim_N[sim,],collapse=","), paste0(")\""),
                 "-d", paste0("tc=\"c("), paste(times_of_change_forw,collapse=","), paste0(")\""),
                 "-d", paste0("ts=\"c("), paste(sim_sample_time$slim_ts,collapse=","), paste0(")\""),
                 "-d", paste0("ss=\"c("), paste(rev(sim_sample_time$sample_sizes),collapse=","), paste0(")\""),
                 "-d", paste0("i=", sim),
                 "-d", paste0("batch_ID=",argv$batch_ID),
                 "-d", paste0("project=\"'",argv$project_name,"'\""),
                 "-d", paste0("np=", argv$num_of_periods_forw),
                 "-d", paste0("na=", sim_sample_time$na),
                 "-d", paste0("L=", Genome$L),
                 "-d", paste0("ends=\"c("), paste(Genome$rec_map_SLiM[,1],collapse=","), paste0(")\""),
                 "-d", paste0("rates=\"c("), paste(Genome$rec_map_SLiM[,2],collapse=","), paste0(")\""),
                 "-s", seed_slim,
                 "src/model.demography.slim > /tmp/slimout.txt")
  system2(command="slim", args=args_slim)
  write(paste("\nslim",paste(args_slim,collapse=" "),"\n"), stdout())

  # python (RECAPITATION + MUTATION + SEQUENCING + SUMSTATS)
  seed_pyslim <- round(runif(1,0,2^32-1))
  args_pyslim <- c("timeadapt.py",
                   "-i", argv$sample_info_file,
                   "-g", argv$genome_info_file,
                   "-s", sim,
                   "-b", argv$batch_ID,
                   "-p", argv$project_name,
                   "-t", sim_sample_time$msprime_ts,
                   "-z", sim_sample_time$sample_sizes,
                   "-o", sim_sample_time$chrono_order,
                   "-d", seed_pyslim,
                   "-n", sim_N[sim,1],
                   "-u", sim_u[sim])
  system2(command = "python3", args = args_pyslim)
  write(paste("\npython3",paste(args_pyslim,collapse=" "),"\n"), stdout())
  
}# end loop simulations     
