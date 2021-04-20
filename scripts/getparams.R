print_info_simulation <- function(){
  write("#########################################",stdout())
  write("#                                       #",stdout())
  write("#      TimeAdapt - getparams            #",stdout())
  write("#      by Miguel de Navascués           #",stdout())
  write("#      Uppsala universitet & INRAE      #",stdout())
  write("#      miguel.navascues@inrae.fr        #",stdout())
  write("#                                       #",stdout())
  write("#########################################",stdout())
}

library(ini, quietly=TRUE)
library(extraDistr, quietly=TRUE)
source("scripts/timeadapt.R")

# read options from file
if (interactive()){ # if interactive uses test values
  options_file <- 'tests/input/config_project.ini'
}else{
  args = commandArgs(trailingOnly=TRUE)
  if (length(args)==0){
    stop("Options file (format *.ini) needs to be provided")
  }else{
    options_file <- args[1]
  }
}
options <- read.ini(options_file)
options$Settings$quiet            <- as.logical(options$Settings$quiet)
options$Model$generations_forward <- as.integer(options$Model$generations_forward) 
options$Model$periods_forward     <- as.integer(options$Model$periods_forward)
options$Model$periods_coalescence <- as.integer(options$Model$periods_coalescence)
options$Settings$num_of_sims      <- as.integer(options$Settings$num_of_sims)

# print script info to screen
if (!interactive() & !options$Settings$quiet){print_info_simulation()}

# set seed for random number generator
set.seed(options$Settings$seed)

# create results directory
dir.create("results", showWarnings = FALSE)
project_dir <- paste("results",options$Settings$project,sep="/")
dir.create(project_dir, showWarnings = FALSE)
batch_dir <- paste("results",options$Settings$project,options$Settings$batch,sep="/")
dir.create(batch_dir, showWarnings = FALSE)

# read sample and genome information from tables in text files
Sample <- read_sample_info(options$Settings$sample_file)
if (!options$Settings$quiet) {cat("Sample:\n");(Sample)}
Genome <- read_genome_info(options$Settings$genome_file)
if (!options$Settings$quiet) {cat("Genome:\n");(Genome)}


# number of generations to simulate in forward (in SLiM)
times_of_change_forw     <- seq(from = options$Model$generations_forward/options$Model$periods_forward,
                                to   = options$Model$generations_forward-1,
                                by   = options$Model$generations_forward/options$Model$periods_forward)
if (!options$Settings$quiet) {cat("Periods in forward:\n");(times_of_change_forw)}

# calculate probability distribution curves for calibrated age of ancient samples
# (from 14C ages)
if(any(Sample$is_ancient)){
  if (file.exists(paste(project_dir,"cal_age_PDF.RDS",sep="/"))){
    if (!options$Settings$quiet) cat("RDS file exists: importing calibrated age PDF\n")
    cal_age_PDF <- readRDS(file = paste(project_dir,"cal_age_PDF.RDS",sep="/"))
  }else{
    if (!options$Settings$quiet) cat("RDS file does not exists: calculating calibrated age PDF\n")
    cal_age_PDF <- get_sample_cal_age_PDF(Sample)
    saveRDS(cal_age_PDF,file = paste(project_dir,"cal_age_PDF.RDS",sep="/"))
  }  
}else{
  cal_age_PDF <- NULL
}


# verify that samples cannot be older than the number of generations simulated in forward
max_sample_age <- maximum_age_of_sample(Sample,
                                        cal_age_PDF,
                                        as.numeric(options$Priors$gen_len_prior_min))
if(max_sample_age+2 < options$Model$generations_forward){
  if (!options$Settings$quiet) cat("Length of simulation forward OK.\n")
}else{
  message("Insufficient length of forward time simulation. ",
          "It is necessary to simulate more than ", max_sample_age+2," generations ",
          "in forward to include all possible ages of samples.")
  quit(save="no",status=30)
}

# SAMPLE FROM PRIORS  
# Generation length = generation interval in years
sim_gen_length  <- rnsbeta(options$Settings$num_of_sims,
                           as.numeric(options$Priors$gen_len_prior_sh1),
                           as.numeric(options$Priors$gen_len_prior_sh2),
                           as.numeric(options$Priors$gen_len_prior_min),
                           as.numeric(options$Priors$gen_len_prior_max))
# census population size
sim_N <- sample_demography_from_prior(options$Settings$num_of_sims,
                                      options$Model$periods_forward+options$Model$periods_coalescence,
                                      as.numeric(options$Priors$pop_size_prior_min),
                                      as.numeric(options$Priors$pop_size_prior_max))
# sim_N_coal <- sim_N[,(1:options$Model$periods_coalescence)]
# sim_N_forw <- sim_N[,-(1:options$Model$periods_coalescence)]

# mutation rate
sim_u <- rlnorm(options$Settings$num_of_sims,
                log(as.numeric(options$Priors$mut_rate_prior_mean)),
                as.numeric(options$Priors$mut_rate_prior_sd))

for (sim in seq_len(options$Settings$num_of_sims)){
  if (!options$Settings$quiet) cat(paste("\n\nSimulation",sim,"\n----------------------------------\n"))
  # simulate ages of aDNA from their calibrated age distribution
  sim_sample_time <- sample_ages_from_prior(Sample,
                                            options$Model$generations_forward,
                                            cal_age_PDF,
                                            gen_length=sim_gen_length[sim])
  # sample seeds for msprime and SLiM
  seed_msprime <- round(runif(1,0,2^32-1))
  seed_slim    <- round(runif(1,0,2^32-1))

  # write sim *.ini file
  sim_ini <- list()
  sim_ini[["Simulation"]] <- list(project  = options$Settings$project,
                                  batch    = options$Settings$batch,
                                  i        = sim)
  sim_ini[["Sample"]] <- list(ss           = paste(sim_sample_time$sample_sizes, collapse=" "),
                              slim_ts      = paste(sim_sample_time$slim_ts, collapse=" "),
                              msprime_ts   = paste(sim_sample_time$msprime_ts, collapse=" "),
                              chrono_order = paste(sim_sample_time$chrono_order, collapse=" ")) 
  sim_ini[["Demography"]] <- list(N  = paste(sim_N[sim,], collapse=" "),
                                  tc = paste(times_of_change_forw, collapse=" ")) 
  sim_ini[["Genome"]] <- list(L         = Genome$L,
                              ends      = paste(Genome$rec_map_SLiM[,1], collapse=" "),
                              rates     = paste(Genome$rec_map_SLiM[,2], collapse=" "),
                              mu        = sim_u[sim],
                              ttratio   = 2.0,
                              seq_error = 0.005)
  sim_ini[["Seeds"]] <- list(seed_slim    = seed_slim,
                             seed_msprime = seed_msprime)
  sim_ini_file <- paste0(batch_dir,"/",options$Settings$project,"_sim_",sim,".ini")
  write.ini(sim_ini, sim_ini_file)
  
  # write command line for SLiM
  command_slim <- paste0("slim ",
                         " -d ", paste0("N=\"c("), paste(sim_N[sim,-(1:options$Model$periods_coalescence)],collapse=","), paste0(")\""),
                         " -d ", paste0("tc=\"c("), paste(times_of_change_forw,collapse=","), paste0(")\""),
                         " -d ", paste0("ts=\"c("), paste(sim_sample_time$slim_ts,collapse=","), paste0(")\""),
                         " -d ", paste0("ss=\"c("), paste(rev(sim_sample_time$sample_sizes),collapse=","), paste0(")\""),
                         " -d ", paste0("i=", sim),
                         " -d ", paste0("batch=",options$Settings$batch),
                         " -d ", paste0("project=\"'",options$Settings$project,"'\""),
                         " -d ", paste0("np=", options$Model$periods_forward),
                         " -d ", paste0("na=", sim_sample_time$na),
                         " -d ", paste0("L=", Genome$L),
                         " -d ", paste0("ends=\"c("), paste(Genome$rec_map_SLiM[,1],collapse=","), paste0(")\""),
                         " -d ", paste0("rates=\"c("), paste(Genome$rec_map_SLiM[,2],collapse=","), paste0(")\""),
                         " -s ", seed_slim,
                         " slim/simulation.slim > /tmp/slimout.txt")
  write(command_slim, file = paste0(batch_dir,"/slim_",sim,".sh"))

  # write command line for pyslim
  # command_pyslim <- paste("python", "python/msprimeNstats.py", options_file, sim_ini_file)
  # write(command_pyslim, file = paste0(batch_dir,"/pyslim_",sim,".sh"))
  
}

