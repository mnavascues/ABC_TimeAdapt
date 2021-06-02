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
options <- set_options_type(read.ini(options_file))

# print script info to screen
if (!interactive() & !options$Settings$quiet){print_info_simulation()}

# set seed for random number generator
set.seed(options$Settings$seed)

# create results directory
dir.create("results", showWarnings = FALSE)
project_dir <- paste("results", options$Settings$project, sep="/")
dir.create(project_dir, showWarnings = FALSE)
batch_dir <- paste(project_dir, options$Settings$batch, sep="/")
dir.create(batch_dir, showWarnings = FALSE)

# read sample and genome information from tables in text files
Sample <- read_sample_info(options$Settings$sample_file)
if (!options$Settings$quiet) {cat("Sample:\n");(Sample)}
Genome <- read_genome_info(options$Settings$genome_file)
if (!options$Settings$quiet) {cat("Genome:\n");(Genome)}


# number of generations to simulate in forward (in SLiM)
times_of_change_forw     <- as.integer(seq(from = options$Model$generations_forward/options$Model$periods_forward,
                                           to   = options$Model$generations_forward-1,
                                           by   = options$Model$generations_forward/options$Model$periods_forward))
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
                                        options$Priors$gen_len_prior_min)
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
                           options$Priors$gen_len_prior_sh1,
                           options$Priors$gen_len_prior_sh2,
                           options$Priors$gen_len_prior_min,
                           options$Priors$gen_len_prior_max)
# census population size
sim_N <- sample_demography_from_prior(options$Settings$num_of_sims,
                                      options$Model$periods_forward+options$Model$periods_coalescence,
                                      options$Priors$pop_size_prior_min,
                                      options$Priors$pop_size_prior_max)
# sim_N_coal <- sim_N[,(1:options$Model$periods_coalescence)]
# sim_N_forw <- sim_N[,-(1:options$Model$periods_coalescence)]

# mutation rate
sim_u <- rlnorm(options$Settings$num_of_sims,
                log(options$Priors$mut_rate_prior_mean),
                options$Priors$mut_rate_prior_sd)

for (sim in seq_len(options$Settings$num_of_sims)){
  if (!options$Settings$quiet) cat(paste("\n\nSimulation",sim,"\n----------------------------------\n"))
  # simulate ages of aDNA from their calibrated age distribution
  sim_sample_time <- sample_ages_from_prior(Sample,
                                            options$Model$generations_forward,
                                            cal_age_PDF,
                                            gen_length=sim_gen_length[sim])
  # sample seeds for msprime and SLiM
  seed_coal <- round(runif(1,0,2^32-1))
  seed_forw <- round(runif(1,0,2^32-1))
  seed_mut  <- round(runif(1,0,2^32-1))

  # write sim *.ini file
  sim_ini <- list()
  sim_ini[["Simulation"]] <- list(#project  = options$Settings$project,
                                  #batch    = options$Settings$batch,
                                  sim      = sim)
  sim_ini[["Sample"]] <- list(ss           = paste(sim_sample_time$sample_sizes, collapse=" "),
                              # slim_ts      = paste(sim_sample_time$slim_ts, collapse=" "),
                              msprime_ts   = paste(sim_sample_time$msprime_ts, collapse=" "),
                              chrono_order = paste(sim_sample_time$chrono_order, collapse=" ")) 
  sim_ini[["Demography"]] <- list(N  = paste(sim_N[sim,], collapse=" "),
                                  tc = paste(times_of_change_forw, collapse=" ")) 
  sim_ini[["Genome"]] <- list(#L         = Genome$L,
                              #ends      = paste(Genome$rec_map_SLiM$ends, collapse=" "),
                              #rates     = paste(Genome$rec_map_SLiM$rates, collapse=" "),
                              mu        = sim_u[sim],
                              ttratio   = 2.0,
                              seq_error = 0.005)
  sim_ini[["Seeds"]] <- list(seed_coal = seed_coal,
                             #seed_forw = seed_forw,
                             seed_mut  = seed_mut)

  sim_ini_file <- paste0(batch_dir,"/sim_",sim,".ini")
  write.ini(sim_ini, sim_ini_file)
  
  # write source file for SLiM
  source4slim <- paste0("setSeed(",seed_forw,");\n",
                        "defineConstant(\"L\",",Genome$L,");\n",
                        "defineConstant(\"N\",c(",paste(sim_N[sim,-(1:options$Model$periods_coalescence)],collapse=","),"));\n",
                        "defineConstant(\"tc\",c(",paste(times_of_change_forw,collapse=","),"));\n",
                        "defineConstant(\"ts\",c(",paste(sim_sample_time$slim_ts,collapse=","),"));\n",
                        "defineConstant(\"ss\",c(",paste(rev(sim_sample_time$sample_sizes),collapse=","),"));\n",
                        "defineConstant(\"np\",",options$Model$periods_forward,");\n",
                        "defineConstant(\"na\",",sim_sample_time$na,");\n",
                        "defineConstant(\"i\",",sim,");\n",
                        "defineConstant(\"batch\",",options$Settings$batch,");\n",
                        "defineConstant(\"project\",\"",options$Settings$project,"\");\n",
                        "defineConstant(\"ends\",c(",paste(Genome$rec_map_SLiM$ends,collapse=","),"));\n",
                        "defineConstant(\"rates\",c(",paste(Genome$rec_map_SLiM$rates,collapse=","),"));\n")
  write(source4slim, file = paste0(batch_dir,"/sim_",sim,".eidos"))
  
}

