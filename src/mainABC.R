# Approximate Bayesian Computation inference 
# from ancient and modern DNA (NGS data)
# assuming a single continuous population
# using SLiM (and msprime via pyslim)

# Miguel Navascués
# Uppsala universitet & INRAE
# 2020

library(abcrf)
library(extraDistr)
library(rcarbon)
source("src/funABC.R")
set.seed(1234567890)

dir.create("results", showWarnings = FALSE)
project_name <- "test"
project_dir <- paste("results",project_name,sep="/")
dir.create(project_dir, showWarnings = FALSE)
batch_ID  <- 1 # an integer
batch_dir <- paste("results",project_name,batch_ID,sep="/")
dir.create(batch_dir, showWarnings = FALSE)

num_of_sims <- 1

# READ DATA INFO FROM FILE
sample_info_file <- "data/SampleInfoTest.txt"
Sample           <- read_sample_info(sample_info_file)

# GENOME INFO
# remove centromeres ?
# length of chromosomes ?
# transition/transversion ratio ?
chromosomes_limits <- read.table(file="data/chromosomes.txt")
L <- chromosomes_limits[22,2] # total genome length

# SET PRIORS
# rescaled beta distribution as prior for generation length:
prior_gen_length_shape1  <- 2
prior_gen_length_shape2  <- 1.465967
prior_gen_length_min     <- 26
prior_gen_length_max     <- 30
# log uniform prior for N (give values in natural scale)
prior_N_min             <- 10
prior_N_max             <- 200
# number of generations to simulate in forward (in SLiM)
num_of_gen_in_forw_sim   <- 400
num_of_periods_forw      <- 8
times_of_change_forw     <- seq(from = num_of_gen_in_forw_sim/num_of_periods_forw,
                                to   = num_of_gen_in_forw_sim-1,
                                by   = num_of_gen_in_forw_sim/num_of_periods_forw)

# write header of files
N_header <- c("sim",paste0("N",seq_len(num_of_periods_forw)))
Ne_header <- paste("sim",paste0("Ne",seq_len(num_of_periods_forw),collapse = " "))
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
check_ts_lower_gen_in_for_sim(num_of_gen_in_forw_sim,
                              Sample,
                              cal_age_PDF,
                              prior_gen_length_min)

# SAMPLE FROM PRIORS  
# Generation length = generation interval in years
sim_gen_length  <- rnsbeta(num_of_sims,
                           prior_gen_length_shape1,
                           prior_gen_length_shape2,
                           prior_gen_length_min,
                           prior_gen_length_max)
# census population size
sim_N <- sample_demography_from_prior(num_of_sims,
                                      num_of_periods_forw,
                                      prior_N_min,
                                      prior_N_max)

ref_table_N <- cbind(seq_len(num_of_sims),sim_N)
colnames(ref_table_N) <- N_header
saveRDS(ref_table_N,file = paste(batch_dir,"ref_table_N.RDS",sep="/"))
rm(ref_table_N);gc()

# mutation rate
sim_u <- rep(1.25e-08,num_of_sims)

sim<-1
#for (sim in seq_len(num_of_sims)){
  
  # simulate ages of aDNA from their calibrated age distribution
  sim_sample_time <- sample_ages_from_prior(Sample,
                                            num_of_gen_in_forw_sim,
                                            cal_age_PDF,
                                            gen_length=sim_gen_length[sim])
  #plot(times_of_change_forw,rep(1,length(times_of_change_forw)),xlim=c(0,num_of_gen_in_forw_sim))
  #points(sim_sample_time$slim_ts,rep(1,length(sim_sample_time$slim_ts)),col="red")
  
  # SLiM (forward simulation of last generations) 
  seed.slim <- round(runif(1,0,2^32-1))
  system2(command="slim",
          args=c("-d", paste0("N=\"c("), paste(sim_N[sim,],collapse=","), paste0(")\""),
                 "-d", paste0("tc=\"c("), paste(times_of_change_forw,collapse=","), paste0(")\""),
                 "-d", paste0("ts=\"c("), paste(sim_sample_time$slim_ts,collapse=","), paste0(")\""),
                 "-d", paste0("ss=\"c("), paste(rev(sim_sample_time$sample_sizes),collapse=","), paste0(")\""),
                 "-d", paste0("i=", sim),
                 "-d", paste0("batch_ID=",batch_ID),
                 "-d", paste0("project=\"'",project_name,"'\""),
                 "-d", paste0("np=", num_of_periods_forw),
                 "-d", paste0("na=", sim_sample_time$na),
                 "-d", paste0("L=", L),
                 "-s", seed.slim,
                 "src/model.demography.slim > /tmp/slimout.txt"))

    #}



# python (RECAPITATION + MUTATION + SEQUENCING + SUMSTATS)
seed.pyslim <- round(runif(1,0,2^32-1))
#system2(command="python3",
#        args=c("src/add.mutations.py",
#               i, N[i], mu[i], seed.pyslim, n0, na, ts,
#               " > /tmp/pyslimout.txt"))
  

system2(command="python3",
        args=c("src/add.mutations.py",
               "-i", sample_info_file,
               "-s", sim,
               "-b", batch_ID,
               "-p", project_name,
               "-t", sim_sample_time$msprime_ts,
               "-z", sim_sample_time$sample_sizes,
               "-d", seed.pyslim,
               "-n", sim_N[sim,1],
               "-u", sim_u[sim]))
               
