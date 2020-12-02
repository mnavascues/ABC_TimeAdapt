#!/usr/bin/env Rscript

# TimeAdapt - simulation.R
# Miguel Navascués
# Uppsala universitet & INRAE
# 2020

library(argparser, quietly=TRUE)
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
if (!argv$quiet) cat("Sample:\n");(Sample)
Genome <- read_genome_info(argv$genome_info_file)
if (!argv$quiet) cat("Genome:\n");(Genome)

# number of generations to simulate in forward (in SLiM)
times_of_change_forw     <- seq(from = argv$num_of_gen_in_forw_sim/argv$num_of_periods_forw,
                                to   = argv$num_of_gen_in_forw_sim-1,
                                by   = argv$num_of_gen_in_forw_sim/argv$num_of_periods_forw)
if (!argv$quiet) cat("Periods in forward:\n");(times_of_change_forw)

# calculate probability distribution curves for calibrated age of ancient samples
# (from 14C ages)
if(any(Sample$is_ancient)){
  if (file.exists(paste(project_dir,"cal_age_PDF.RDS",sep="/"))){
    message("RDS file exists: importing calibrated age PDF")
    cal_age_PDF <- readRDS(file = paste(project_dir,"cal_age_PDF.RDS",sep="/"))
  }else{
    message("RDS file does not exists: calculating calibrated age PDF")
    cal_age_PDF <- get_sample_cal_age_PDF(Sample)
    saveRDS(cal_age_PDF,file = paste(project_dir,"cal_age_PDF.RDS",sep="/"))
  }  
}

