# TimeAdapt - myfunctions.R
# Miguel Navascués
# Uppsala universitet & INRAE
# 2020

library(rcarbon, quietly=TRUE)

get_arguments <- function(){
  ap <- arg_parser(description = paste("Approximate Bayesian computation analysis",
                                       "for joint inference of demography and selection",
                                       "from temporal population genomic data.",
                                       "Reading sample info and sampling parameters from priors"))
  ap <- add_argument(parser = ap,
                     arg = "--quiet",
                     help = paste0("Run on quiet mode."),
                     default = FALSE,
                     short = "-q")
  ap <- add_argument(parser = ap,
                     arg = "--seed",
                     help = "Seed for random number generator ",
                     type = "numeric",
                     short = "-d")
  ap <- add_argument(parser = ap,
                     arg = "--project_name",
                     help=paste0("Name of the project analysis. It is used as directory ",
                                 "in the path to the output files"),
                     default = "test",
                     type = "character",
                     short = "-p")
  ap <- add_argument(parser = ap,
                     arg = "--batch_ID",
                     help=paste0("Number used to identify the batch of simulations. ",
                                 "It is used as part of the output file names"),
                     type = "integer",
                     short = "-b")
  ap <- add_argument(parser = ap,
                     arg = "--sample_info_file",
                     help = paste0("Text file with sample information organised as in ",
                                   "the example file in data folder (sample_info_test.txt)"),
                     type = "character",
                     short = "-i")
  ap <- add_argument(parser = ap,
                     arg = "--genome_info_file",
                     help = paste0("Text file with genome information organised as in ",
                                   "the example file in data folder (genome_info_test.txt)"),
                     type = "character",
                     short = "-g")
  ap <- add_argument(parser = ap,
                     arg = "--num_of_gen_in_forw_sim",
                     help = paste0("Number of generations to run in forward (SLiM). ",
                                   "It must be big enough as to cover the whole range of ",
                                   "possible ages for samples."),
                     type = "integer",
                     short = "-f")
  ap <- add_argument(parser = ap,
                     arg = "--num_of_periods_forw",
                     help = paste0("Number of periods with different parameter values ",
                                   "during the simulation forward (SLiM)."),
                     type = "integer",
                     short = "-w")
  ap <- add_argument(parser = ap,
                     arg = "--generation_length_prior_params",
                     help = paste0("Shape 1, shape 2, minimum and maximum, for a rescaled ",
                                   "beta distribution used as a prior for generation length"),
                     nargs = 4,
                     type = "numeric",
                     short = "-l")
  ap <- add_argument(parser = ap,
                     arg = "--num_of_sims",
                     help = paste0("Number of simulations to run."),
                     type = "integer",
                     short = "-s")
  ap <- add_argument(parser = ap,
                     arg = "--population_size_prior_params",
                     help = paste0("Minimum and maximum, for a loguniform ",
                                   "distribution used as a prior for population size"),
                     nargs = 2,
                     type = "numeric",
                     short = "-n")
  if(! interactive()){
    f_argv <- parse_args(ap)
  }else{
    f_argv <- parse_args(ap, c("-q", "FALSE",                     # quiet
                               "-d", "1234567890",                # seed
                               "-p", "test",                      # project_name
                               "-b", "1",                         # batch_ID
                               "-i", "data/sample_info_test.txt", # sample_info_file
                               "-g", "data/genome_info_test.txt", # genome_info_file
                               "-f", "400",                       # num_of_gen_in_forw_sim
                               "-w", "8",                         # num_of_periods_forw
                               "-l", "2","1.465967","26","30",    # generation_length_prior_params
                               "-s", "3",                         # num_of_sims
                               "-n", "10","200"                   # population_size_prior_params
                               ))
  }
  if (!f_argv$quiet){
    print(ap)
    print_arguments(f_argv)
  }
  return(f_argv)
}

print_arguments <- function(f_argv){
  write(paste0("Quiet: ",f_argv$quiet), stdout())
  write(paste0("Seed: ",f_argv$seed), stdout())
  write(paste0("Project: ",f_argv$project_name), stdout())
  write(paste0("Batch: ",f_argv$batch_ID), stdout())
  write(paste0("Sample file: ",f_argv$sample_info_file), stdout())
  write(paste0("Genome file: ",f_argv$genome_info_file), stdout())
  write(paste0("Generations in forward: ",f_argv$num_of_gen_in_forw_sim), stdout())
  write(paste0("Number of periods in forward: ",f_argv$num_of_periods_forw), stdout())
  write(paste0("Parameters for generation length prior: ",f_argv$generation_length_prior_params), stdout())
  write(paste0("Number of simulations: ",f_argv$num_of_sims), stdout())
  write(paste0("Parameters for poplation size prior: ",f_argv$population_size_prior_params), stdout())
}


check_file_header <- function(expected_header,file_header){
  missing=!is.element(expected_header,file_header)
  if (any(missing)){
    stop(paste("Missing columns in input file:",expected_header[missing]))
  }
  return(TRUE)
}

read_sample_info <- function(file="data/sample_info_test.txt"){
  info <- read.table(file,header=T,stringsAsFactors=F,strip.white=T)
  expected_header <- c("sampleID","age14C","age14Cerror","year","coverage","damageRepair","groups")
  if(check_file_header(expected_header, file_header = colnames(info))){
    if(all(is.character(info$sampleID),
           is.numeric(info$age14C),
           is.numeric(info$age14Cerror),
           is.numeric(info$year),
           is.numeric(info$coverage),
           is.logical(info$damageRepair))){
      return( list(id            = info$sampleID,
                   age14C        = info$age14C,
                   age14Cerror   = info$age14Cerror,
                   ageBCAD       = info$year, 
                   coverage      = info$coverage,
                   is_modern     = !is.na(info$year), 
                   is_ancient    = !is.na(info$age14C),
                   is_dr         = info$damageRepair,
                   total_ancient = sum(!is.na(info$age14C)),
                   size          = nrow(info),
                   t0            = max(info$year,na.rm=T) ) )
    }
    stop(paste("Wrong data type in file:",file))
  }else{
    quit("no",status=10)
  }
}


read_genome_info <- function(file="data/genome_info_test.txt"){
  info <- read.table(file,header=T)
  #expected_header <- c("ID","chromosome_start","chromosome_end","centromere_start","centromere_end","recombination_rate")
  expected_header <- c("chromosome_end","recombination_rate") # for the moment these are the only columns used
  if(check_file_header(expected_header, file_header = colnames(info))){
    number_of_chromosomes <- nrow(info)
    rec_map_SLiM_rates <- numeric()
    rec_map_SLiM_ends <- numeric()
    for (chr in seq_len(number_of_chromosomes)){
      rec_map_SLiM_rates <- c(rec_map_SLiM_rates, info$recombination_rate[chr], 0.5)
      rec_map_SLiM_ends  <- c(rec_map_SLiM_ends, 
                              info$chromosome_end[chr], 
                              info$chromosome_end[chr]+1)
    }
    return( list(number_of_chromosomes = number_of_chromosomes,
                 L = info$chromosome_end[number_of_chromosomes], # total genome length
                 rec_map_SLiM = cbind(ends=rec_map_SLiM_ends,
                                      rates=rec_map_SLiM_rates) ))
  }else{
    quit("no",status=20)
  }
}


get_sample_cal_age_PDF <- function(Sample,calibration_curve='shcal13'){
  cal_age_PDF <- vector("list",Sample$size)
  for (k in which(Sample$is_ancient)){
    cal_age_dist <- calibrate(x         = Sample$age14C[k],
                              error     = Sample$age14Cerror[k],
                              calCurves = calibration_curve,
                              verbose   = F)
    cal_age_PDF[[k]] <-  cal_age_dist$grids$`1`
  }
  return(cal_age_PDF)
}


check_ts_lower_gen_in_for_sim <- function(num_of_gen_in_for_sim,
                                          Sample,
                                          cal_age_PDF,
                                          prior_gen_length_min){
  # TODO: modify for the case of only modern samples (from different years)
  ts <- 0
  for (k in which(Sample$is_ancient)){
    ts <- max(ts,cal_age_PDF[[k]]$calBP)
  }
  ts <- BPtoBCAD(ts)
  ts <- min(c(ts,Sample$ageBCAD),na.rm=TRUE) # check if "modern" samples are not older, just in case
  ts <- round(abs(ts-Sample$t0)/prior_gen_length_min)
  ts <- ts+2
  if (ts > num_of_gen_in_for_sim){
    stop(paste("Insufficient length of forward time simulation.",
               "It is necessary to simulate more generations",
               "in forward to include all possible ages of samples"))
  }
  return(TRUE)
}

sample_N_trajectory <- function(num_of_periods_forw,
                                prior_N_min,
                                prior_N_max){
  row_sim_N<-array(NA,num_of_periods_forw)
  row_sim_N[1] <-exp(runif(1,log(prior_N_min),log(prior_N_max)))
  for (i in 2:num_of_periods_forw){
    alpha <- runif(1,-1,1)
    row_sim_N[i] <- 10^max(min(log10(row_sim_N[i-1])+alpha,log10(prior_N_max)),log10(prior_N_min))
  }
  return(round(row_sim_N))
}

sample_demography_from_prior <- function(num_of_sims,
                                         num_of_periods_forw,
                                         prior_N_min,
                                         prior_N_max){
  
  sim_N <- t(replicate(num_of_sims,sample_N_trajectory(num_of_periods_forw,
                                                       prior_N_min,
                                                       prior_N_max)))
  return(sim_N)
}

