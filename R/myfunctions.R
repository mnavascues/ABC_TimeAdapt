# TimeAdapt - myfunctions.R
# Miguel de Navascués
# Uppsala universitet & INRAE
# 2021

library(rcarbon, quietly=TRUE)

get_arguments <- function(test=FALSE, test_arg=NULL){
  ap <- arg_parser(description = paste("TimeAdapt. Approximate Bayesian computation analysis",
                                       "for joint inference of demography and selection",
                                       "from temporal population genomic data."))
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
                                   "the example file in tests folder (sample_info_test.txt)"),
                     type = "character",
                     short = "-i")
  ap <- add_argument(parser = ap,
                     arg = "--genome_info_file",
                     help = paste0("Text file with genome information organised as in ",
                                   "the example file in tests folder (genome_info_test.txt)"),
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
  ap <- add_argument(parser = ap,
                     arg = "--mutation_rate_prior_params",
                     help = paste0("Mean and SD, for a lognormal ",
                                   "distribution used as a prior for mutation rate"),
                     nargs = 2,
                     type = "numeric",
                     short = "-u")
  if(!test){
    f_argv <- parse_args(ap)
  }else{
    f_argv <- parse_args(ap, test_arg)
  }
  if (!f_argv$quiet){
    print(ap)
    #print_arguments(f_argv)
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
  write(paste0("Parameters for mutation rate prior: ",f_argv$mutation_rate_prior_params), stdout())
}





check_file_header <- function(expected_header,file_header){
  missing=!is.element(expected_header,file_header)
  if (any(missing)){
    stop(paste("Missing columns in input file:",expected_header[missing]))
  }
  return(TRUE)
}

read_sample_info <- function(file){
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
    } else {
      stop(paste("Wrong data type in file:",file))
    }
  }else{
    quit("check_file_header error on read_sample_info",status=10)
  }
}


read_genome_info <- function(file){
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
    quit("check_file_header error on read_genome_info",status=20)
  }
}


get_sample_cal_age_PDF <- function(f_Sample,calibration_curve='shcal13'){
  cal_age_PDF <- vector("list",f_Sample$size)
  for (k in which(f_Sample$is_ancient)){
    cal_age_dist <- rcarbon::calibrate(x         = f_Sample$age14C[k],
                                       error     = f_Sample$age14Cerror[k],
                                       calCurves = calibration_curve,
                                       verbose   = F)
    cal_age_PDF[[k]] <-  cal_age_dist$grids$`1`
  }
  return(cal_age_PDF)
}

maximum_age_of_sample <- function(f_Sample,
                                  f_cal_age_PDF,
                                  prior_gen_length_min){
  ts <- 0
  for (k in which(f_Sample$is_ancient)){
    ts <- max(ts,f_cal_age_PDF[[k]]$calBP)
  }
  ts <- rcarbon::BPtoBCAD(ts)
  ts <- min(ts,f_Sample$ageBCAD[which(f_Sample$is_modern)])
  ts <- round(abs(ts-f_Sample$t0)/prior_gen_length_min)
  return(ts)
}


# function no longer in use substituted by maximum_age_of_sample()
check_ts_lower_gen_in_for_sim <- function(num_of_gen_in_for_sim,
                                          f_Sample,
                                          f_cal_age_PDF,
                                          prior_gen_length_min){
  ts <- 0
  for (k in which(f_Sample$is_ancient)){
    ts <- max(ts,f_cal_age_PDF[[k]]$calBP)
  }
  ts <- BPtoBCAD(ts)
  ts <- min(c(ts,f_Sample$ageBCAD),na.rm=TRUE) # check if "modern" samples are not older, just in case
  ts <- round(abs(ts-f_Sample$t0)/prior_gen_length_min)
  ts <- ts+2
  if (ts > num_of_gen_in_for_sim){
    message("Insufficient length of forward time simulation. ",
            "It is necessary to simulate more than ", ts," generations ",
            "in forward to include all possible ages of samples.")
    return(FALSE)
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




sample_ages_from_prior <- function(f_Sample,
                                   num_of_gen_in_for_sim,
                                   f_cal_age_PDF,
                                   gen_length){
  
  # TODO: modify so it can get samples with only modern DNA
  
  ages_sim <- array(NA,f_Sample$size)
  t_max <- num_of_gen_in_for_sim
  for (k in which(f_Sample$is_ancient)){
    # sample from PDF of calibrated ages (BP) and transform to year (BC or AD)
    ages_sim[k] <- BPtoBCAD(sample(f_cal_age_PDF[[k]]$calBP,1,prob = f_cal_age_PDF[[k]]$PrDens))
  }
  ages_sim[f_Sample$is_modern] <- f_Sample$ageBCAD[f_Sample$is_modern]
  # transform to generations before "present" (present=year of most recent sample)
  ages_sim <- round(abs(ages_sim-f_Sample$t0)/gen_length)
  
  chrono_order <- order(ages_sim)-1 # f_Sample$id[order(ages_sim)]
  a <- sort(as.numeric(levels(as.factor(ages_sim))))
  s <- numeric()
  for (age in a){
    s<-c(s,sum(ages_sim==age))
  }
  
  return(list(na           = sum(a!=0),       # number of non-present sample times
              t_max        = t_max,
              slim_ts      = sort(t_max-a),   # time of samples in SLiM
              msprime_ts   = a,               # time of samples in tree sequences  
              sample_sizes = s,               # sample size per generation
              chrono_order = chrono_order ) ) # chronological order of samples
  
  
}
