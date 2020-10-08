# functions for
# Approximate Bayesian Computation inference 
# from ancient and modern DNA (NGS data)
# assuming a single continuous population
# using SLiM (and msprime via pyslim)

# Miguel Navascués
# Uppsala universitet & INRAE
# 2020

get_arguments <- function(){
  ap <- arg_parser(description=paste0("Approximate Bayesian computation analysis ",
                                      "for joint inference of demography and selection ",
                                      "from temporal population genomic data."))
  ap <- add_argument(parser = ap,
                     arg = "--project_name",
                     help=paste0("Name of the project analysis. It is used as directory ",
                                 "in the path to the output files"),
                     default = "test",
                     type = "character",
                     short = "-p")
  ap <- add_argument(parser = ap,
                     arg = "--quiet",
                     help=paste0("Run on quiet mode."),
                     default = FALSE,
                     short = "-q")
  ap <- add_argument(parser = ap,
                     arg = "--seed",
                     help="Seed for random number generator",
                     type = "numeric",
                     short = "-d")
  ap <- add_argument(parser = ap,
                     arg = "--batch_ID",
                     help=paste0("Number used to identify the batch of simulations. ",
                                 "It is used as part of the output file names"),
                     type = "integer",
                     short = "-b")
  ap <- add_argument(parser = ap,
                     arg = "--num_of_sims",
                     help = paste0("Number of simulations to run."),
                     type = "integer",
                     short = "-s")
  ap <- add_argument(parser = ap,
                     arg = "--sample_info_file",
                     help = paste0("Text file with sample information organised as in ",
                                   "the example file in data folder (SampleInfoTest.txt)"),
                     type = "character",
                     short = "-i")
  ap <- add_argument(parser = ap,
                     arg = "--genome_info_file",
                     help = paste0("Text file with genome information organised as in ",
                                   "the example file in data folder (genome_test.txt)"),
                     type = "character",
                     short = "-g")
  ap <- add_argument(parser = ap,
                     arg = "--generation_length_prior_params",
                     help = paste0("Shape 1, shape 2, minimum and maximum, for a rescaled ",
                                   "beta distribution used as a prior for generation length"),
                     nargs = 4,
                     type = "numeric",
                     short = "-l")
  ap <- add_argument(parser = ap,
                     arg = "--population_size_prior_params",
                     help = paste0("Minimum and maximum, for a loguniform ",
                                   "distribution used as a prior for population size"),
                     nargs = 2,
                     type = "numeric",
                     short = "-n")
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
  
  
  
  if(! interactive()){
    argv <- parse_args(ap)
  }else{
    argv <- parse_args(ap, c("-d", "1234567890",              # seed
                             "-b", "1",                       # batch_ID
                             "-s", "3",                       # num_of_sims
                             "-i", "data/SampleInfoTest.txt", # sample_info_file
                             "-g", "data/genome_test.txt",    # genome_info_file
                             "-l", "2","1.465967","26","30",  # generation_length_prior_params
                             "-n", "10","200",                # population_size_prior_params 
                             "-f", "400",                     # num_of_gen_in_forw_sim
                             "-w", "8"                       # num_of_periods_forw
                             ))
  }
  if (!argv$quiet){
    print(ap)
    print_arguments(argv)
  }
  return(argv)
}


print_arguments <- function(argv){
  write(paste0("Seed: ",argv$seed), stdout())
  write(paste0("Project: ",argv$project_name), stdout())
  write(paste0("Batch: ",argv$batch_ID), stdout())
  write(paste0("Number of simulations: ",argv$num_of_sims), stdout())
  write(paste0("Sample input file: ",argv$sample_info_file), stdout())
  write(paste0("Genome input file: ",argv$genome_info_file), stdout())
  write(paste0("Generation lenght prior parameters (sh1 sh2 min max): ",
               paste0(argv$generation_length_prior_params, collapse=" ")), stdout())
  write(paste0("Population size prior parameters (min max): ",
               paste0(argv$population_size_prior_params, collapse=" ")), stdout())
  write(paste0("Generations in forward simulation: ",argv$num_of_gen_in_forw_sim), stdout())
  write(paste0("Number of periods in forward simulation: ",argv$num_of_periods_forw), stdout())
  #write(paste0("Parameter: ",argv$param), stdout())
  #write(paste0("Parameter: ",argv$param), stdout())
  #write(paste0("Parameter: ",argv$param), stdout())
  #write(paste0("Parameter: ",argv$param), stdout())
  #write(paste0("Parameter: ",argv$param), stdout())
}








read_sample_info <- function(file="data/SampleInfoTest.txt"){
  info <- read.table(file,header=T,stringsAsFactors=F,strip.white=T)
  # TODO add verification of header
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

read_genome_info <- function(file="data/genome_test.txt"){
  info <- read.table(file,header=T)
  # TODO add verification of header
  number_of_chromosomes <- nrow(info)
  rec_map_SLiM_rates <- numeric()
  rec_map_SLiM_ends <- numeric()
  #genome_interval_start <- numeric()
  #genome_interval_end <- numeric()
  # genome_interval_rate <- numeric()
  for (chr in seq_len(number_of_chromosomes)){
    rec_map_SLiM_rates <- c(rec_map_SLiM_rates, info$recombination_rate[chr], 0.5)
    rec_map_SLiM_ends  <- c(rec_map_SLiM_ends, 
                            info$chromosome_end[chr], 
                            info$chromosome_end[chr]+1)
    #genome_interval_start <- c(genome_interval_start,
    #                           info$chromosome_start[chr],
    #                           info$centromere_end[chr])
    #genome_interval_end <- c(genome_interval_end,
    #                         info$centromere_start[chr],
    #                         info$chromosome_end[chr])
    #genome_interval_rate <- c(genome_interval_rate,
    #                          info$recombination_rate[chr],
    #                          info$recombination_rate[chr])
  }
  return( list(number_of_chromosomes = number_of_chromosomes,
               L = info$chromosome_end[number_of_chromosomes], # total genome length
               rec_map_SLiM = cbind(ends=rec_map_SLiM_ends,
                                    rates=rec_map_SLiM_rates)
               #,genome_intervals = cbind(start=genome_interval_start,
               #                          end=genome_interval_end,
               #                         rate=genome_interval_rate)
               ) )
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
  ts <- 0
  for (k in which(Sample$is_ancient)){
    ts <- max(ts,cal_age_PDF[[k]]$calBP)
  }
  ts <- BPtoBCAD(ts)
  ts <- round(abs(ts-Sample$t0)/prior_gen_length_min)
  ts <- ts+2
  if (ts > num_of_gen_in_for_sim){
    write(paste0("ERROR: It is necessary to simulate more generations ",
                 "in forward to include all possible ages of samples"), stderr())
    quit("no", status=1)  }
}


sample_ages_from_prior <- function(Sample,
                                   num_of_gen_in_for_sim,
                                   cal_age_PDF,
                                   gen_length){

  ages_sim <- array(NA,Sample$size)
  t_max <- num_of_gen_in_for_sim
  for (k in which(Sample$is_ancient)){
    # sample from PDF of calibrated ages (BP) and transform to year (BC or AD)
    ages_sim[k] <- BPtoBCAD(sample(cal_age_PDF[[k]]$calBP,1,prob = cal_age_PDF[[k]]$PrDens))
  }
  ages_sim[Sample$is_modern] <- Sample$ageBCAD[Sample$is_modern]
  # transform to generations before "present" (present=year of most recent sample)
  ages_sim <- round(abs(ages_sim-Sample$t0)/gen_length)

  chrono_order <- order(ages_sim)-1 # Sample$id[order(ages_sim)]
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

sample_Ne_trajectory <- function(num_of_periods_forw,
                                 prior_Ne_min,
                                 prior_Ne_max){
  row_sim_Ne<-array(NA,num_of_periods_forw)
  row_sim_Ne[1] <-exp(runif(1,log(prior_Ne_min),log(prior_Ne_max)))
  for (i in 2:num_of_periods_forw){
    alpha <- runif(1,-1,1)
    row_sim_Ne[i] <- 10^max(min(log10(row_sim_Ne[i-1])+alpha,log10(prior_Ne_max)),log10(prior_Ne_min))
  }
  return(round(row_sim_Ne))
}

sample_demography_from_prior <- function(num_of_sims,
                                         num_of_periods_forw,
                                         prior_Ne_min,
                                         prior_Ne_max){
  
  sim_Ne <- t(replicate(num_of_sims,sample_Ne_trajectory(num_of_periods_forw,
                                                         prior_Ne_min,
                                                         prior_Ne_max)))
  return(sim_Ne)
}
