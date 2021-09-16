# TimeAdapt - timeadapt.R
# Miguel de Navascués
# Uppsala universitet & INRAE
# 2021

library(rcarbon, quietly=TRUE)

set_options_type <- function(options){
  
  options$Settings$seed              <- as.integer(options$Settings$seed)
  options$Settings$quiet             <- as.logical(options$Settings$quiet)
  options$Settings$num_of_sims       <- as.integer(options$Settings$num_of_sims)

  options$Model$generations_forward  <- as.integer(options$Model$generations_forward) 
  options$Model$periods_forward      <- as.integer(options$Model$periods_forward)
  options$Model$periods_coalescence  <- as.integer(options$Model$periods_coalescence)

  options$Priors$gen_len_prior_min   <- as.numeric(options$Priors$gen_len_prior_min)
  options$Priors$gen_len_prior_sh1   <- as.numeric(options$Priors$gen_len_prior_sh1)
  options$Priors$gen_len_prior_sh2   <- as.numeric(options$Priors$gen_len_prior_sh2)
  options$Priors$gen_len_prior_min   <- as.numeric(options$Priors$gen_len_prior_min)
  options$Priors$gen_len_prior_max   <- as.numeric(options$Priors$gen_len_prior_max)
  options$Priors$pop_size_prior_min  <- as.numeric(options$Priors$pop_size_prior_min)
  options$Priors$pop_size_prior_max  <- as.numeric(options$Priors$pop_size_prior_max)  
  options$Priors$mut_rate_prior_mean <- as.numeric(options$Priors$mut_rate_prior_mean)  
  options$Priors$mut_rate_prior_sd   <- as.numeric(options$Priors$mut_rate_prior_sd)

  return(options)
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
  header_ok <- FALSE
  header_ok <- check_file_header(expected_header, file_header = colnames(info))
  if(header_ok){
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
    quit(save="no",status=10)
  }
}


read_genome_info <- function(file){
  info <- read.table(file,header=T)
  #expected_header <- c("ID","chromosome_start","chromosome_end","centromere_start","centromere_end","recombination_rate")
  expected_header <- c("Chromosome","Position","Recombination_rate") # for the moment these are the only columns used
  header_ok <- FALSE
  header_ok <- check_file_header(expected_header, file_header = colnames(info))
  if(header_ok){
    nchr <- nlevels(as.factor(info$Chromosome))
    chr_ends_index <- c(which(diff(info$Chromosome)!=0), length(info$Chromosome))
    rescaling_values <- c(0,cumsum(info$Position[chr_ends_index]))  
    chromo = 1
    for (i in seq_along(info$Position)){
      info$Position[i] = info$Position[i] + rescaling_values[chromo]
      if (any(i == chr_ends_index)) chromo = chromo + 1
    }
    rec_map_SLiM_ends <- numeric()
    rec_map_SLiM_rates <- numeric()
    for (chromo in seq_len(nchr)){
      rec_map_SLiM_ends <- c(rec_map_SLiM_ends,info$Position[which(info$Chromosome==chromo)]-1)
      rec_map_SLiM_rates <- c(rec_map_SLiM_rates,info$Recombination_rate[which(info$Chromosome==chromo)])
      if (chromo!=max(info$Chromosome)){
        rec_map_SLiM_ends <- c(rec_map_SLiM_ends,info$Position[chr_ends_index[chromo]])
        rec_map_SLiM_rates <- c(rec_map_SLiM_rates,0.5)
      }
    }
    return( list(number_of_chromosomes = nchr,
                 L = as.integer(rec_map_SLiM_ends[length(rec_map_SLiM_ends)]), # total genome length
                 rec_map_SLiM = data.frame(ends=as.integer(rec_map_SLiM_ends),
                                           rates=rec_map_SLiM_rates) ))
  }else{
    quit(save="no",status=20)
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
