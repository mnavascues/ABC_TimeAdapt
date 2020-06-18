# functions for
# Approximate Bayesian Computation inference 
# from ancient and modern DNA (NGS data)
# assuming a single continuous population
# using SLiM (and msprime via pyslim)

# Miguel Navascués
# Uppsala universitet & INRAE
# 2020

read_sample_info <- function(file="data/SampleInfoTest.txt"){
  info <- read.table(file,header=T,stringsAsFactors=F)
  return( list(id            = info$sample,
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
  if(ts > num_of_gen_in_for_sim) stop("It is necessary to simulate more generations in
forward to include all possible ages of samples")
  return(ts <= num_of_gen_in_for_sim)
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

  chrono_order <- Sample$id[order(ages_sim)]
  a <- sort(as.numeric(levels(as.factor(ages_sim))))
  s <- numeric()
  for (age in a){
    s<-c(s,sum(ages_sim==age))
  }
  
  return(list(na           = sum(a!=0),       # number of ancient samples
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
