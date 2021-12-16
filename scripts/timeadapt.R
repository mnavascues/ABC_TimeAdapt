#    TimeAdapt: joint inference of demography and selection
#    Copyright (C) 2021  Miguel de Navascués, Uppsala universitet, INRAE
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <https://www.gnu.org/licenses/>.
#

library(rcarbon, quietly = TRUE)

### PRINT INFO ######################################################################################
print_info = function(script_name, verbose, project = NA, batch = NA, sim = NA)
{
  if (verbose >= 1) write("#########################################",stdout())
  if (verbose >= 0)
  {
    info_message = paste0("TimeAdapt - ", script_name)
    if (!is.na(project)) info_message = paste0(info_message, " - project ", project)
    if (!is.na(batch)) info_message = paste0(info_message, " - batch ", batch)
    if (!is.na(sim)) info_message = paste0(info_message, " - simulation ", sim)
    write(info_message, stdout())
  }
  if (verbose >= 1)
  {
    write("by Miguel de Navascués", stdout())
    write("INRAE & Uppsala universitet", stdout())
    write("miguel.navascues@inrae.fr", stdout())
    write("#########################################", stdout())
  } 
}

### GET PROJECT OPTIONS ###############################################################################
get_project_options = function(options_file){
  options = read.ini(options_file)
  # SETTINGS
  #options$Settings$project        = options$Settings$project
  options$Settings$num_of_batches = as.integer(options$Settings$num_of_batches)
  #options$Settings$genome_file    = options$Settings$genome_file
  #options$Settings$sample_file    = options$Settings$sample_file
  options$Settings$verbose        = as.integer(options$Settings$verbose)
  options$Settings$num_of_sims    = as.integer(options$Settings$num_of_sims)
  options$Settings$seed           = as.integer(options$Settings$seed)
  # MODEL
  options$Model$generations_forward     = as.integer(options$Model$generations_forward)
  options$Model$periods_forward         = as.integer(options$Model$periods_forward)
  options$Model$generations_coalescence = as.integer(options$Model$generations_coalescence)
  options$Model$periods_coalescence     = as.integer(options$Model$periods_coalescence)
  #options$Model$calibration_curve       = options$Model$calibration_curve
  # PRIORS
  options$Priors$gen_len_sh1   = as.numeric(options$Priors$gen_len_sh1)
  options$Priors$gen_len_sh2   = as.numeric(options$Priors$gen_len_sh2)
  options$Priors$gen_len_min   = as.numeric(options$Priors$gen_len_min)
  options$Priors$gen_len_max   = as.numeric(options$Priors$gen_len_max)
  options$Priors$pop_size_min  = as.numeric(options$Priors$pop_size_min)
  options$Priors$pop_size_max  = as.numeric(options$Priors$pop_size_max)  
  options$Priors$mut_rate_mean = as.numeric(options$Priors$mut_rate_mean)  
  options$Priors$mut_rate_sd   = as.numeric(options$Priors$mut_rate_sd)
  return(options)
}

### GET TIMES OF CHANGE ######################################################################################
get_times_of_change = function(total_length, number_of_periods, mode = "regular", base = 2){
  if (number_of_periods == 1) return(NA)
  if (total_length <= number_of_periods) stop("number of periods must be lower than length of simulation")
  if (mode == "exponential"){
    t0 = total_length / ((base^(number_of_periods - 1) - 1) / (base - 1))
    if (t0 < 1) stop("incompatible value, try less periods or longer simulation")
    times_of_change = round(cumsum(t0 * base^(seq_len(number_of_periods - 1) - 1)))
  }else{
    if (mode != "regular") warning("Wrong value for mode in get_times_of_change(); using mode='regular' instead")
    times_of_change = round(total_length / number_of_periods * seq_len(number_of_periods - 1))
  }
  if (length(times_of_change) != number_of_periods - 1) stop("wrong number of times")
  if (length(unique(times_of_change)) < length(times_of_change)) stop('repeated values in times of change')
  return(times_of_change)
}

### CHECK FILE HEADER ######################################################################################
check_file_header = function(expected_header, file_header){
  missing =! is.element(expected_header, file_header)
  if (any(missing)) return(list(ok = FALSE, missing = expected_header[missing]))
  else              return(list(ok = TRUE,  missing = NA)) 
}

### READ SAMPLE INFO ######################################################################################
read_sample_info = function(file){
  info = read.table(file,header=T,stringsAsFactors=F,strip.white=T)
  expected_header = c("sampleID", "age14C", "age14Cerror", "calCurves", "year", "coverage", "damaged", "groups")
  header_check = check_file_header(expected_header, file_header = colnames(info))
  if(header_check$ok){
    if(all(all(is.na(info$age14C)) | is.numeric(info$age14C),
           all(is.na(info$age14Cerror)) | is.numeric(info$age14Cerror),
           all(is.na(info$calCurves)) | is.character(info$calCurves),
           all(is.na(info$year)) | is.numeric(info$year),
           is.numeric(info$coverage),
           !any(is.na(info$damaged)) & is.logical(info$damaged),
           is.character(info$groups))){
      if (any(is.na(info$year) == is.na(info$age14C))) stop("Samples with undefined ages")
      return( list(id            = info$sampleID,
                   age14C        = info$age14C,
                   age14Cerror   = info$age14Cerror,
                   calCurves     = info$calCurves,
                   ageBCAD       = info$year, 
                   coverage      = info$coverage,
                   # is_modern     = !is.na(info$year), 
                   is_ancient    = !is.na(info$age14C),
                   is_damaged    = info$damaged,
                   total_ancient = sum(!is.na(info$age14C)),
                   size          = nrow(info),
                   groups        = info$groups) )
    } else {
      stop(paste("Wrong data type in file:",file))
    }
  }else{
    stop(paste("Missing columns in input file:",header_check$missing))
  }
}

### READ GENOME INFO ######################################################################################
read_genome_info = function(file){
  info = read.table(file,header=T)
  expected_header = c("Chromosome","Position","Recombination_rate")
  header_check = check_file_header(expected_header, file_header = colnames(info))
  if(header_check$ok){
    rates            = info$Recombination_rate
    positions        = info$Position
    nchr             = nlevels(as.factor(info$Chromosome))
    chr_ends_index   = c(which(diff(info$Chromosome)!=0), length(info$Chromosome))
    rescaling_values = c(0,cumsum(info$Position[chr_ends_index]))  
    chromo = 1
    for (i in seq_along(info$Position)){
      positions[i] = info$Position[i] + rescaling_values[chromo]
      if (any(i == chr_ends_index)) chromo = chromo + 1
    }
    chr_ends = as.integer(positions[chr_ends_index])

    slim_positions =  numeric()
    slim_rates =  numeric()
    msprime_positions =  0
    msprime_rates =  numeric()
    for (chromo in seq_len(nchr)){
      slim_positions    = c(slim_positions,    positions[which(info$Chromosome==chromo)]-1)
      slim_rates        = c(slim_rates,        rates[which(info$Chromosome==chromo)])
      msprime_positions = c(msprime_positions, positions[which(info$Chromosome==chromo)])
      msprime_rates     = c(msprime_rates,     rates[which(info$Chromosome==chromo)])
      if (chromo!=max(info$Chromosome)){
        slim_positions    = c(slim_positions,    positions[chr_ends_index[chromo]])
        slim_rates        = c(slim_rates,        0.5)
        msprime_positions = c(msprime_positions, positions[chr_ends_index[chromo]]+1)
        msprime_rates     = c(msprime_rates,     log(2))
      }
    } 
    slim_positions = as.integer(slim_positions)
    msprime_positions = as.integer(msprime_positions)
    L = as.integer(slim_positions[length(slim_positions)]) # total genome length
    return( list(nchr     = nchr,
                 chr_ends = chr_ends,
                 L        = L,
                 msprime_r_map = list(rates     = msprime_rates,
                                      positions = msprime_positions),
                 slim_r_map    = list(rates     = slim_rates,
                                      positions = slim_positions) ))
  }else{
    stop(paste("Missing columns in input file:",header_check$missing))
  }
}

### GET AGE PDF ######################################################################################
get_age_pdf = function(Sample){
  age_pdf = vector("list",Sample$size)
  if(any(Sample$is_ancient)){
    for (i in which(Sample$is_ancient)){
      calibrated = rcarbon::calibrate(Sample$age14C[i],
                                      Sample$age14Cerror[i],
                                      calCurves = Sample$calCurves[i],
                                      verbose = F)$grids[[1]]
      age_pdf[[i]] = data.frame(ageBCAD = rcarbon::BPtoBCAD(calibrated$calBP),
                                PrDens  = calibrated$PrDens)
    }
  }
  return(age_pdf)
}

### GET SAMPLE AGE INTERVAL ######################################################################################
get_sample_age_interval = function(Sample){
  oldest_sample_age = as.integer(format(Sys.Date(), "%Y"))
  if (any(!Sample$is_ancient)){
    oldest_sample_age = min(c(oldest_sample_age,Sample$ageBCAD), na.rm = TRUE)
  }
  if (any(Sample$is_ancient)){
    for (i in which(Sample$is_ancient)){
      oldest_sample_age = min(c(oldest_sample_age, Sample$age_pdf[[i]]$ageBCAD), na.rm = TRUE)
    }
  }
  youngest_sample_age = oldest_sample_age
  if (any(!Sample$is_ancient)){
    youngest_sample_age = max(c(youngest_sample_age,Sample$ageBCAD), na.rm = TRUE)
  }
  if (any(Sample$is_ancient)){
    for (i in which(Sample$is_ancient)){
      youngest_sample_age = max(c(youngest_sample_age, Sample$age_pdf[[i]]$ageBCAD), na.rm = TRUE)
    }
  }
  return(list(oldest_sample_age   = oldest_sample_age,
              youngest_sample_age = youngest_sample_age))
}
  
### SAMPLE PARAMETER TRAJECTORY THROUGH TIME ###############################################################
sample_param_trajectory = function(num_of_periods,
                                   prior_min,
                                   prior_max){
  trajectory    = array(NA,num_of_periods)
  trajectory[1] = exp(runif(1,log(prior_min),log(prior_max)))
  if (num_of_periods>1){
    for (i in 2:num_of_periods){
      alpha = runif(1,-1,1)
      trajectory[i] = 10^max(min(log10(trajectory[i-1])+alpha,log10(prior_max)),log10(prior_min))
    }
  }
  return(trajectory)
}

### GET SAMPLE AGES ######################################################################################
get_sample_ages = function(Sample, gen_len){
  sample_ages = rep(NA,Sample$size)
  for (i in seq_len(Sample$size)){
    if (!Sample$is_ancient[i]) sample_ages[i] = Sample$ageBCAD[i]
    if (Sample$is_ancient[i]){
      sample_ages[i] = sample(Sample$age_pdf[[i]]$ageBCAD, 1, prob = Sample$age_pdf[[i]]$PrDens)
    } 
  }
  sample_ages = round(abs((sample_ages-Sample$t0)/gen_len))
  return(sample_ages)
}

