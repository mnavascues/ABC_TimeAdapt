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

library(ini, quietly = TRUE)
source("scripts/timeadapt.R")

# read command line arguments
args = commandArgs(trailingOnly = TRUE)
if (length(args) == 0){
  stop("One positional argument necessary in the command (options file)")
  quit(save="no")
} else options_file = args[1] # options_file = "tests/config_project.ini"

# read options file
opts     = get_project_options(options_file)
Settings = opts$Settings
Model    = opts$Model
Priors   = opts$Priors

# print script info to screen
print_info("setproject.R", Settings$verbose, project = Settings$project)

# set seed for random number generator
set.seed(Settings$seed)
#batch_seeds = as.integer(round(runif(Settings$num_of_batches, 0, 2^31)))
#if (Settings$verbose >= 10) write(paste("Seed batch", seq_len(Settings$num_of_batches), ":", batch_seeds), stdout())

# create results directory
dir.create("results", showWarnings = FALSE)
project_dir <- paste("results", Settings$project, sep="/")
dir.create(project_dir, showWarnings = FALSE)
if (Settings$verbose >= 10) write(paste("Project folder:", project_dir), stdout())

# read sample and genome information from tables in text files
Sample <- read_sample_info(Settings$sample_file)
if (Settings$verbose >= 10) write(paste("Sample size:", Sample$size), stdout())  
Genome <- read_genome_info(Settings$genome_file)
if (Settings$verbose >= 10) write(paste("Genome length:", Genome$L+1), stdout()) 

# calculate times of change
times_of_change_forw = get_times_of_change(Model$generations_forward,
                                           Model$periods_forward)
times_of_change_back = get_times_of_change(Model$generations_coalescence,
                                           Model$periods_coalescence,
                                           mode="exponential")
if (Settings$verbose >= 10){
  write("Times of change:", stdout())
  write(times_of_change_forw, stdout(), ncolumns = length(times_of_change_forw))
  write(times_of_change_back, stdout(), ncolumns = length(times_of_change_back))
}

# calculate probability distribution curves for calibrated age of ancient samples
# (from 14C ages)
age_pdf = get_age_pdf(Sample)

# get the oldest possible age of samples and verify that falls within forward time simulation period
sample_age_interval = get_sample_age_interval(Sample, age_pdf)
forward_simulation_rage = 2 + round(abs(sample_age_interval$oldest_sample_age - sample_age_interval$youngest_sample_age) / Priors$gen_len_min)
if (forward_simulation_rage >= Model$generations_forward){
  stop("Larger period in forward simulation necessary to include oldest sample")
  quit(save="no")
}
if (Settings$verbose >= 10) write(paste("Youngest sample age:", sample_age_interval$youngest_sample_age), stdout()) 
if (Settings$verbose >= 10) write(paste("Oldest sample age:", sample_age_interval$oldest_sample_age), stdout()) 

# groups
group_levels = nchar(Sample$groups[1])


# save options list in RDS file (for R)
Settings = modifyList(Settings, list(project_dir          = project_dir))
Model    = modifyList(Model,    list(times_of_change_forw = times_of_change_forw,
                                     times_of_change_back = times_of_change_back))
Sample   =  modifyList(Sample,  list(age_pdf              = age_pdf,
                                     t0                   = sample_age_interval$youngest_sample_age))
opts = list(Settings = Settings,
            Model    = Model,
            Priors   = Priors,
            Sample   = Sample,
            Genome   = Genome)  
saveRDS(opts, file = paste0(project_dir,"/project_options.RDS"))

# project options for Pyhton (.ini file)
project_ini <- list()
project_ini[["Settings"]] = list(project     = Settings$project,
                                 verbose     = Settings$verbose,
                                 project_dir = project_dir)
project_ini[["Model"]]    = list(periods_coalescence  = Model$periods_coalescence,
                                 times_of_change_back = paste(as.integer(times_of_change_back), collapse=" "))
project_ini[["Sample"]]   = list(size         = Sample$size,
                                 coverage     = paste(Sample$coverage, collapse=" "),
                                 is_damaged   = paste(Sample$is_damaged, collapse=" "),
                                 group_levels = group_levels,
                                 groups       = paste(Sample$groups, collapse=" "))
project_ini[["Genome"]]   = list(nchr                    = Genome$nchr,
                                 chr_ends                = paste(Genome$chr_ends, collapse=" "),
                                 msprime_r_map_positions = paste(Genome$msprime_r_map$positions, collapse=" "),
                                 msprime_r_map_rates     = paste(Genome$msprime_r_map$rates, collapse=" "))



write.ini(project_ini, file = paste0(project_dir,"/project_options.ini"))
