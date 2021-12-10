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

source("scripts/timeadapt.R")

# read command line arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0){
  stop("Two positional arguments necessary in the command (RDS options file and batch number)")
  quit(save="no")
} else {
  options_file = args[1]             # options_file = "results/test/project_options.RDS"
  batch        = as.integer(args[2]) # batch = 1
}

opts = readRDS(options_file)
Settings = opts$Settings
Model    = opts$Model
Priors   = opts$Priors
Sample   = opts$Sample
Genome   = opts$Genome

# print script info to screen
print_info("reftable.R", Settings$verbose, project = Settings$project, batch = batch)
if (Settings$verbose>10) quiet = FALSE else quiet = TRUE

# set empty data frames for sumstats and latent variables
sumstats_header = scan(paste0(Settings$project_dir,"/",batch,"/sumstats_1.txt"),
                       what = character(), nlines = 1, quiet = quiet)
ref_table_sumstats = data.frame(matrix(nrow=Settings$num_of_sims,
                                       ncol=length(sumstats_header)))
names(ref_table_sumstats) = sumstats_header

latent_variables_header = scan(paste0(Settings$project_dir,"/",batch,"/latent_variables_1.txt"),
                               what = character(), nlines = 1, quiet = quiet)
ref_table_latent_variables = data.frame(matrix(nrow=Settings$num_of_sims,
                                               ncol=length(latent_variables_header)))
names(ref_table_latent_variables) = latent_variables_header

for (sim in seq_len(Settings$num_of_sims)){
  ref_table_latent_variables[sim,] = read.table(paste0(Settings$project_dir,"/",batch,"/latent_variables_",sim,".txt"), header = T)
  ref_table_sumstats[sim,] = read.table(paste0(Settings$project_dir,"/",batch,"/sumstats_",sim,".txt"), header = T)
}

saveRDS(ref_table_latent_variables, paste0(Settings$project_dir,"/",batch,"/ref_table_latent_variables.RDS"))
saveRDS(ref_table_sumstats, paste0(Settings$project_dir,"/",batch,"/ref_table_sumstats.RDS"))
                