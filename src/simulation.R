#!/usr/bin/env Rscript

# TimeAdapt
# Miguel Navascués
# Uppsala universitet & INRAE
# 2020

library(argparser, quietly=TRUE)
source("src/myfunctions.R")

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

