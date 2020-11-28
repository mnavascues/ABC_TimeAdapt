# TimeAdapt
# Miguel Navascués
# Uppsala universitet & INRAE
# 2020

get_arguments <- function(){
  ap <- arg_parser(description = paste0("Approximate Bayesian computation analysis ",
                                        "for joint inference of demography and selection ",
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
  if(! interactive()){
    argv <- parse_args(ap)
  }else{
    argv <- parse_args(ap, c("-q", "FALSE",      # quiet
                             "-d", "1234567890", # seed
                             "-p", "test",       # project
                             "-b", "1"           # batch_ID
                             ))
  }
  if (!argv$quiet){
    print(ap)
    print_arguments(argv)
  }
  return(argv)
}

print_arguments <- function(argv){
  write(paste0("Quiet: ",argv$quiet), stdout())
  write(paste0("Seed: ",argv$seed), stdout())
  write(paste0("Project: ",argv$project_name), stdout())
  write(paste0("Batch: ",argv$batch_ID), stdout())
}
