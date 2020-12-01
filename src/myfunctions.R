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
  if(! interactive()){
    argv <- parse_args(ap)
  }else{
    argv <- parse_args(ap, c("-q", "FALSE",                     # quiet
                             "-d", "1234567890",                # seed
                             "-p", "test",                      # project_name
                             "-b", "1",                         # batch_ID
                             "-i", "data/sample_info_test.txt", # sample_info_file
                             "-g", "data/genome_info_test.txt"  # genome_info_file
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
  write(paste0("Sample file: ",argv$sample_info_file), stdout())
  write(paste0("Genome file: ",argv$genome_info_file), stdout())
}
