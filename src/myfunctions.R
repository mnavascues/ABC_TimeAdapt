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
  if(! interactive()){
    argv <- parse_args(ap)
  }else{
    argv <- parse_args(ap, c("-q", "FALSE", # quiet
                             "-d", "1234567890" # seed
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
}
