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
    f_argv <- parse_args(ap)
  }else{
    f_argv <- parse_args(ap, c("-q", "FALSE",                     # quiet
                             "-d", "1234567890",                # seed
                             "-p", "test",                      # project_name
                             "-b", "1",                         # batch_ID
                             "-i", "data/sample_info_test.txt", # sample_info_file
                             "-g", "data/genome_info_test.txt"  # genome_info_file
                             ))
  }
  if (!f_argv$quiet){
    print(ap)
    print_arguments(f_argv)
  }
  return(f_argv)
}

print_arguments <- function(f_argv){
  write(paste0("Quiet: ",f_argv$quiet), stdout())
  write(paste0("Seed: ",f_argv$seed), stdout())
  write(paste0("Project: ",f_argv$project_name), stdout())
  write(paste0("Batch: ",f_argv$batch_ID), stdout())
  write(paste0("Sample file: ",f_argv$sample_info_file), stdout())
  write(paste0("Genome file: ",f_argv$genome_info_file), stdout())
}


verify_file_header <- function(expected_header,file_header){
  missing=!is.element(expected_header,file_header)
  if (any(missing)){
    stop(paste("Missing columns in input file:",expected_header[missing]))
  } 
}

read_sample_info <- function(file="data/sample_info_test.txt"){
  info <- read.table(file,header=T,stringsAsFactors=F,strip.white=T)
  expected_header <- c("sampleID","age14C","age14Cerror","year","coverage","damageRepair","groups")
  verify_file_header(expected_header, file_header = colnames(info))
  
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

