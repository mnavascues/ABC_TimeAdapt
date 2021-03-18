library(argparser, quietly=TRUE)
source("../../R/myfunctions.R")

context("Command line arguments")

passed <- TRUE
send_error <- FALSE

# get_arguments()
passed <- test_that("Check for incomplete set of arguments", {
  source("../testarg.R")
  argv <- get_arguments(test=TRUE, test_arg=testarg)

  expect_equal(argv$help,F)
  expect_equal(argv$quiet,T)
  expect_equal(argv$seed,12345678909)
  expect_equal(argv$project_name,"test")
  expect_equal(argv$batch_ID,1)
  expect_equal(argv$sample_info_file,"tests/sample_info_test.txt")
  expect_equal(argv$genome_info_file,"tests/genome_info_test.txt")
  expect_equal(argv$num_of_gen_in_forw_sim,400)
  expect_equal(argv$num_of_periods_forw,8)
  expect_equal(argv$generation_length_prior_params,c(2.000000,1.465967,26.000000,30.000000))
  expect_equal(argv$num_of_sims,3)
  expect_equal(argv$population_size_prior_params,c(10,200))
  expect_equal(argv$mutation_rate_prior_params,c(5e-08,5e-01))

  expect_error(get_arguments(test=TRUE, test_arg=testarg_error))
})
if (!passed) send_error <- TRUE
# print_arguments()

  
context("Read info files")

# check_file_header()
passed <- test_that("Check header sends error when missing coloumn", {
  expect_error(check_file_header(c("a","b"),c("a")),"Missing columns in input file: b")
  expect_error(check_file_header(c("a"),c("b")),"Missing columns in input file: a")
  expect_error(check_file_header(c("a"),character()),"Missing columns in input file: a")
  expect_equal(check_file_header(c("a","b"),c("a","b","c")),T)
  expect_equal(check_file_header(c("a","b"),c("a","b")),T)
})
if (!passed) send_error <- TRUE
# read_sample_info("../sample_info_test.txt")
passed <- test_that("Readings sample info file sends error when wrong data type", {
  fname <- tempfile()
  write("sampleID age14C age14Cerror year coverage damageRepair groups", fname)
  write("1 a a a a a a", fname, append=T)
  expect_error(read_sample_info(fname),paste("Wrong data type in file:",fname))
  unlink(fname)
})
if (!passed) send_error <- TRUE

# read_genome_info()

context("Ages of samples")

# get_sample_cal_age_PDF()
# check_ts_lower_gen_in_for_sim()

context("Sample from priors")

# sample_N_trajectory()
# sample_demography_from_prior()
# sample_ages_from_prior()

if (send_error) stop('Some tests did not pass')