source("../../R/myfunctions.R")

context("Command line arguments")

# get_arguments()
# print_arguments()

context("Read info files")

# check_file_header()
test_that("Check header sends error when missing coloumn", {
  expect_error(check_file_header(c("a","b"),c("a")))
  expect_error(check_file_header(c("c"),c("a")))
  expect_error(check_file_header(c("c"),character()))
})
# read_sample_info("../../data/sample_info_test.txt")
# read_genome_info()

context("Ages of samples")

# get_sample_cal_age_PDF()
# check_ts_lower_gen_in_for_sim()

context("Sample from priors")

# sample_N_trajectory()
# sample_demography_from_prior()

