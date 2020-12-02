source("../../R/myfunctions.R")

context("Read info files")

# read_sample_info("../../data/sample_info_test.txt")

test_that("Check header sends error when missing coloumn", {
  expect_error(check_file_header(c("a","b"),c("a")))
  expect_error(check_file_header(c("c"),c("a")))
  expect_error(check_file_header(c("c"),character()))
})
