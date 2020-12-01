source("../../R/myfunctions.R")

context("Read info files")

test_that("Verify type in sample info file", {
  expect_true(is.character(read_sample_info("../../data/sample_info_test.txt")$id))
  expect_true(is.numeric(read_sample_info("../../data/sample_info_test.txt")$age14C))
  expect_true(is.numeric(read_sample_info("../../data/sample_info_test.txt")$age14Cerror))
  expect_true(is.numeric(read_sample_info("../../data/sample_info_test.txt")$ageBCAD))
  expect_true(is.numeric(read_sample_info("../../data/sample_info_test.txt")$coverage))
  # expect_true(is.logical(read_sample_info()$is_modern))
  # expect_true(is.logical(read_sample_info()$is_ancient))
  # expect_true(is.logical(read_sample_info()$is_dr))
  expect_true(is.numeric(read_sample_info("../../data/sample_info_test.txt")$total_ancient))
  expect_true(is.numeric(read_sample_info("../../data/sample_info_test.txt")$size))
  expect_true(is.numeric(read_sample_info("../../data/sample_info_test.txt")$t0))
})
