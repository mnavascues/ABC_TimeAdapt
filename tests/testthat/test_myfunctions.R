source("../../R/myfunctions.R")

context("Read info files")

test_that("Verify type in sample info file", {
  info <- read_sample_info("../../data/sample_info_test.txt")
  expect_true(is.character(info$id))
  expect_true(is.numeric(info$age14C))
  expect_true(is.numeric(info$age14Cerror))
  expect_true(is.numeric(info$ageBCAD))
  expect_true(is.numeric(info$coverage))
  expect_true(is.logical(info$is_dr))
  expect_true(is.numeric(info$total_ancient))
  expect_true(is.numeric(info$size))
  expect_true(is.numeric(info$t0))
})
})
