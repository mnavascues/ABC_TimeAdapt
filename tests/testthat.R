#    TimeAdapt: joint inference of demography and selection
#    Copyright (C) 2021  Miguel de Navascués, Uppsala universitet, INRAE
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <https://www.gnu.org/licenses/>.
#
source("../scripts/timeadapt.R")
#source("scripts/timeadapt.R")

context("Read info files")

# check_file_header()
test_that("check_file_header() returns FALSE and missing elements if given incomplete header", {
  expect_false(check_file_header(c("a","b"),"a")$ok)
  expect_equal(check_file_header(c("a","b"),"a")$missing,"b")
  expect_false(check_file_header("a","b")$ok)
  expect_equal(check_file_header("a","b")$missing,"a")
  expect_false(check_file_header("a",character())$ok)
  expect_equal(check_file_header("a",character())$missing,"a")
})
test_that("check_file_header() returns TRUE when all expected elements are present, additional elements are OK", {
  expect_true(check_file_header(c("a","b"),c("a","b","c"))$ok)
  expect_true(check_file_header(c("a","b"),c("a","b"))$ok)
})

# read_sample_info()
test_that("read_sample_info() returns an error when wrong data type: coverage should be numeric", {
  fname <- tempfile()
  write("sampleID age14C age14Cerror calCurves year coverage damaged groups", fname)
  write("1 NA NA NA NA a TRUE a", fname, append=T)
  expect_error(read_sample_info(fname),paste("Wrong data type in file:",fname))
  unlink(fname)
})
test_that("read_sample_info() returns an error when wrong data type: groups should be character", {
  fname <- tempfile()
  write("sampleID age14C age14Cerror calCurves year coverage damaged groups", fname)
  write("1 NA NA NA NA 40 TRUE 1", fname, append=T)
  expect_error(read_sample_info(fname),paste("Wrong data type in file:",fname))
  unlink(fname)
})
test_that("read_sample_info() returns an error when wrong data type: damaged should be logical with no missing data", {
  fname <- tempfile()
  write("sampleID age14C age14Cerror calCurves year coverage damaged groups", fname)
  write("1 1000 10 intcal20 NA 40 NA a", fname, append=T)
  expect_error(read_sample_info(fname),paste("Wrong data type in file:",fname))
  unlink(fname)
})
test_that("read_sample_info() returns an error when missing columns", {
  fname <- tempfile()
  write("sampleID age14C age14Cerror calCurves year damaged groups", fname)
  write("1 1000 10 intcal20 NA TRUE a", fname, append=T)
  expect_error(read_sample_info(fname),paste("Missing columns in input file:","coverage"))
  unlink(fname)
})
test_that("read_sample_info() gets right values for samples", {
  Sample <- read_sample_info("test_sample.txt")
  expect_equal(c(Sample$total_ancient,Sample$size), c(13,17))
})

# read_genome_info()
test_that("read_genome_info() gets right values for genome", {
  Genome <- read_genome_info("test_genome.txt")
  expect_true(is.integer(Genome$L))
  expect_true(is.integer(Genome$nchr))
  expect_true(is.integer(Genome$slim_r_map$positions))
  expect_true(is.double(Genome$slim_r_map$rates))
  expect_true(is.integer(Genome$msprime_r_map$positions))
  expect_true(is.double(Genome$msprime_r_map$rates))
})




