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
source("timeadapt.R")

context("Read info files")

# check_file_header()
test_that("Check header sends error when missing coloumn", {
  expect_error(check_file_header(c("a","b"),c("a")),"Missing columns in input file: b")
  expect_error(check_file_header(c("a"),c("b")),"Missing columns in input file: a")
  expect_error(check_file_header(c("a"),character()),"Missing columns in input file: a")
  expect_true(check_file_header(c("a","b"),c("a","b","c")))
  expect_true(check_file_header(c("a","b"),c("a","b")))
})

# read_sample_info("../sample_info_test.txt")
test_that("Reading sample info file sends error when wrong data type", {
  fname <- tempfile()
  write("sampleID age14C age14Cerror year coverage damageRepair groups", fname)
  write("1 a a a a a a", fname, append=T)
  expect_error(read_sample_info(fname),paste("Wrong data type in file:",fname))
  unlink(fname)
})
test_that("Reading sample info file gets right values", {
  fname <- tempfile()
  write("sampleID age14C age14Cerror year coverage damageRepair groups", fname)
  write("modern   NA     NA          2010 30.03    TRUE         0", fname, append=T)
  write("ancient  1980   20          NA   10.01    TRUE         1", fname, append=T)
  Sample <- read_sample_info(fname)
  expect_equal(Sample$id,            c("modern", "ancient"))
  expect_equal(Sample$age14C,        c(NA, 1980))
  expect_equal(Sample$age14Cerror,   c(NA, 20))
  expect_equal(Sample$ageBCAD,       c(2010, NA))
  expect_equal(Sample$coverage,      c(30.03, 10.01))
  expect_equal(Sample$is_modern,     c(T,F))
  expect_equal(Sample$is_ancient,    c(F,T))
  expect_equal(Sample$is_dr,         c(T,T))
  expect_equal(Sample$total_ancient, 1)
  expect_equal(Sample$size,          2)
  expect_equal(Sample$t0,            2010)
  unlink(fname)
})

# read_genome_info()
test_that("Reading genome info file gets right values", {
  fname <- tempfile()
  write("Chromosome Position Recombination_rate", fname)
  write("1          500000	 1.11E-09", fname, append=T)
  write("2          1000000  1.01E-09", fname, append=T)
  Genome <- read_genome_info(fname)
  expect_equal(Genome$number_of_chromosomes, 2)
  expect_equal(Genome$L,                     1500000-1)
  expect_equal(Genome$rec_map_SLiM[,1],      c(500000-1, 500000, 1500000-1))
  expect_equal(Genome$rec_map_SLiM[,2],      c(1.11E-09, 0.5, 1.01E-09))
  expect_true(is.integer(Genome$number_of_chromosomes))
  expect_true(is.integer(Genome$L))
  expect_true(is.integer(Genome$rec_map_SLiM$ends))
  expect_true(is.double(Genome$rec_map_SLiM$rates))
  unlink(fname)
})




context("Ages of samples")

# get_sample_cal_age_PDF()
# maximum_age_of_sample()

test_that("Gets correct maximum age", {
  test_cal_age_PDF <- list(NULL,data.frame(calBP=c(2001,2000,1999),PrDens=c(1.6-05,1.6-05,1.6e-05)))
  test_Sample <- list(is_ancient=c(F,T),is_modern=c(T,F),ageBCAD=c(2020,NA),t0=2020)
  expect_equal(maximum_age_of_sample(test_Sample,test_cal_age_PDF,25), 83)
})



context("Sample from priors")

# sample_N_trajectory()
test_that("Output is correct size", {
  a=100;b=1;c=10
  test_N_trajectory <- sample_N_trajectory(a,b,c)
  expect_length(test_N_trajectory,a)
  expect_true(all(test_N_trajectory>=b))
  expect_true(all(test_N_trajectory<=c))
})

# sample_demography_from_prior()
# sample_ages_from_prior()

