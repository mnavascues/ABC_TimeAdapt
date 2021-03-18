library(testthat)
res <- test_file("tests/testthat/test_myfunctions.R")
failed <- sum(print(res)$failed)
if (failed > 0) stop("Failed test(s) present")