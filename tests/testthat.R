library("gcTensor")
library("rTensor")
library("einsum")
library("testthat")

options(testthat.use_colours = FALSE)

source("testthat/setting.R")
test_file("testthat/test_InputObjectType.R")
test_file("testthat/test_OutputObjectType.R")
test_file("testthat/test_NonNegative.R")
test_file("testthat/test_Consistency.R")
test_file("testthat/test_DecreaseError.R")
# test_file("testthat/test_Estimates.R") # comment out (Heavy)
test_file("testthat/test_FixZ.R")
