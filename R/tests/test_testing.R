library(testthat)
context("Testing csr matrix generation")

source("../stan_helpers.R")

# test case for division of two numbers
test_that("get_sparse", {
  expect_equal(4/2, 2)
})


