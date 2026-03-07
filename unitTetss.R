# Test functions to ensure accuracy

# Load library
library(testthat)

# Load functions
source("R\\init.R")

# Test bcgd
# Case 1:
# Basis function = x
# m = 1
# theta_0 = (0,0)
# max_iters = 1
# x_0 = (1,1)
# x_1 = (2,2)
# lambda = 1

testthat("Check that the bcgd function works properly on some minimal examples", {
  # Set variables
  x_test = c(1,1,2,2)
  n_samples_test = c(2,2)
  n_total_test = length(n_samples_test)
  n_total = sum(n_samples_test)
  m_test = length(n_samples_test) - 1
  d_test = 1
  model_test = 1
  lambda_test = 1
  max_iters_test = 1
  theta_0_test_0 = c(0,0)
  theta_0_test_1 = c(0,1)
  
  # Do some preliminary tests to validate inputs
  expect_equal(
    # Validate LDL
    negLDL(theta_0_test_0, x_test, n_total_test, n_samples_test, m_test, model_test, d_test),
    0
  )
  expect_equal(
    # Validate Gradient
    negLDLGr(theta_0_test_0, x_test, n_total_test, n_samples_test, m_test, model_test, d_test),
    c(0,-1)
  )
  expect_equal(
    # Validate Hessian
    diag(negLDLHessian(theta_0_test_0, x_test, n_total_test, n_samples_test, m_test, model_test, d_test)),
    c(1,2.5)
  )
  
  # Now, validate the full algorithm for one iteration
  expect_equal(
      bcgd(
        theta_0=theta_0_test_0, x=x_test, n_total=n_total_test, n_samples=n_samples_test,
        m=m_test, d=d_test, model=model_test, lambda=lambda_test, max_iters=max_iters_test
      )$par,
      c(0,0)
  )
})