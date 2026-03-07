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

testthat("Check that the bcgd function works properly on a minimal example.", {
  # Set variables
  x_test = c(1,1,2,2)
  n_samples_test = c(2,2)
  n_total_test = length(x_test)
  m_test = length(n_samples_test) - 1
  d_test = 1
  model_test = 1
  lambda_test = 0.5
  max_iters_test = 1
  theta_test_0 = c(0,0)
  
  # Do some preliminary tests to validate inputs
  expect_equal(
    # Validate LDL
    negLDL(theta_test_0, x_test, n_total_test, n_samples_test, m_test, model_test, d_test),
    0
  )
  expect_equal(
    # Validate Gradient
    negLDLGr(theta_test_0, x_test, n_total_test, n_samples_test, m_test, model_test, d_test),
    c(0,-1)
  )
  expect_equal(
    # Validate Hessian
    diag(negLDLHessian(theta_test_0, x_test, n_total_test, n_samples_test, m_test, model_test, d_test)),
    c(1,2.5)
  )
  
  # In the first iteration, the intercept will not be moved as d_g = 0
  # In the second iteration, d_g = 1/5
  # Validate some intermediary inputs
  expect_equal(
    # Validate objective function value
    negLDLGL(theta_test_0 + c(0,1/5), x_test, n_total_test, n_samples_test, m_test, model_test, d_test, lambda_test),
    (-7/10)+2*(log((1/4)*(1+exp(1/5))*(1+exp(2/5))))
  )
  
  # Now we should immediately satisfy the amijo rule
  # Thus, we will move 1/5 for beta 1, and this should be the output
  # Now, validate the full algorithm for one iteration
  expected_bcgd_output = list(obj=negLDLGL(theta_test_0 + c(0,1/5), x_test, n_total_test, n_samples_test, m_test, model_test, d_test, lambda_test), iters=1, par=c(0,1/5))
  expect_equal(
    bcgd(
      theta_0=theta_test_0, x=x_test, n_total=n_total_test, n_samples=n_samples_test,
      m=m_test, d=d_test, model=model_test, lambda=lambda_test, max_iters=max_iters_test
    ),
    expected_bcgd_output
  )
  
  # So the objective function moved from 0 to -0.05028. So, if we increase iterations, but change the tol to 0.1, we should return the same value
  expect_equal(
    bcgd(
      theta_0=theta_test_0, x=x_test, n_total=n_total_test, n_samples=n_samples_test,
      m=m_test, d=d_test, model=model_test, lambda=lambda_test, max_iters=100, threshold=0.1
    ),
    expected_bcgd_output
  )
  
  # But if we drop this back down to 0.01, we should get a different answer
  expect_false(
    bcgd(
      theta_0=theta_test_0, x=x_test, n_total=n_total_test, n_samples=n_samples_test,
      m=m_test, d=d_test, model=model_test, lambda=lambda_test, max_iters=100, threshold=0.01
    ) == expected_bcgd_output
  )
})


# Test bcgd in a slightly more complex context
# Case 1:
# Basis function = (x, x^2)
# m = 2
# theta_0 = (0, 0, 0, 0, 0, 0)
# max_iters = 1
# x_0 = (1,1)
# x_1 = (2,2)
# x_2 = (3,3)
# lambda = 1

testthat("Check a slightly more complex version for one iteration.", {
  # Set variables
  x_test = c(1, 1, 2, 2, 3, 3)
  n_samples_test = c(2, 2, 2)
  n_total_test = length(x_test)
  m_test = length(n_samples_test) - 1
  d_test = 2
  model_test = 5
  lambda_test = 1
  max_iters_test = 1
  theta_test_0 = c(0, 0, 0, 0, 0, 0)
  omega_test = 1
  psi_test = 0.5
  sigma_test = 0.1
  # Setup indices for groups
  idx0 = c(1,4)
  idx1 = c(2,5)
  idx2 = c(3,6)
  # Setup function for finding w
  weight_finder = function(theta, d_g, x, n_total, n_samples, m, model, d, lambda, omega, psi, sigma, obj, delta){
    w = omega
    while(negLDLGL(theta + w*d_g, x, n_total, n_samples, m, model, d, lambda) 
          > 
          (obj + w*sigma*delta)){
      w = w*psi
    }
    return(w)
  }
  
  # Infer what the algorithm should be doing
  
  # When g == 0 (intercept)
  obj0 = negLDLGL(theta_test_0, x_test, n_total_test, n_samples_test, m_test, model_test, d_test, lambda_test)
  grad0 = negLDLGr(theta_test_0, x_test, n_total_test, n_samples_test, m_test, model_test, d_test)[idx0]
  hess0 = max(diag(negLDLHessian(theta_0_test_0, x_test, n_total_test, n_samples_test, m_test, model_test, d_test))[idx0])
  d_g_0 = (1/hess0)*(-1)*grad0 # basically 0s
  delta_0 = sum(grad0*d_g_0) + negLDLGL(theta_test_0 + d_g_0, x_test, n_total_test, n_samples_test, m_test, model_test, d_test, lambda_test) - obj0 
  w_0 = weight_finder(theta_test_0, d_g_0, x_test, n_total_test, n_samples_test, m_test, model_test, d_test, lambda_test, omega_test, psi_test, sigma_test, obj0, delta_0)
  # Update theta - basically 0
  theta_test_1 = theta_test_0
  theta_test_1[idx0] = theta_test_1[idx0] + w_0*d_g_0
  
  # When g == 1 (x coefficients)
  obj1 = negLDLGL(theta_test_1, x_test, n_total_test, n_samples_test, m_test, model_test, d_test, lambda_test)
  grad1 = negLDLGr(theta_test_1, x_test, n_total_test, n_samples_test, m_test, model_test, d_test)[idx1]
  hess1 = max(diag(negLDLHessian(theta_test_1, x_test, n_total_test, n_samples_test, m_test, model_test, d_test))[idx1])
  if (norm((-1)*grad1+hess1*theta_test_1[idx1], type="2") > lambda_test){
    d_g_1 = (1/hess1)*((-1)*grad1 - lambda_test*((hess1*theta_test_1[idx1] - grad1)/norm((-1)*grad1+hess1*theta_test_1[idx1], type="2")))
  } else {
    d_g_1 = (-1)*theta_test_1[idx1]
  }
  delta_1 = sum(grad1*d_g_1) + negLDLGL(theta_test_1 + d_g_1, x_test, n_total_test, n_samples_test, m_test, model_test, d_test, lambda_test) - obj1 
  w_1 = weight_finder(theta_test_1, d_g_1, x_test, n_total_test, n_samples_test, m_test, model_test, d_test, lambda_test, omega_test, psi_test, sigma_test, obj1, delta_1)
  # Update theta - moving slightly now
  theta_test_2 = theta_test_1
  theta_test_2[idx1] = theta_test_2[idx1] + w_1*d_g_1
  
  # When g == 2 (x^2 coefficients)
  obj2 = negLDLGL(theta_test_2, x_test, n_total_test, n_samples_test, m_test, model_test, d_test, lambda_test)
  grad2 = negLDLGr(theta_test_2, x_test, n_total_test, n_samples_test, m_test, model_test, d_test)[idx2]
  hess2 = max(diag(negLDLHessian(theta_test_2, x_test, n_total_test, n_samples_test, m_test, model_test, d_test))[idx2])
  if (norm((-1)*grad2+hess2*theta_test_2[idx2], type="2") > lambda_test){
    d_g_2 = (1/hess2)*((-1)*grad2 - lambda_test*((hess2*theta_test_2[idx2] - grad2)/norm((-1)*grad2+hess2*theta_test_2[idx2], type="2")))
  } else {
    d_g_2 = (-1)*theta_test_2[idx2]
  }
  delta_2 = sum(grad2*d_g_2) + negLDLGL(theta_test_2 + d_g_2, x_test, n_total_test, n_samples_test, m_test, model_test, d_test, lambda_test) - obj2 
  w_2 = weight_finder(theta_test_2, d_g_2, x_test, n_total_test, n_samples_test, m_test, model_test, d_test, lambda_test, omega_test, psi_test, sigma_test, obj2, delta_2)
  # Update theta - moving slightly now
  theta_test_final = theta_test_2
  theta_test_final[idx2] = theta_test_final[idx2] + w_2*d_g_2
  obj_final = negLDLGL(theta_test_final, x_test, n_total_test, n_samples_test, m_test, model_test, d_test, lambda_test)
  
  # Now validate the results
  expect_equal(
    # Validate Gradient
    bcgd(
      theta_0=theta_test_0, x=x_test, n_total=n_total_test, n_samples=n_samples_test,
      m=m_test, d=d_test, model=model_test, lambda=lambda_test, max_iters=max_iters_test
    ),
    list(obj = obj_final, iters = 1, par = theta_test_final)
  )
})

