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

test_that("Check that the bcgd function works properly on a minimal example.", {
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
      m=m_test, model=model_test, d=d_test, lambda=lambda_test, max_iters=max_iters_test
    ),
    expected_bcgd_output
  )
  
  # So the objective function moved from 0 to -0.05028. So, if we increase iterations, but change the tol to 0.1, we should return the same value
  expect_equal(
    bcgd(
      theta_0=theta_test_0, x=x_test, n_total=n_total_test, n_samples=n_samples_test,
      m=m_test, model=model_test, d=d_test, lambda=lambda_test, max_iters=100, threshold=0.1
    ),
    expected_bcgd_output
  )
  
  # But if we drop this back down to 0.01, we should get a different answer
  expect_false(
    identical(bcgd(
      theta_0=theta_test_0, x=x_test, n_total=n_total_test, n_samples=n_samples_test,
      m=m_test, model=model_test, d=d_test, lambda=lambda_test, max_iters=100, threshold=0.01
    ), expected_bcgd_output)
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

test_that("Check a slightly more complex version for one iteration.", {
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
  hess0 = max(diag(negLDLHessian(theta_test_0, x_test, n_total_test, n_samples_test, m_test, model_test, d_test))[idx0])
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
      m=m_test, model=model_test, d=d_test, lambda=lambda_test, max_iters=max_iters_test
    ),
    list(obj = obj_final, iters = 1, par = theta_test_final)
  )
})


test_that("Test the AIC/BIC function", {
  x_test = c(1, 1, 2, 2, 3, 3)
  n_samples_test = c(2, 2, 2)
  n_total_test = length(x_test)
  m_test = length(n_samples_test) - 1
  d_test = 2
  model_test = 5
  
  theta_test_1 = c(0, 0, 0, 0, 0, 0)
  # This should just return 2*negLDL for both because there will be no parameter for the penalties
  expect_equal(
    aic_bic_drm(theta_test_1, x_test, n_total_test, n_samples_test, m_test, model_test, d_test),
    rep(2*negLDL(theta_test_1, x_test, n_total_test, n_samples_test, m_test, model_test, d_test), 2)
  )
  
  theta_test_2 = c(1, 0, 0, 1, 0, 0)
  # These should also equal because we don't penalize the intercept
  expect_equal(
    aic_bic_drm(theta_test_2, x_test, n_total_test, n_samples_test, m_test, model_test, d_test),
    rep(2*negLDL(theta_test_2, x_test, n_total_test, n_samples_test, m_test, model_test, d_test), 2)
  )
  
  fn = function(x){
    return(c(x,x^2))
  }
  # These should return the same, the function will just use the user function in the second case
  expect_equal(
    aic_bic_drm(theta_test_2, x_test, n_total_test, n_samples_test, m_test, model_test, d_test),
    aic_bic_drm(theta_test_2, x_test, n_total_test, n_samples_test, m_test, fn, d_test)
  )
  
  theta_test_3 = c(1, 1, 0, 1, 1, 0)
  ldlVal3 = negLDL(theta_test_3, x_test, n_total_test, n_samples_test, m_test, model_test, d_test)
  # we should be adding 2*2 = 4 to AIC, 2*log(n) to BIC
  expect_equal(
    aic_bic_drm(theta_test_3, x_test, n_total_test, n_samples_test, m_test, model_test, d_test),
    c(2*ldlVal3 + 4, 2*ldlVal3 + 2*log(n_total_test))
  )
  
  # If the values are less than the set tolerance, we should again get the LDL values
  theta_test_4 = c(1, 1e-13, 0, 1, 1e-13, 0)
  # we should be adding 2*2 = 4 to AIC, 2*log(n) to BIC
  expect_equal(
    aic_bic_drm(theta_test_4, x_test, n_total_test, n_samples_test, m_test, model_test, d_test),
    rep(2*negLDL(theta_test_4, x_test, n_total_test, n_samples_test, m_test, model_test, d_test), 2)
  )
})


test_that("Test the SolutionPath function", {
  # Set variables
  x_test = c(rep(1,5), rep(2,5), rep(3,5))
  n_samples_test = rep(5,3)
  n_total_test = length(x_test)
  m_test = length(n_samples_test) - 1
  d_test = 5
  model_test = 12
  
  # Confirm that function returns -1 if we pass a negative lambda value
  lambda_vals_0 = c(-1,1)
  expect_equal(
    solutionPath(x_test, n_total_test, n_samples_test, m_test, model_test, d_test, lambdaVals = lambda_vals_0),
    -1
  )
  
  # Confirm that if we don't pass 0, we get the same result as passing 0
  lambda_vals_1_a = c(0,1)
  lambda_vals_1_b = c(1)
  expect_equal(
    solutionPath(x_test, n_total_test, n_samples_test, m_test, model_test, d_test, lambdaVals = lambda_vals_1_a),
    solutionPath(x_test, n_total_test, n_samples_test, m_test, model_test, d_test, lambdaVals = lambda_vals_1_b)
  )
  
  # Actually confirm what the values should be
  lambda_vals_2 = c(1,2)
  drm = drmdel(x=x_test, n_samples=n_samples_test, basis_func=model_test)
  aic_bic_mele = aic_bic_drm(theta = drm$mele, x=x_test, n_total=n_total_test, n_samples=n_samples_test, m=m_test, basis_func=model_test, d=d_test)
  bcgd1 = bcgd(
    theta_0=drm$mele, x=x_test, n_total=n_total_test, n_samples=n_samples_test, m=m_test, model=model_test, d=d_test, lambda=lambda_vals_2[1]
  )
  aic_bic_1 = aic_bic_drm(theta = bcgd1$par, x=x_test, n_total=n_total_test, n_samples=n_samples_test, m=m_test, basis_func=model_test, d=d_test)
  bcgd2 = bcgd(
    theta_0=bcgd1$par, x=x_test, n_total=n_total_test, n_samples=n_samples_test, m=m_test, model=model_test, d=d_test, lambda=lambda_vals_2[2]
  )
  aic_bic_2 = aic_bic_drm(theta = bcgd2$par, x=x_test, n_total=n_total_test, n_samples=n_samples_test, m=m_test, basis_func=model_test, d=d_test)
  expected_return_matrix = matrix(
    c(0, 0, drm$negldl, aic_bic_mele, drm$mele,
      lambda_vals_2[1], bcgd1$iters, bcgd1$obj, aic_bic_1, bcgd1$par,
      lambda_vals_2[2], bcgd2$iters, bcgd2$obj, aic_bic_2, bcgd2$par
    ), nrow=length(lambda_vals_2)+1, ncol=length(drm$mele) + 5, byrow=TRUE
  )
  colnames(expected_return_matrix) = c("Lambda", "Iters", "Obj", "AIC", "BIC", names(drm$mele))
  
  # Run the test
  expect_equal(
    solutionPath(x_test, n_total_test, n_samples_test, m_test, model_test, d_test, lambdaVals = lambda_vals_2),
    expected_return_matrix
  )
  
  # Test that the lambda vals are sorting as they should be
  lambda_vals_3_a = c(1,2)
  lambda_vals_3_b = c(2,1)
  expect_equal(
    solutionPath(x_test, n_total_test, n_samples_test, m_test, model_test, d_test, lambdaVals = lambda_vals_3_a),
    solutionPath(x_test, n_total_test, n_samples_test, m_test, model_test, d_test, lambdaVals = lambda_vals_3_b)
  )
})