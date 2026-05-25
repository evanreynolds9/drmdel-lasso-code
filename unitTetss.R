# Test functions to ensure accuracy

# Load library
library(testthat)

# Load functions
source("R\\init.R")

### TEST DRMDELLASSO FUNCTIONS ###

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

# Test BCGD for data that was causing NAs prior to fix... ensure no error is throw
test_that("Test that the previously problematic output works okay, doesn't create NAs", {
  n_test = 50
  model_test = 12
  lambda_val_test = (n_test*3)^(1/3)*0.25
  d_test = 5
  m_test = 2
  n_samples_test = rep(n_test, m_test+1)
  model_test = 12
  
  # Here is the x-value - it happens on the second simulation with distribution=normal and paramSetup=2
  x_test = c(
    6.6309574, 3.1413579, 5.6474208, 6.6009971, 3.8154782, 2.9604379, 5.8827453, 
    4.3491227, 7.2508840, 1.6393554, 3.2628853, 4.7829593, 6.6365459, 4.3216228, 
    4.9365632, 6.7583933, 5.6616546, 4.8618242, 7.1794476, 5.2004390, 3.4153271, 
    4.9647271, 6.5926146, 3.3438906, 4.9524755, 4.9366501, 3.9358898, 5.0967965, 
    1.6581449, 6.9975067, 3.9228392, 4.2679283, 5.8949320, 3.9245702, 5.8373946, 
    7.3050936, 4.4784620, 5.7996959, 4.4287142, 3.9586621, 5.5625431, 2.0383891, 
    6.6154937, 6.4673656, 6.1316003, 2.4202314, 4.0600167, 6.5076318, 4.7709201, 
    7.3070676, 4.2019153, 4.6962754, 3.2233401, 5.4424402, 1.5348571, 4.9644339, 
    5.1357341, 3.7727272, 3.4366335, 3.2885508, 4.5677271, 5.4378566, 4.4858408, 
    5.0498646, 8.1368716, 3.7136574, 4.1579265, 5.9728063, 5.5325910, 3.9254638, 
    3.8152936, 2.4688471, 4.2907580, 3.8362365, 3.2533921, 3.4642051, 5.5902771, 
    3.8684889, 3.1718230, 3.5471559, 4.2786596, 3.4513952, 4.7384459, -0.1690255, 
    5.0261968, 1.6038025, 5.8403471, 5.2039663, 6.3906567, 4.5570211, 5.0176926, 
    4.9436399, 3.3677537, 5.8492259, 5.0942571, 2.8353539, 3.5413431, 5.9582463, 
    5.7291097, 4.4729758, 5.5601496, 6.5055580, 6.9798187, 4.6717211, 4.9190521, 
    5.5147621, 3.2505284, 4.8750979, 6.9785663, 4.5708385, 5.9111681, 5.8213275, 
    6.3228406, 5.9904783, 4.8181384, 6.2666558, 6.1426041, 4.5560148, 5.1078446, 
    5.2219760, 7.4200820, 6.8896429, 4.9345390, 4.0773469, 5.9655801, 5.5186582, 
    5.8661555, 6.8201747, 5.2725275, 6.1032118, 4.2611269, 3.1034228, 5.2048233, 
    5.5435162, 5.4554876, 4.9292365, 5.0635803, 5.1970922, 5.6578684, 6.6902415, 
    6.8124955, 4.2336403, 5.6004375, 6.5272864, 3.6415381, 3.7916745, 4.8084081, 
    3.4779163, 4.3930668, 5.8975375
  )
  
  # This is the parameter - returned by using the prior lambda value in our simulationRun
  theta = c(56.479963, 17.019179, 88.699394, 59.529174, -103.891429, 2.916249,-28.729046, 12.454162,-34.787001, -24.271241, 35.745495, -1.016633)
  
  # If there is an NA in the output of theta, an error will be thrown after running aic_bic_drm
  bcgd_problem = bcgd(theta_0 = theta, x = x_test, n_total = sum(n_samples_test), n_samples = n_samples_test, m = m_test, 
                      model = model_test, d = d_test, lambda = lambda_val_test)
  
  expect_no_error(aic_bic_drm(bcgd_problem$par, x_test, sum(n_samples_test), n_samples_test, m_test, model_test, d_test))
  
  # Try another case
  x_test_2 = c(
    5.08921501, 3.48895830, 4.10847868, 4.10070928, 6.31930075, 6.93327126, 5.98730380, 5.22725814, 2.23979247, 5.66423004, 5.96186074,
    3.11019455, 7.50072359, 2.77321179, 3.66800263, 5.85462413, 9.15618101, 4.12664830, 6.49964522, 5.39262095, 5.36686728, 4.68764314,
    5.64029205, 3.92765935, 5.56366088, 6.31668749, 6.22840515, 5.90480492, 6.02842074, 4.96615932, 4.02439362, 2.84294004, 4.75750274,
    5.49106221, 7.58223028, 5.48184413, 4.94249951, 2.37037457, 4.83748578, 4.04621983, 2.14868022, 3.84411765, 3.50198047, 3.22384735,
    5.18805203, 4.75995562, 5.41625396, 4.87584495, 2.83483231, 6.87051341, 5.25774785, 4.62752413, 5.14595067, 3.93844564, 4.67767032,
    3.87574114, 6.54660701, 4.97940653, 4.95170261, 1.95191979, 3.27220583, 5.19902487, 4.76591975, 4.12728271, 3.78907579, 4.95271583,
    0.06669263, 2.57965712, 3.78773122, 3.77236590, 4.32019831, 5.07886402, 5.01379912, 5.60866737, 8.87100559, 3.93686185, 4.63802233,
    5.12390967, 3.64986274, 4.23010373, 6.33887196, 5.80591074, 5.48947116, 4.83061755, 5.86497328, 3.03758350, 4.47151695, 5.63328336,
    5.22804092, 2.36825142, 5.70496762, 4.85276837, 3.57159435, 5.88168620, 5.34927824, 2.98098584, 3.74667964, 6.20513768, 2.80651358,
    3.80285838, 6.07885092, 7.05664401, 4.62162409, 5.15511362, 5.64876581, 5.64146158, 5.29983445, 6.50590245, 4.81449234, 7.50787595,
    6.17799486, 3.64160223, 7.09616238, 4.25209901, 6.91906517, 4.90423486, 5.22720090, 4.67816856, 7.43723780, 4.31475436, 5.70611733,
    4.88029070, 5.29327931, 3.62545503, 3.59481104, 3.79065323, 5.09221493, 5.66168246, 5.86207124, 3.75941521, 6.13068332, 6.39153640,
    5.83456851, 6.99305586, 6.64673158, 6.52450535, 5.38599599, 6.27122000, 3.88361344, 6.66454002, 6.33334507, 5.40377832, 6.31784487,
    3.89582589, 4.32587796, 5.54726012, 6.16724324, 4.72855890, 6.73962367, 5.07581358
  )
  theta_2 = c(56.479963, 17.019179, 88.699394, 59.529174, -103.891429, 2.916249,
              -28.729046, 12.454162,-34.787001, -24.271241, 35.745495, -1.016633)
  bcgd_problem = bcgd(theta_0 = theta_2, x = x_test_2, n_total = sum(n_samples_test), n_samples = n_samples_test, m = m_test, 
                      model = model_test, d = d_test, lambda = lambda_val_test)
  
  expect_no_error(aic_bic_drm(bcgd_problem$par, x_test_2, sum(n_samples_test), n_samples_test, m_test, model_test, d_test))
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
  
  # Confirm that function throws error if we pass a negative lambda value
  lambda_vals_0 = c(-1,1)
  expect_error(
    solutionPath(x_test, n_total_test, n_samples_test, m_test, model_test, d_test, lambdaVals = lambda_vals_0)
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


### TEST INIT FUNCTIONS ###

test_that("Test createRandomGenerator function that sets seed and returns a random number generator",{
  distribution_test = "gamma"
  paramSetup_test = 1
  n_test = 50
  
  # Test for some basic errors
  # Expect error when invalid distribution is passed
  expect_error(createRandomGenerator(distribution = "thisisnotadistn", paramSetup = paramSetup_test, 
                             n = n_test))
  # Expect error when passing an invalid paramsetup code
  expect_error(createRandomGenerator(distribution = distribution_test, paramSetup = 99999, 
                             n = n_test))
  
  # Test that seed setting works as expected
  test_generator = createRandomGenerator(distribution = distribution_test, paramSetup = paramSetup_test, 
                                         n = n_test)
  x1_test = test_generator()
  x2_test = test_generator()
  
  # Compute expected results
  set.seed(18)
  rate_0_test = 2
  rate_1_test = 1.4
  rate_2_test = 1.2
  shape_0_test = 1.8
  shape_1_test = 1.2
  shape_2_test = 1
  
  x11_expected = rgamma(n_test, shape=shape_0_test, rate=rate_0_test)
  x12_expected = rgamma(n_test, shape=shape_1_test, rate=rate_1_test)
  x13_expected = rgamma(n_test, shape=shape_2_test, rate=rate_2_test)
  x1_expected = c(x11_expected, x12_expected, x13_expected)
  
  x21_expected = rgamma(n_test, shape=shape_0_test, rate=rate_0_test)
  x22_expected = rgamma(n_test, shape=shape_1_test, rate=rate_1_test)
  x23_expected = rgamma(n_test, shape=shape_2_test, rate=rate_2_test)
  x2_expected = c(x21_expected, x22_expected, x23_expected)
  
  expect_equal(x1_test, x1_expected)
  expect_equal(x2_test, x2_expected)
})


test_that("Test the runSimulation function",{
  n_test = 50
  model_test = 6
  d_test = 2
  lambdaVals_test = c(0,1)
  n_samples_test = rep(n_test, 3)
  m_test = length(n_samples_test) - 1
  n_total_test = sum(n_samples_test)
  runs_test = 1
  
  # Load x data from mock file
  df = read.csv("mock_runSimulation.csv")
  x = as.matrix(df)
  colnames(x) = NULL
  
  x_test = x[1,]
  
  # Compute a second set of random variables - we will test two runs
  x_test_b = x[2,]
  
  # confirm function throws error if n and number of columns in x are incompatible
  expect_error(runSimulation(x = x, n=1, model = model_test, d = d_test, lambdaVals = lambdaVals_test))
  
  # Compute expected output
  simPath_test = solutionPath(x=x_test,n_total=n_total_test,n_samples=n_samples_test,
                              m=m_test,model=model_test,d=d_test,lambdaVals=lambdaVals_test)
  simPath_test = cbind(c(rep(1, nrow(simPath_test))), simPath_test)
  colnames(simPath_test) = c("Run", colnames(simPath_test[, -1]))
  
  # Compute output
  simPathOutput = runSimulation(x = matrix(x[1,], nrow=1), n=n_test, 
                                model=model_test, d=d_test,lambdaVals = lambdaVals_test)
  
  # Confirm equality
  expect_equal(simPath_test, simPathOutput)
  
  # Check that output is the same when 0 isn't in the lambda values
  simPathOutput_1b = runSimulation(x = matrix(x[1,], nrow=1), n=n_test, 
                                model=model_test, d=d_test,lambdaVals = c(1)) 
  
  expect_equal(simPathOutput, simPathOutput_1b)
  
  
  # Confirm for the adaptive case as well

  # Run computations
  simPath_test_adap = solutionPath(x=x_test,n_total=n_total_test,n_samples=n_samples_test,
                              m=m_test,model=model_test,d=d_test,lambdaVals=lambdaVals_test, adaptive = TRUE)
  simPath_test_adap = cbind(c(rep(1, nrow(simPath_test_adap))), simPath_test_adap)
  colnames(simPath_test_adap) = c("Run", colnames(simPath_test_adap[, -1]))
  
  # Compute output
  simPathOutputAdap = runSimulation(x = matrix(x[1,], nrow=1), n=n_test, 
                                model=model_test, d=d_test,lambdaVals = lambdaVals_test, adaptive = TRUE)
  # Expect equality again
  expect_equal(simPath_test_adap, simPathOutputAdap)
  
  # Run a test on multiple simulations to confirm matrix is being setup properly
  # Construct two run simulation from first run
  simPath_test_b = solutionPath(x=x_test_b,n_total=n_total_test,n_samples=n_samples_test,
                              m=m_test,model=model_test,d=d_test,lambdaVals=lambdaVals_test)
  simPath_test_b = cbind(c(rep(2, nrow(simPath_test_b))), simPath_test_b)
  simPath_test_2 = rbind(simPath_test, simPath_test_b) # This should inherit the colnames
  
  # Compute output
  simPathOutput_2 = runSimulation(x = x, n=n_test, 
                                model=model_test, d=d_test,lambdaVals = lambdaVals_test)
  
  # Confirm equality
  expect_equal(simPath_test_2, simPathOutput_2)
})


# Define helper function to expose actual basis functions for the next two tests
# This is needed because the helper will use the built-in C basis functions when available
expose_C_basis_function = function(model){
  if(is.double(model)){
    model_function = switch(
      as.character(model),
      `1` = function(x) x,
      `2` = function(x) log(abs(x)),
      `3` = function(x) sqrt(abs(x)),
      `4` = function(x) x^2,
      `5` = function(x) c(x, x^2),
      `6` = function(x) c(x, log(abs(x))),
      `7` = function(x) c(log(abs(x)), sqrt(abs(x)), x),
      `8` = function(x) c(log(abs(x)), sqrt(abs(x)), x^2),
      `9` = function(x) c(log(abs(x)), x, x^2),
      `10` = function(x) c(sqrt(abs(x)), x, x^2),
      `11` = function(x) c(log(abs(x)), sqrt(abs(x)), x, x^2),
      `12` = function(x) c(log(abs(x)), (log(abs(x)))^2, sqrt(abs(x)), x, x^2)
    )
  }else{
    model_function = model
  }
  return(model_function)
}


test_that("Test the basis function generator for the AIC/BIC simulations",{
  # The function produces every sub-function of this function: (log(x), (log(x))^2, sqrt(x), x, x^2)
  ids = 1:31
  
  x_test = 2
  # Create full test vector
  full_test_vector = c(log(x_test), log(x_test)^2, sqrt(x_test), x_test, x_test^2)
  subset_list_expected = list()
  
  # Create every possible subvector
  for(i in 1:5){
    subsets = combn(full_test_vector, i)
    for(j in 1:length(subsets[1,])){
      subset_list_expected = append(subset_list_expected, list(c(subsets[,j])))
    }
  }
  
  subset_list_test = list()
  for(i in ids){
    basis_func_or_id = subBasisFunc(i)
    if(is.double(basis_func_or_id)){ 
      if(basis_func_or_id == 6){
        # Need to hard code this as built-in C reverses this order
        func = function(x) c(log(abs(x)), x)
      }else{
        func = expose_C_basis_function(basis_func_or_id)
      }
    }else{
      func = expose_C_basis_function(basis_func_or_id)
    }
    subset_list_test = append(subset_list_test, list(func(x_test)))
  }
  
  expect_true(setequal(subset_list_expected, subset_list_test))
  
})


test_that("Test the simulation run function for AIC/BIC",{
  # The function produces a number of simulation runs using the AIC/BIC methodology
  # This computes the AIC/BIC of every sub-model of (log(x), log(x)^2, sqrt(x), x, x^2)
  # simulate data
  
  # Setup variables
  n_test = 50
  n_samples_test = rep(n_test, 3)
  n_total_test = sum(n_samples_test)
  ids = 1:31
  
  # These are fixed inputs of the function - will need to change if test changes
  max_basis_function_length = 5
  m_test = 2
  
  # Load x data from mock file
  df = read.csv("mock_runSimulation.csv")
  x = as.matrix(df)
  colnames(x) = NULL
  
  x_test_1 = x[1,]
  
  # Compute a second set of random variables - we will test two runs
  x_test_2 = x[2,]
  
  # confirm function throws error if n and number of columns in x are incompatible
  expect_error(simAICBIC(x = x, n=1))

  # Compute expected outputs
  total_cols = 5+(m_test*(max_basis_function_length+1))
  expected_output_1 = matrix(0, nrow=length(ids), ncol=total_cols)
  colnames(expected_output_1) = c(c("Run", "subFuncID", "d", "AIC", "BIC"), 
                                  paste("par", 1:(m_test*(max_basis_function_length+1)), sep="_"))
  
  # Compute outputs
  sim_output_1 = simAICBIC(x = matrix(x[1,], nrow=1), n=n_test)
  sim_output_2 = simAICBIC(x = x, n=n_test)
  
  for(id in ids){
    basis_func_for_id = expose_C_basis_function(subBasisFunc(id))
    d = length(basis_func_for_id(0))
    mele = drmdel(x = x_test_1, n_samples = n_samples_test, basis_func = basis_func_for_id)$mele
    aicbic = aic_bic_drm(theta=mele, x=x_test_1, n_total=n_total_test, n_samples=n_samples_test, 
                         m=m_test, basis_func=basis_func_for_id, d=d)
    
    expected_output_1[id, 1] = 1
    expected_output_1[id, 2] = id
    expected_output_1[id, 3] = d
    expected_output_1[id, 4:5] = aicbic
    expected_output_1[id, 6:(total_cols)] = c(mele, rep(0, total_cols-5-length(mele)))
  }
  
  # Compare to output
  expect_equal(expected_output_1, sim_output_1)
  
  # Compute second simulation
  expected_output_1_b = matrix(0, nrow=length(ids), ncol=total_cols)
  
  for(id in ids){
    basis_func_for_id = expose_C_basis_function(subBasisFunc(id))
    d = length(basis_func_for_id(0))
    mele = drmdel(x = x_test_2, n_samples = n_samples_test, basis_func = basis_func_for_id)$mele
    aicbic = aic_bic_drm(theta=mele, x=x_test_2, n_total=n_total_test, n_samples=n_samples_test, 
                         m=m_test, basis_func=basis_func_for_id, d=d)
    
    expected_output_1_b[id, 1] = 2
    expected_output_1_b[id, 2] = id
    expected_output_1_b[id, 3] = d
    expected_output_1_b[id, 4:5] = aicbic
    expected_output_1_b[id, 6:(total_cols)] = c(mele, rep(0, total_cols-5-length(mele)))
  }
  
  expected_output_2 = rbind(expected_output_1, expected_output_1_b)
  
  expect_equal(expected_output_2, sim_output_2)
})


test_that("Test the summariseSim function on a test csv files",{
  results1a = summariseSim(distribution="normal", file_name="mock_summariseSim.csv", basis_func=12, tol=0.0001)
  expected_results1a = c(runs=1, AIC=1, AIC_sub=1, BIC=0.5, BIC_sub=1)
  expect_equal(results1a, expected_results1a)
  
  results1b = summariseSim(distribution="gamma", file_name="mock_summariseSim.csv", basis_func=12, tol=0.0001)
  expected_results1b = c(runs=0, AIC=0, AIC_sub=0, BIC=0, BIC_sub=0.5)
  expect_equal(results1b, expected_results1b)
  
  results1c = summariseSim(distribution="normal", file_name="mock_summariseSim.csv", basis_func=12, tol=0.005)
  expected_results1c = c(runs=0.5, AIC=0.5, AIC_sub=0.5, BIC=0, BIC_sub=0.5)
  expect_equal(results1c, expected_results1c)
})

