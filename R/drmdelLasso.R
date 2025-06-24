# ##############################
# This software is written by Song Cai and published under GPLv3.
#
# Version 1.3.1, December 31, 2014.
#
# Lasso functionality added by Evan Reynolds, 2024
# ##############################

bcgd = function(theta_0, x ,n_total, n_samples, m, d, model, lambda, pen_g = rep(1,d),
                max_iters = 1000, threshold = 1e-6, omega_0 = 1, psi = 0.5, sigma = 0.1){
  # Optimize the drm group lasso objective function with bcgd,
  # starting at a given value of parameter.
  #
  # double *restrict n_total, /* total pooled sample size */
  # double *restrict n_samples, /* vector samples sizes*/
  # double *restrict m, double *restrict d, /* number of samples and length of basis function*/
  # double *restrict model, double *restrict x, /*vector of sample values*/
  # double *restrict lambda, double *restrict pen_g, /* penalty values*/
  # double *restrict omega_0, double *restrict psi, double *restrict sigma, /* optimization hyperparameters*/
  # double *threshold, double *max_iters, /*convergence criteria*/
  # double *theta_0, /*input: initial value of parameter vector*/
  # double *opt_val, /*output: the minimized value of the negLDLGL function */
  # double *total_iters /*output: number of iterations to convergence*/
  
  # End function if pen_g has length not equal to d
  if(length(pen_g) != d){
    print("Length of penalty does not equal length of basis function!")
    return(-1)
  }
  
  
  bcgdOutput = .C("bcgd", as.double(n_total), as.double(n_samples), as.double(m),
                  as.double(d), as.double(model), as.double(x), as.double(lambda), as.double(pen_g),
                  as.double(omega_0), as.double(psi), as.double(sigma),
                  as.double(threshold), as.double(max_iters),
                  theta_f = as.double(theta_0), opt_val = as.double(0), total_iters = as.double(0))
  
  return(list(obj = bcgdOutput$opt_val, iters = bcgdOutput$total_iters, par = bcgdOutput$theta_f))
}

# Create aic/bic function
aic_bic_drm = function(theta, x, n_total ,n_samples, m, basis_func, d){
  
  # Compute the AIC and BIC of a given DRM parameter estimate
  negLDLVal = negLDL(theta, x, n_total, n_samples, m, basis_func, d)
  theta_mat = matrix(theta, nrow = m, ncol = d+1, byrow = TRUE)[,-1] # Remove first column for constants
  group_sums = colSums(theta_mat^2)
  
  # Count the number of non-zero groups for the penalty
  pen = sum(group_sums>1e-12)*m # Set tolerance in case values are extremely close to 0
  return(c(2*negLDLVal+2*pen,2*negLDLVal+log(n_total)*pen))
}

solutionPath = function(x ,n_total, n_samples, m, d, model, lambdaVals, adaptive = FALSE){
  # Generate a solution path for given data
  # Can compute multiple solution paths to accommodate simulation studies
  # By default, the solution path will always
  # 1. Use the default hyper-parameters of the bcgd function
  # 2. Start with the initial parameter as the MELE, which will also be the first value in the solution path
  # The arguments are the usual drmdel arguments, except for lambdaVals which should be the grid of lambda values
  # to be evaluated over in the solution path, and a flag for whether to use the adaptive lasso (default is FALSE)
  
  # Abort process if any lambda's are not positive
  if(sum(lambdaVals < 0) != 0){
    print("Negative lambda values detected in solution path, aborting simulation.")
    return(-1)
  }
  
  # Check if 0 is in lamdba_vals, add it if not
  if(!(0 %in% lambdaVals)){
    lambdaVals = c(0,lambdaVals)
  }
  
  # Sort lambdaVals
  lambdaVals = sort(lambdaVals)
    
  # Create matrix to store simulation results, initialize with 0s
  totalCols = m*(d+1) + 5
  sim_results = matrix(0, nrow = length(lambdaVals), ncol = totalCols)
  # compute the mele
  mele_sol = drmdel(x=x, n_samples=n_samples, basis_func=model)
  mele = mele_sol$mele
  mele_obj = mele_sol$negldl
  # Set names of matrix columns
  colnames(sim_results) = c("Lambda", "Iters", "Obj", "AIC", "BIC",
                                names(mele))
    
  # If adaptive is true, we set the inverse of the groupwise mele to be pen_G
  if(adaptive){
    theta_mat = matrix(mele, nrow = m, ncol = d+1, byrow = TRUE)
    pen_g = 1/(sqrt(colSums(theta_mat^2))[-1]) # Delete column for constants
  } else{
    pen_g = rep(1,d)
  }
  # Set value for initial guess of theta to be the mele
  init_theta = mele
  # Iterate over lambda vals to populate matrix for each run
  # The matrix columns are set as follows
  # 1. The lambda value
  # 2. The number of iteration to convergence
  # 3. The objective function value
  # 4. the AIC value
  # 5. The BIC value
  # 6 to 6+m*(d+1): The parameter values
  for(i in 1:length(lambdaVals)){
    lambda = lambdaVals[i]
    sim_results[i, 1] = lambda
    if(lambda == 0){ # Compute mele
      sim_results[i, 2] = 0
      sim_results[i, 3] = mele_obj
      aic_bic = aic_bic_drm(mele, x, n_total ,n_samples, m, model, d)
      sim_results[i, 4] = aic_bic[1]
      sim_results[i, 5] = aic_bic[2]
      sim_results[i, 6:totalCols] = mele
    } else {
      # Compute bcgd optimization
      bcgdOpts = bcgd(init_theta, x ,n_total, n_samples, m, d, model, lambda, pen_g = pen_g)
      sim_results[i, 2] = bcgdOpts$iters
      sim_results[i, 3] = bcgdOpts$obj
      aic_bic = aic_bic_drm(bcgdOpts$par, x, n_total ,n_samples, m, model, d)
      sim_results[i, 4] = aic_bic[1]
      sim_results[i, 5] = aic_bic[2]
      sim_results[i, 6:totalCols] = bcgdOpts$par
      if(i < length(lambdaVals)){
        init_theta = bcgdOpts$par
      }
    }
  }
  # Return matrix
  return(sim_results)
}

negLDL <- function(par, x, n_total, n_samples, m, model, d) {
# Calculate negative log dual empirical likelihood for a
# given value of parameter.
#
# logDualLWrapper prototype:
# void logDualLWrapper( double * restrict n_total, /*inputs*/
#     double * restrict n_samples, /*inputs*/
#     double * restrict m, double * restrict d,
#     double * restrict par, /*inputs*/
#     double * restrict model, double * restrict x, /*inputs*/
#     double * restrict ldl_val /*output*/)

  logDL <- .C("logDualLWrapper", as.double(n_total),
              as.double(n_samples), as.double(m), as.double(d),
              as.double(par), as.double(model), as.double(x),
              ldl_val=double(1))

  return(-logDL$ldl_val)

}

negLDLGL <- function(par, x, n_total, n_samples, m, model, d, lambda, pen_g = rep(1,d)) {
  # Calculate negative log dual empirical likelihood for a
  # given value of parameter, with a GL penalty applied across
  # the coefficients of each basis function term and a specified
  # penalty threshold.
  #
  # logDualLGLWrapper prototype:
  # void logDualLGLWrapper( double * restrict n_total, /*inputs*/
  #     double * restrict n_samples, /*inputs*/
  #     double * restrict m, double * restrict d,
  #     double * restrict par, /*inputs*/
  #     double * restrict model, double * restrict x, /*inputs*/
  #     double * restrict lambda, double * restrict pen_g, /*inputs*/
  #     double * restrict ldl_val /*output*/)
  
  logDL <- .C("logDualLGLWrapper", as.double(n_total),
              as.double(n_samples), as.double(m), as.double(d),
              as.double(par), as.double(model), as.double(x), 
              as.double(lambda), as.double(pen_g),
              ldlGL_val=double(1))
  
  return(-logDL$ldlGL_val)
  
}

# User specified basis function version
negLDLUf <- function(par, x, n_total, n_samples, m, basis_func, d) {
# Calculate negative log dual empirical likelihood for a
# given value of parameter.
#
# basis_func must be a function.
#
# logDualLUfWrapper prototype:
# SEXP logDualLUfWrapper(SEXP n_total, SEXP n_samples, SEXP m, SEXP d, SEXP par,
#     SEXP h_func, SEXP env, SEXP x)


  val <- .Call("logDualLUfWrapper", as.double(n_total), as.double(n_samples),
        as.double(m), as.double(d), as.double(par),
        basis_func, new.env(), as.double(x))

  return(-val)

}


negLDLGr <- function(par, x, n_total, n_samples, m, model, d) {
# Calculate the gradient of the negative log dual empirical
# likelihood for a given value of parameter.
#
# logDualLDLGrWrapper prototype:
# void logDualLGrWrapper( double * restrict n_total, /*inputs*/
#     double * restrict n_samples, /*inputs*/
#     double * restrict m, double * restrict d,
#     double * restrict par, /*inputs*/
#     double * restrict model, double * restrict x, /*inputs*/
#     double * restrict ldl_gradient /*output*/)


  logDLGr <- .C("logDualLGrWrapper", as.double(n_total),
              as.double(n_samples), as.double(m), as.double(d),
              as.double(par), as.double(model), as.double(x),
              ldl_gr=double(m*(d+1)))

  return(-logDLGr$ldl_gr)

}

# User specified basis function version
negLDLGrUf <- function(par, x, n_total, n_samples, m, basis_func, d) {
# Calculate the gradient of the negative log dual empirical
# likelihood for a given value of parameter.
#
# basis_func must be a function.
#
# logDualLGrUfWrapper prototype:
# SEXP logDualLGrUfWrapper(SEXP n_total, SEXP n_samples, SEXP m, SEXP d, SEXP par,
#     SEXP h_func, SEXP env, SEXP x /*input*/)

  ldl_gr <- .Call("logDualLGrUfWrapper", as.double(n_total),
                  as.double(n_samples), as.double(m), as.double(d),
                  as.double(par), basis_func, new.env(), as.double(x))

  return(-ldl_gr)

}

gen_par_pos <- function(m, d){

  par_pos_alpha <- seq(from=1, by=(d+1), length.out=m) 
  par_pos_beta <- 1:(m*d) + rep(1:m, each=d)

  return(list(alpha=par_pos_alpha, beta=par_pos_beta))

}

negLDL_null <- function(par_null_full, g_null, g_null_jac=NULL,
                        par_pos, par_dim, par_dim_null=NULL,
                        x, n_total, n_samples, m, model, d) {
# g_null is a null mapping from gamma to beta, but
# par_null_full must be the full null parameter that
# includes alpha, i.e. par_null_full=(alpha, gamma).
#   alpha: the first m elements of par_null_full
#   gamma: the rest (from (m+1)th to the last) elements of
#     par_null_full
#
# par_dim is the dimension of DRM parameter $theta$, which
#   equals m*(d+1).
# par_dim_null is the dimension of the null mapping, i.e.
#   dim(gamma).

  par <- numeric(par_dim)
  par[par_pos$alpha] <- par_null_full[1:m]  # get alpha
  par[par_pos$beta] <- g_null(par_null_full[-(1:m)])  # get beta

  return(negLDL(par=par, x=x, n_total=n_total, n_samples=n_samples, m=m,
         model=model, d=d))

}

# User specified basis function version
negLDLUf_null <- function(par_null_full, g_null, g_null_jac=NULL,
                          par_pos, par_dim, par_dim_null=NULL,
                          x, n_total, n_samples, m, basis_func, d) {
# basis_func must be a function.
#
# g_null is a null mapping from gamma to beta, but
# par_null_full must be the full null parameter that
# includes alpha, i.e. par_null_full=(alpha, gamma).
#   alpha: the first m elements of par_null_full
#   gamma: the rest (from (m+1)th to the last) elements of
#     par_null_full
#
# par_dim is the dimension of DRM parameter $theta$, which
#   equals m*(d+1).
# par_dim_null is the dimension of the null mapping, i.e.
#   dim(gamma).

  par <- numeric(par_dim)
  par[par_pos$alpha] <- par_null_full[1:m]  # get alpha
  par[par_pos$beta] <- g_null(par_null_full[-(1:m)])  # get beta

  return(negLDLUf(par=par, x=x, n_total=n_total, n_samples=n_samples, m=m,
         basis_func=basis_func, d=d))

}

negLDLGr_null <- function(par_null_full, g_null, g_null_jac,
                          par_pos, par_dim, par_dim_null,
                          x, n_total, n_samples, m, model, d) {
# par_null_full must be the full null parameter that
# includes alpha, i.e. par_null_full=(alpha, gamma).
#
# g_null is a null mapping from gamma to beta.
# g_null_jac is the jacobian matrix of g_null, i.e. a matrix
#   of dimension m*d by dim(gamma).
#
# par_dim is the dimension of DRM parameter $theta$, which
#   equals m*(d+1).
# par_dim_null is the dimension of the null mapping, i.e.
#   dim(gamma).

  par <- numeric(par_dim)
  par[par_pos$alpha] <- par_null_full[1:m]  # get alpha
  par_gamma <- par_null_full[-(1:m)]  # get gamma
  par[par_pos$beta] <- g_null(par_gamma)

  par_full_gr <- negLDLGr(par=par, x=x, n_total=n_total, n_samples=n_samples,
                          m=m, model=model, d=d)

  # Jacobian matrix for g_null_full (null mapping from (alpha,
  # gammma) to (theta_1, theta_2, ..., theta_m))
  par_dim_null_full <- m + par_dim_null
  g_null_full_jac <- matrix(rep(0, par_dim*par_dim_null_full),
                                nrow=par_dim,
                                ncol=par_dim_null_full)

  if (par_dim_null > 0) {
    g_null_full_jac[par_pos$beta, (m+1):par_dim_null_full] <- g_null_jac(par_gamma)
  }
  g_null_full_jac[par_pos$alpha, 1:m] <- diag(m) 

  return(as.numeric(t(g_null_full_jac) %*% par_full_gr))

}

# User specified basis function version
negLDLGrUf_null <- function(par_null_full, g_null, g_null_jac, 
                            par_pos, par_dim, par_dim_null,
                            x, n_total, n_samples, m, basis_func, d) {
# basis_func must be a function.
#
# par_null_full must be the full null parameter that
# includes alpha, i.e. par_null_full=(alpha, gamma).
#
# g_null is a null mapping from gamma to beta.
# g_null_jac is the jacobian matrix of g_null, i.e. a matrix
#   of dimension m*d by dim(gamma).
#
# par_dim is the dimension of DRM parameter $theta$, which
#   equals m*(d+1).
# par_dim_null is the dimension of the null mapping, i.e.
#   dim(gamma).

  par <- numeric(par_dim)
  par[par_pos$alpha] <- par_null_full[1:m]  # get alpha
  par_gamma <- par_null_full[-(1:m)]  # get gamma
  par[par_pos$beta] <- g_null(par_gamma)

  par_full_gr <- negLDLGrUf(par=par, x=x, n_total=n_total, n_samples=n_samples,
                            m=m, basis_func=basis_func, d=d)

  # Jacobian matrix for g_null_full (null mapping from (alpha,
  # gammma) to (theta_1, theta_2, ..., theta_m))
  par_dim_null_full <- m + par_dim_null
  g_null_full_jac <- matrix(rep(0, par_dim*par_dim_null_full),
                                nrow=par_dim,
                                ncol=par_dim_null_full)

  if (par_dim_null > 0) {
    g_null_full_jac[par_pos$beta, (m+1):par_dim_null_full] <- g_null_jac(par_gamma)
  }
  g_null_full_jac[par_pos$alpha, 1:m] <- diag(m) 

  return(as.numeric(t(g_null_full_jac) %*% par_full_gr))

}

negLDLHessian <- function(par, x, n_total, n_samples, m, model, d) {
# Calculate the hessian of the negative log dual empirical
# likelihood for a given value of parameter.
#
# logDualLDLHessianWrapper prototype:
# void logDualLHessianWrapper( double * restrict n_total, /*inputs*/
#     double * restrict n_samples, /*inputs*/
#     double * restrict m, double * restrict d,
#     double * restrict par, /*inputs*/
#     double * restrict model, double * restrict x, /*inputs*/
#     double * restrict ldl_hessian /*output*/)


  logDLHessian <- .C("logDualLHessianWrapper", as.double(n_total),
              as.double(n_samples), as.double(m), as.double(d),
              as.double(par), as.double(model), as.double(x),
              ldl_hessian=double(m*(d+1)*m*(d+1)))

  return( matrix(-(logDLHessian$ldl_hessian), m*(d+1), m*(d+1),
                 byrow=TRUE) )

}

# User specified basis function version
negLDLHessianUf <- function(par, x, n_total, n_samples, m, basis_func, d) {
# Calculate the hessian of the negative log dual empirical
# likelihood for a given value of parameter.
#
# basis_func must be a function.
#
# logDualLDLHessianUfWrapper prototype:
# SEXP logDualLHessianUfWrapper(SEXP n_total, SEXP n_samples, SEXP m, SEXP d,
#     SEXP par, SEXP h_func, SEXP env, SEXP x /*input*/)

  ldl_hessian <- .Call("logDualLHessianUfWrapper", as.double(n_total),
                       as.double(n_samples), as.double(m), as.double(d),
                       as.double(par), basis_func, new.env(), as.double(x))

  return(-ldl_hessian)

}

# User specified basis function version


pEst <- function(x, n_total, n_samples, m, model, d, mele) {
# Extract estimated probabilities given MELE
#
# probEstWrapper prototype:
# void probEstWrapper( double * restrict n_total, /*inputs*/
#     double * restrict n_samples, /*inputs*/
#     double * restrict m, double * restrict d,
#     double * restrict par, /*inputs*/
#     double * restrict model, double * restrict x, /*inputs*/
#     double * restrict normalize,
#     double * restrict pEst /*output*/)

  #n_total <- sum(n_samples)
  #m <- length(n_samples) - 1

  pEstimates <- .C("probEstWrapper", as.double(n_total),
                      as.double(n_samples), as.double(m), as.double(d),
                      as.double(mele), as.double(model), as.double(x),
                      1.0,
                      p_est=double((m+1)*n_total))

  population_indicator <- rep(x=seq(0, m, by=1), times=rep(n_total, m+1))

  return(data.frame(k=population_indicator, x=x, p_est=pEstimates$p_est))

}

# User specified basis function version
pEstUf <- function(x, n_total, n_samples, m, basis_func, d, mele) {
# Extract estimated probabilities given MELE
#
# basis_func must be a function.
#
# probEstUfWrapper prototype:
# SEXP probEstUfWrapper(SEXP n_total, SEXP n_samples, SEXP m, SEXP d,
#     SEXP par, SEXP h_func, SEXP env, SEXP x, SEXP normalize /*input*/)

  #n_total <- sum(n_samples)
  #m <- length(n_samples) - 1

  p_est <- .Call("probEstUfWrapper", as.double(n_total),
                    as.double(n_samples), as.double(m), as.double(d),
                    as.double(mele), basis_func, new.env(), as.double(x),
                    1.0)

  population_indicator <- rep(x=seq(0, m, by=1), times=rep(n_total, m+1))

  return(data.frame(k=population_indicator, x=x, p_est=p_est))

}


cdfEst <- function(x, n_samples, p_est) {
# Estimate distribution functions F_k's, k = 0, 1, ..., m.

  # order of x
  x_order <- order(x)
  x_sort <- x[x_order]

  # extract p_est for k^{th} sample
  m <- length(n_samples) - 1
  n_total <- sum(n_samples)
  cum_p <- NULL
  for (i in 1:(m+1)) {
    pe_tmp <- p_est[p_est$k==(i-1), 3]
    cum_p <- c(cum_p, cumsum(pe_tmp[x_order]))
  }

  population_indicator <- rep(x=seq(0, m, by=1), times=rep(n_total, m+1))
  x_data <- rep(x=x_sort, m+1)

  return(data.frame(k=population_indicator, x=x_data, cdf_est=cum_p))

}

quantEst <- function(k, p, cdf_est, n_samples,
                     interpolation=TRUE, adj=FALSE,
                     adj_val=NULL) {
# Estimate the quantile of the k[i]^th, k[i] = 0, 1, ..., m, population at
#   probability p[i] without giving covariance estimates.

  # Arguments handling and checking
  if ((length(k) != length(p)) && (length(k) != 1 ) && (length(p) != 1)) {
    stop("The lengths of vector 'k' and 'p' must be the same, or one of the vector must have length one!")
  }
  nK <- max(length(k), length(p))
  if (length(k)==1) k <- rep(k, nK)
  if (length(p)==1) p <- rep(p, nK)

  if (!is.logical(interpolation)) {
    stop("The argument 'interpolation' must be a logical variable (either TRUE or FALSE)!")
  }

  if (!is.logical(adj)) {
    stop("The argument 'adj' must be a logical variable (either TRUE or FALSE)!")
  }

  if (adj==TRUE) {

    if (is.null(adj_val)) {

      adj_val <- -1/(2*n_samples[(k+1)])

    } else {

      if (!is.numeric(adj_val)) {
        stop("The argument 'adj_val' must either be NULL or a numerical value vector!")
      } else if (length(adj_val) != 1 && length(adj_val) != nK) {
        stop("The length of the numerical argument 'adj_val' must either be 1 or the same as the length of the argument 'k' or, when the length of 'k' is 1, the same as the length of the argument 'p')!")
      }

      if (length(adj_val) == 1) adj_val <- rep(adj_val, nK)

    }

  }

  # Extract sorted data and useful information about DRM 
  x_sort <- cdf_est[cdf_est$k==0,]$x

  qe <- numeric(nK)
  for (i in 1:nK) {
    cdf_est_tmp <- cdf_est[cdf_est$k==k[i], 3]
    if (adj==TRUE) {
      cdf_est_tmp <- cdf_est_tmp + adj_val[i]
    }
    if (p[i] >= tail(cdf_est_tmp, 1)) {
      pos_tmp <- length(x_sort) - 1
    } else {
      pos_tmp <- length(which(cdf_est_tmp < p[i]))
    }

    if (interpolation==TRUE) {

      # linear interpolation
      qe[i] <- x_sort[pos_tmp] + (x_sort[pos_tmp + 1] - x_sort[pos_tmp]) *
        (p[i] - cdf_est_tmp[pos_tmp]) /
          (cdf_est_tmp[pos_tmp+1] - cdf_est_tmp[pos_tmp])

    } else {

      # without linear interpolation
      qe[i] <- x_sort[pos_tmp + 1]

    }

  }

  return(qe)

}

summaryDRMFEst <- function(n_samples, p_est, cdf_est, interpolation=TRUE) {
# Estimated mean, variance, IQR, of distribution functions F_k's, k = 0, 1, ..., m.

  m <- length(n_samples) - 1
  n_total <- sum(n_samples)

  # Estimated means
  mu_est <- numeric(m+1)
  var_est <- numeric(m+1)
  x <- p_est[p_est$k==0,2]
  for (i in 1:(m+1)) {
    pe_tmp <- p_est[p_est$k==(i-1),3]
    mu_est[i] <- sum(x * pe_tmp)
    var_est[i] <- sum(x^2 * pe_tmp) - mu_est[i]^2
  }
  sd_est <- sqrt(var_est)

  k <- 0:m
  q1 <- quantEst(k, 0.25, cdf_est, n_samples,
                 interpolation=interpolation, adj=FALSE,
                 adj_val=NULL)
  q3 <- quantEst(k, 0.75, cdf_est, n_samples,
                 interpolation=interpolation, adj=FALSE,
                 adj_val=NULL)
  iqr <- q3 - q1


  result <- data.frame(mean=mu_est, var=var_est, sd=sd_est,
                       Q1=q1, Q3=q3, IQR=iqr)

  rnames <- c("F0")
  for (i in 1:m) {
    rnames <- c(rnames, paste("F", i, sep=""))
  }
  row.names(result) = rnames 

  return(result)
}

displayPar <- function(par, m) {

  d <- length(par)/m - 1

  par_out <- matrix(par, d+1, m) 
  par_out <- t(par_out)
  par_out <- data.frame(par_out)
  rnames <- NULL
  for (i in 1:m) {
    #rnames <- c(rnames, paste("theta[", i, "]", sep=""))
    rnames <- c(rnames, paste("F", i, sep=""))
  }
  cnames <- paste("alpha[]")
  for (i in 1:d) {
    cnames <- c(cnames, paste("beta[,", i, "]", sep=""))
  }
  row.names(par_out) <- rnames
  colnames(par_out) <- cnames

  return(par_out)

}

drmdel <- function(x, n_samples, basis_func, g_null=NULL,
                   g_null_jac=NULL, par_dim_null=NULL,
                   ...)
{
  # setting useful constants
  n_total <- sum(n_samples)
  if (length(x) != n_total) {
    stop("Length of data vector 'x' is not consistent with sample sizes 'n_samples'!")
  }  # check data lengths
  rho <- n_samples/n_total
  m <- length(n_samples) - 1

  if (is.function(basis_func)) {
    d <- length(basis_func(x[1]))
  } else if (basis_func %in% 1:4) {
    d <- 1
  } else if (basis_func %in% 5:6) {
    d <- 2
  } else if (basis_func %in% 7:10) {
    d <- 3
  } else if (basis_func == 11){
    d <- 4
  } else if (basis_func == 12){
    d <- 5
  }else {
    stop("Parameter 'basis_func' must be a function of a single variable or an integer between 1 and 11!")
  }

  drm_info <- list(m=m, d=d, basis_func=basis_func,
                   n_samples=n_samples, n_total=n_total, rho=rho)

  par_dim <- m*(d+1)

  par_pos <- gen_par_pos(m, d)

  # extracting ... arguments for function optim()
  dot_args <- list(...)
  # setting default method as "BFGS" for function optim()
  if (is.null(dot_args$method)) dot_args <- c(dot_args, list(method="BFGS"))
  # setting default maximum iteration number to 1000 for function optim()
  if (is.null(dot_args$control)) {
    #dot_args <- c(dot_args, list(control=list(maxit=10000,
                                              #reltol=.Machine$double.eps)))
    dot_args <- c(dot_args, list(control=list(maxit=10000)))
  } else {

    if (is.null(dot_args$control$maxit)) {
      dot_args$control <- c(dot_args$control, list(maxit=10000))
    }

    #if (is.null(dot_args$control$reltol)) {
      #dot_args$control <- c(dot_args$control,
                            #list(reltol=.Machine$double.eps))
    #}

  }
  # safeguard one-dimensional case
  if (par_dim==1) {
    dot_args$method="Brent"
    if (is.null(dot_args$lower)) dot_args <- c(dot_args, list(lower=-100))
    if (is.null(dot_args$upper)) dot_args <- c(dot_args, list(upper=100))
  }

  # calculating MELE
  # About initial values: it is better to always set initial values to zeros
  # for all parameters to ensure that, at the initial value, the dual
  # log empirical likelihood is finite!
  par_init <- rep(0, par_dim)  # 

  if (is.function(basis_func)) {

    drm_opt <- do.call(optim, c(list(par=par_init, fn=negLDLUf,
                                     gr=negLDLGrUf),
                                dot_args,
                                list(x=x, n_total=n_total,
                                     n_samples=n_samples, m=m,
                                     basis_func=basis_func,
                                     d=d)))
    mele <- drm_opt$par
    mnames <- NULL
    for (i in 1:m) {
      mnames <- c(mnames, paste("alpha[", i, "]", sep=""))
      for (j in 1:d) {
        mnames <- c(mnames, paste("beta[", i, ",", j, "]", sep=""))
      }
    }
    names(mele) <- mnames

    negldl <- drm_opt$value

    # ##### Estimate information matrix #####
    logDLHessianMdele <- .Call("logDualLHessianUfWrapper", as.double(n_total),
                               as.double(n_samples), as.double(m),
                               as.double(d), as.double(mele), basis_func,
                               new.env(), as.double(x))

    # Observed information matrix
    info_mat <- matrix(-logDLHessianMdele/n_total, par_dim, par_dim,
                       byrow=TRUE)

    rm(logDLHessianMdele)
    # ##########

    # estimate dF_{k}(x_kj)
    p_est <- pEstUf(x=x, n_total=n_total, n_samples=n_samples, m=m,
                    basis_func=basis_func, d=d, mele=mele)

    # estimate F_{k}(x_kj)
    cdf_est <- cdfEst(x=x, n_samples=n_samples, p_est=p_est)

    if (!is.null(g_null)) {

      if (is.null(par_dim_null)) {
        stop("One must provide 'par_dim_null', the dimension of the null parameter!")
      }

      par_init_null_full <- rep(0, (m+par_dim_null))
      par_dim_null_full <- par_dim_null + m

      dot_args_null <- dot_args
      # safeguard one-dimensional case
      if (par_dim_null_full==1) {
        dot_args_null$method="Brent"
        if (is.null(dot_args_null$lower)) {
          dot_args_null <- c(dot_args_null, list(lower=-100))
        }
        if (is.null(dot_args_null$upper)) {
          dot_args_null <- c(dot_args_null, list(upper=100))
        }
      }

      if (!is.null(g_null_jac) || par_dim_null==0) {
        drm_opt_null <- do.call(optim, c(list(par=par_init_null_full,
                                              fn=negLDLUf_null,
                                              gr=negLDLGrUf_null),
                                         dot_args_null,
                                         list(g_null=g_null,
                                              g_null_jac=g_null_jac,
                                              par_pos=par_pos,
                                              par_dim=par_dim,
                                              par_dim_null=par_dim_null,
                                              x=x, n_total=n_total,
                                              n_samples=n_samples, m=m,
                                              basis_func=basis_func, d=d)))
      } else {
        drm_opt_null <- do.call(optim, c(list(par=par_init_null_full,
                                              fn=negLDLUf_null),
                                         dot_args_null,
                                         list(g_null=g_null,
                                              par_pos=par_pos,
                                              par_dim=par_dim,
                                              x=x, n_total=n_total,
                                              n_samples=n_samples, m=m,
                                              basis_func=basis_func, d=d)))
      }

      negldl_null <- drm_opt_null$value
      mele_null <- list(alpha=drm_opt_null$par[1:m],
                        gamma=drm_opt_null$par[-(1:m)])
      delr <- -2*(negldl - negldl_null)  # DEL ratio statistics
      df <- par_dim - par_dim_null_full
      p_val <- 1-pchisq(delr, df)

      return(list(drm_info=drm_info, mele=mele,
                  info_mat=info_mat, negldl=negldl,
                  mele_null=mele_null, negldl_null=negldl_null,
                  delr=delr, df=df, p_val=p_val,
                  p_est=p_est, cdf_est=cdf_est))

    } else {

      delr <- -2*negldl  # DEL ratio statistics
      df <- m*d
      p_val <- 1-pchisq(delr, df)

      return(list(drm_info=drm_info, mele=mele,
                  info_mat=info_mat, negldl=negldl,
                  delr=delr, df=df, p_val=p_val,
                  p_est=p_est, cdf_est=cdf_est))

    }

  } else {

    drm_opt <- do.call(optim, c(list(par=par_init, fn=negLDL, gr=negLDLGr),
                                dot_args,
                                list(x=x, n_total=n_total,
                                     n_samples=n_samples, m=m,
                                     model=basis_func,
                                     d=d)))
    mele <- drm_opt$par
    mnames <- NULL
    for (i in 1:m) {
      mnames <- c(mnames, paste("alpha[", i, "]", sep=""))
      for (j in 1:d) {
        mnames <- c(mnames, paste("beta[", i, ",", j, "]", sep=""))
      }
    }
    names(mele) <- mnames

    negldl <- drm_opt$value

    # ##### Estimate information matrix #####
    logDLHessianMdele <- .C("logDualLHessianWrapper", as.double(n_total),
                as.double(n_samples), as.double(m), as.double(d),
                as.double(mele), as.double(basis_func), as.double(x),
                ldl_hessian=double(par_dim*par_dim))

    # Observed information matrix
    info_mat <- matrix(-(logDLHessianMdele$ldl_hessian)/n_total, par_dim,
                       par_dim, byrow=TRUE)

    rm(logDLHessianMdele)
    # ##########

    # estimate dF_{k}(x_kj)
    p_est <- pEst(x=x, n_total=n_total, n_samples=n_samples, m=m,
                  model=basis_func, d=d, mele=mele)

    # estimate F_{k}(x_kj)
    cdf_est <- cdfEst(x=x, n_samples=n_samples, p_est=p_est)

    if (!is.null(g_null)) {

      if (is.null(par_dim_null)) {
        stop("One must provide 'par_dim_null', the dimension of the null parameter!")
      }

      par_init_null_full <- rep(0, (m+par_dim_null))
      par_dim_null_full <- par_dim_null + m

      dot_args_null <- dot_args
      # safeguard one-dimensional case
      if (par_dim_null_full==1) {
        dot_args_null$method="Brent"
        if (is.null(dot_args_null$lower)) {
          dot_args_null <- c(dot_args_null, list(lower=-100))
        }
        if (is.null(dot_args_null$upper)) {
          dot_args_null <- c(dot_args_null, list(upper=100))
        }
      }

      if (!is.null(g_null_jac) || par_dim_null==0) {
        drm_opt_null <- do.call(optim, c(list(par=par_init_null_full,
                                              fn=negLDL_null,
                                              gr=negLDLGr_null),
                                         dot_args_null,
                                         list(g_null=g_null,
                                              g_null_jac=g_null_jac,
                                              par_pos=par_pos,
                                              par_dim=par_dim,
                                              par_dim_null=par_dim_null,
                                              x=x, n_total=n_total,
                                              n_samples=n_samples, m=m,
                                              model=basis_func, d=d)))
      } else {
        drm_opt_null <- do.call(optim, c(list(par=par_init_null_full,
                                              fn=negLDL_null),
                                         dot_args_null,
                                         list(g_null=g_null,
                                              par_pos=par_pos,
                                              par_dim=par_dim,
                                              x=x, n_total=n_total,
                                              n_samples=n_samples, m=m,
                                              model=basis_func, d=d)))
      }

      negldl_null <- drm_opt_null$value
      mele_null <- list(alpha=drm_opt_null$par[1:m],
                        gamma=drm_opt_null$par[-(1:m)])
      delr <- -2*(negldl - negldl_null)  # DEL ratio statistics
      df <- par_dim - par_dim_null_full
      p_val <- 1-pchisq(delr, df)

      return(list(drm_info=drm_info, mele=mele,
                  info_mat=info_mat, negldl=negldl,
                  mele_null=mele_null, negldl_null=negldl_null,
                  delr=delr, df=df, p_val=p_val,
                  p_est=p_est, cdf_est=cdf_est))

    } else {

      delr <- -2*negldl  # DEL ratio statistics
      df <- m*d
      p_val <- 1-pchisq(delr, df)

      return(list(drm_info=drm_info, mele=mele,
                  info_mat=info_mat, negldl=negldl,
                  delr=delr, df=df, p_val=p_val,
                  p_est=p_est, cdf_est=cdf_est))

    }

  }

}

summaryDRM <- function(drmfit)
{
  space3 <- "   "
  cat("Basic information about the DRM:\n")
  cat(space3, "Number of samples (m+1):", drmfit$drm_info$m+1, "\n")
  if (is.function(drmfit$drm_info$basis_func)) {
    cat(space3, "Basis function:\n")
    print(drmfit$drm_info$basis_func)
  } else {
    cat(space3, "Basis function:", drmfit$drm_info$basis_func, "\n")
  }
  cat(space3, "Dimension of the basis function (d):", drmfit$drm_info$d, "\n")
  cat(space3, "Sample sizes:", drmfit$drm_inf$n_samples, "\n")
  cat(space3, "Sample proportions (rho):",
      format(drmfit$drm_inf$rho, digits=3), "\n")

  cat("\n")
  mele_display <- displayPar(drmfit$mele, drmfit$drm_info$m)
  cat("Maximum empirical likelihood estimator (MELE):\n")
  print(format(mele_display, digits=3))

  cat("\n")
  if (is.null(drmfit$mele_null)) {
    cat("Default hypothesis testing problem:\n")
    cat(space3, "H_0: All distribution functions, F_k's, are equal.\n")
    cat(space3, "H_1: One of the distribution functions is different from the others.\n")
  }
  cat("Dual empirical likelihood ratio statistc (DELR):",
      format(drmfit$delr, nsmall=3), "\n")
  cat("Degree of freedom:", drmfit$df, "\n")
  cat("p-value:", format(drmfit$p_val, digits=3), "\n")

  cat("\n")
  summaryStatF <- summaryDRMFEst(drmfit$drm_info$n_samples, drmfit$p_est,
                                 drmfit$cdf_est, interpolation=TRUE)
  cat("Summary statistics of the estimated F_k's (mean, var -- variance, sd -- standard deviation, Q1 -- first quartile, Q3 -- third quartile, IQR -- inter-quartile range):\n")
  print(format(summaryStatF, digits=3))

}

trueParNormal <- function(m, mu, sigma) {
# m -- number of samples - 1
# mu -- mean vector of length m+1
# sigma -- standard deviation (sd) vector of length m+1

  par <- numeric(m*3)
  for (i in 1:m) {
    par[1 + 3*(i-1)] <- log(sigma[1]/sigma[i+1]) + (mu[1])^2/2/(sigma[1])^2 -
      (mu[i+1])^2/2/(sigma[i+1])^2  # alpha

    par[2 + 3*(i-1)] <- mu[i+1]/(sigma[i+1])^2 - mu[1]/(sigma[1])^2 # beta 1

    par[3 + 3*(i-1)] <- 1/2/(sigma[1])^2 - 1/2/(sigma[i+1])^2 # beta 2
  }

  return(par)

}

#trueParGamma <- function(m, shape, scale) {
trueParGamma <- function(m, shape, rate) {
# m -- number of samples - 1
# shape -- shape vector of length m+1
# scale -- scale vector of length m+1
# rate -- rate vector of length m+1
# gamma density:
  #f(x) = ( 1/(scale^shape * Gamma(shape)) ) * x^(shape) * exp(-x/scale)
  #or
  #f(x) = ( rate^shape / Gamma(shape) ) * x^(shape) * exp(-rate * x)

  par <- numeric(m*3)

  # scale formulation
  #for (i in 1:m) {
    #par[1 + 3*(i-1)] <- (shape[1] * log(scale[1]) + log(gamma(shape[1]))) -
    #(shape[i+1] * log(scale[i+1]) + log(gamma(shape[i+1])))  # alpha

    #par[2 + 3*(i-1)] <- 1/scale[1] - 1/scale[i+1] # beta 1

    #par[3 + 3*(i-1)] <- shape[i+1] - shape[1] # beta 2
  #}

  # rate formulation
  for (i in 1:m) {
    par[1 + 3*(i-1)] <- (shape[i+1] * log(rate[i+1]) - log(gamma(shape[i+1]))) -
      (shape[1] * log(rate[1]) - log(gamma(shape[1]))) # alpha

    par[2 + 3*(i-1)] <- rate[1] - rate[i+1] # beta 1

    par[3 + 3*(i-1)] <- shape[i+1] - shape[1] # beta 2
  }

  return(par)

}
