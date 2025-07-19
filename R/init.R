# Initialize the shared C library, R Scripts for drmdellasso
# This is assumed to be called from the directory root
# Changes to relative paths must be made if this is not the case

# Load .Renviron file to get shared lib name
readRenviron(".Renviron")
shared_lib = Sys.getenv("SHARED_LIB")

# Load dplyr library
library(dplyr)

# Build the shared library
setwd("src")
lib_str = paste0(shared_lib, ".dll")
command_str = paste("R CMD SHLIB -o",lib_str,"drmdelLasso.c utilities.c basisFuncs.c")
system(command_str)

# Load the shared library
dyn.load(lib_str)

# Load the R wrappers from the R folder
source("..\\R\\drmdelLasso.R")

# Define a wrapper function to run simulations for the default distributions and parameter values
# These are gamma and normal distributions
runSimulation = function(distribution, paramSetup, n, d, model, lambdaVals, adaptive = FALSE, runs = 1000){
  # x ,n_total, n_samples, m, d, model, lambda_vals, adaptive = FALSE
  
  # Ensure distribution is normal or gamma
  if(!(distribution %in% c("normal", "gamma", "lognormal"))){
    stop('Invalid distribution string passed. Currently only "normal", "lognormal" and "gamma" are supported.')
  }
  
  # Define parameters for the distributions
  if(distribution == "normal" || distribution == "lognormal"){
    if (paramSetup == 1){
      mu_0 = 2
      mu_1 = 1.8
      mu_2 = 1.9
      sigma_0 = 2
      sigma_1 = 0.89
      sigma_2 = 0.875
    } else if (paramSetup == 2){
      mu_0 = 5
      mu_1 = 4.5
      mu_2 = 5.5
      sigma_0 = 1.5
      sigma_1 = 1.25
      sigma_2 = 1
    } else if (paramSetup == 3){
      mu_0 = 18
      mu_1 = 17
      mu_2 = 19
      sigma_0 = 6
      sigma_1 = 6.5
      sigma_2 = 7
    } else {
      stop("An invalid code was passed to paramSetup. Currently only 1, 2, and 3 are supported.")
    }
  } else if(distribution == "gamma") { 
    if (paramSetup == 1){
      rate_0 = 2
      rate_1 = 1.4
      rate_2 = 1.2
      shape_0 = 1.8
      shape_1 = 1.2
      shape_2 = 1
    } else if (paramSetup == 2){
      rate_0 = 1
      rate_1 = 1.25
      rate_2 = 1.5
      shape_0 = 4
      shape_1 = 3
      shape_2 = 2
    } else if (paramSetup == 3){
      rate_0 = 1
      rate_1 = 0.9
      rate_2 = 0.8
      shape_0 = 6
      shape_1 = 7
      shape_2 = 8
    } else {
      stop("An invalid code was passed to paramSetup. Currently only 1, 2, and 3 are supported.")
    }
  }
  
  # Set the seed based on the sample size
  if(n == 50){
    set.seed(18)
  } else if (n == 100){
    set.seed(19)
  } else if(n == 250){
    set.seed(20)
  } else if (n == 500){
    set.seed(21)
  } else if (n == 1000){
    set.seed(22)
  } else if (n == 2500){
    set.seed(23)
  } else if (n == 5000){
    set.seed(24)
  } else {
    stop("An invalid sample size (individual) was provided. Currently only 250, 500, 1000, 2500 and 5000 are supported.")
  }
  
  # Set m - always setting to two for these simulation
  m = 2
  
  # Set n_total and n_samples
  n_samples = rep(n, m+1)
  n_total = sum(n_samples)
  
  # Compute path length, noting that 0 will be added if not passed
  # It will be added in the solution path function if not
  if(0 %in% lambdaVals){
    pathLength = length(lambdaVals)
  } else{
    pathLength = length(lambdaVals)+1
  }
  
  # Create matrix for simulation results
  pathLength = length(lambdaVals)
  simulationResults = matrix(0, nrow = runs*pathLength, ncol = 2*(d+1) + 6)
  
  # Populate matrix by generating solution paths
  for(i in 1:runs){
    # Simulate data
    if(distribution == "normal"){
      x0_test = rnorm(n,mu_0,sigma_0)
      x1_test = rnorm(n,mu_1,sigma_1)
      x2_test = rnorm(n,mu_2,sigma_2)
    } else if(distribution == "gamma") { 
      x0_test = rgamma(n,shape=shape_0,rate=rate_0)
      x1_test = rgamma(n,shape=shape_1,rate=rate_1)
      x2_test = rgamma(n,shape=shape_2,rate=rate_2)
    } else { # The distribution is log-normal
      x0_test = rlnorm(n,meanlog=mu_0,sdlog=sigma_0)
      x1_test = rlnorm(n,meanlog=mu_1,sdlog=sigma_1)
      x2_test = rlnorm(n,meanlog=mu_2,sdlog=sigma_2)
    }
    x_test = c(x0_test,x1_test,x2_test)
    
    # compute solution path
    solPath = solutionPath(x = x_test, n_total = n_total, n_samples = n_samples, 
                           m = m, d = d, model = model, lambdaVals, adaptive = adaptive)
    
    # Set the values of simulationResults
    simulationResults[(1+pathLength*(i-1)):(pathLength*i), 1] = i # Set first column to run iteration
    # Set the other values to solPath
    simulationResults[(1+pathLength*(i-1)):(pathLength*i), 2:(2*(d+1) + 6)] = solPath # essentially using m = 2
    
    # If we are on the first run, pull the column names from solPath
    colnames(simulationResults) = c("Run", colnames(solPath))
  }
  
  # Return simulationResults
  return(simulationResults)
}

# Define wrapper to return a sub-function of basis function 12 based on integer input
subBasisFunc = function(id) {
  result = switch(as.character(id),
   `1`  = 2, # code for log(x)
   `2`  = function(x) log(abs(x))^2,
   `3`  = function(x) c(log(abs(x)), log(abs(x))^2), # True basis function for lognormal
   `4`  = 3, # code for sqrt(x)
   `5`  = function(x) c(log(abs(x)), sqrt(abs(x))),
   `6`  = function(x) c(log(abs(x))^2, sqrt(abs(x))),
   `7`  = function(x) c(log(abs(x)), log(abs(x))^2, sqrt(abs(x))), # sub for lognormal
   `8`  = 1, # code for x
   `9`  = 6, # code for log(x), x - true basis function for gamma
   `10` = function(x) c(log(abs(x))^2, x),
   `11` = function(x) c(log(abs(x)), log(abs(x))^2, x), # sub for gamma, lognormal
   `12` = function(x) c(sqrt(abs(x)), x),
   `13` = 7, # code for log(x), sqrt(x), x, sub for gamma
   `14` = function(x) c(log(abs(x))^2, sqrt(abs(x)), x),
   `15` = function(x) c(log(abs(x)), log(abs(x))^2, sqrt(abs(x)), x), # sub for gamma, lognormal
   `16` = 4, # code for x^2
   `17` = function(x) c(log(abs(x)), x^2),
   `18` = function(x) c(log(abs(x))^2, x^2),
   `19` = function(x) c(log(abs(x)), log(abs(x))^2, x^2), # sub for lognormal
   `20` = function(x) c(sqrt(abs(x)), x^2),
   `21` = 8, # code for log(x), sqrt(x), x^2
   `22` = function(x) c(log(abs(x))^2, sqrt(abs(x)), x^2),
   `23` = function(x) c(log(abs(x)), log(abs(x))^2, sqrt(abs(x)), x^2), # sub for lognormal
   `24` = 5, # code for x, x^2 - true basis function for normal
   `25` = 9, # code for log(x), x, x^2, sub for gamma
   `26` = function(x) c(log(abs(x))^2, x, x^2),
   `27` = function(x) c(log(abs(x)), log(abs(x))^2, x, x^2),  # sub for gamma, lognormal
   `28` = 10, # code for srqt(x), x, x^2 
   `29` = 11, # code for log(x), sqrt(x), x, x^2, sub for gamma
   `30` = function(x) c(log(abs(x))^2, sqrt(abs(x)), x, x^2),
   `31` = 12  # code for log(x), log(x)^2, sqrt(x), x, x^2, sub for gamma, lognormal
   # To compute sub-basis proportion:
   # IDs 24-31 contain the normal basis function
   # IDs 9, 11, 13, 15, 25, 27, 29, 31 contain the gamma basis function
  )
  
  if (is.null(result)) stop("Invalid ID: must be 1 through 31")
  return(result)
}

# Define a wrapper to run simulations but based purely on the AIC/BIC of the MELEs of the sub-functions,
# as defined by Fokianos (2007)
simAICBIC = function(distribution, paramSetup, n, runs = 1000){
  # Ensure distribution is normal or gamma
  if(!(distribution %in% c("normal", "gamma", "lognormal"))){
    stop('Invalid distribution string passed. Currently only "normal", "lognormal" and "gamma" are supported.')
  }
  
  # Define parameters for the distributions
  if(distribution == "normal" || distribution == "lognormal"){
    if (paramSetup == 1){
      mu_0 = 2
      mu_1 = 1.8
      mu_2 = 1.9
      sigma_0 = 2
      sigma_1 = 0.89
      sigma_2 = 0.875
    } else if (paramSetup == 2){
      mu_0 = 5
      mu_1 = 4.5
      mu_2 = 5.5
      sigma_0 = 1.5
      sigma_1 = 1.25
      sigma_2 = 1
    } else if (paramSetup == 3){
      mu_0 = 18
      mu_1 = 17
      mu_2 = 19
      sigma_0 = 6
      sigma_1 = 6.5
      sigma_2 = 7
    } else {
      stop("An invalid code was passed to paramSetup. Currently only 1, 2, and 3 are supported.")
    }
  } else if(distribution == "gamma") { 
    if (paramSetup == 1){
      rate_0 = 2
      rate_1 = 1.4
      rate_2 = 1.2
      shape_0 = 1.8
      shape_1 = 1.2
      shape_2 = 1
    } else if (paramSetup == 2){
      rate_0 = 1
      rate_1 = 1.25
      rate_2 = 1.5
      shape_0 = 4
      shape_1 = 3
      shape_2 = 2
    } else if (paramSetup == 3){
      rate_0 = 1
      rate_1 = 0.9
      rate_2 = 0.8
      shape_0 = 6
      shape_1 = 7
      shape_2 = 8
    } else {
      stop("An invalid code was passed to paramSetup. Currently only 1, 2, and 3 are supported.")
    }
  }
  
  # Set the seed based on the sample size
  if(n == 50){
    set.seed(18)
  } else if (n == 100){
    set.seed(19)
  } else if(n == 250){
    set.seed(20)
  } else if (n == 500){
    set.seed(21)
  } else if (n == 1000){
    set.seed(22)
  } else if (n == 2500){
    set.seed(23)
  } else if (n == 5000){
    set.seed(24)
  } else {
    stop("An invalid sample size (individual) was provided. Currently only 250, 500, 1000, 2500 and 5000 are supported.")
  }
  
  # Set m - always setting to two for these simulation
  m = 2
  
  # Set n_total and n_samples
  n_samples = rep(n, m+1)
  n_total = sum(n_samples)
  
  # Setup matrix to store results
  # Only setup for basis function 12 right now - loop over all 31 options - so 31 rows per run
  # We will also only store columns with the AIC/BIC values since we are not currently interested in inference
  subFuncs = 31
  simulationResults = matrix(0, nrow = runs*subFuncs, ncol = 4)
    # Column 1: Run
    # Column 2: basisFuncID
    # Column 3: AIC
    # Column 4: BIC
  
  # Set columns names of the results
  colnames(simulationResults) = c("Run", "subFuncID", "AIC", "BIC")
  
  # Start runs
  for(i in 1:runs){
    # Generate random data
    # Simulate data
    if(distribution == "normal"){
      x0_test = rnorm(n,mu_0,sigma_0)
      x1_test = rnorm(n,mu_1,sigma_1)
      x2_test = rnorm(n,mu_2,sigma_2)
    } else if(distribution == "gamma") { 
      x0_test = rgamma(n,shape=shape_0,rate=rate_0)
      x1_test = rgamma(n,shape=shape_1,rate=rate_1)
      x2_test = rgamma(n,shape=shape_2,rate=rate_2)
    } else { # The distribution is log-normal
      x0_test = rlnorm(n,meanlog=mu_0,sdlog=sigma_0)
      x1_test = rlnorm(n,meanlog=mu_1,sdlog=sigma_1)
      x2_test = rlnorm(n,meanlog=mu_2,sdlog=sigma_2)
    }
    x_test = c(x0_test,x1_test,x2_test)
    
    # Iterate over each basis function
    for(j in 1:subFuncs){
      # Get sub function
      model = subBasisFunc(j)
      
      # Compute mele
      mele = drmdel(x=x_test, n_samples=n_samples, basis_func=model)$mele
      
      # Get d from mele
      d = (length(mele)-m)/m
      
      # Compute BIC and AIC
      AICBIC = aic_bic_drm(theta = mele, x = x_test, n_total = n_total, 
                           n_samples = n_samples, m = m, basis_func = model, d = d)
      
      # Set Values in matrix
      rowIdx = (i-1)*subFuncs + j
      
      simulationResults[rowIdx,1] = i
      simulationResults[rowIdx,2] = j
      simulationResults[rowIdx,3:4] = AICBIC
      
    }
  }
  
  # Return results
  return(simulationResults)
}

# Define a simple wrapper function that, for data from a simulation, computes:
#   1. The proportion of runs in the simulation that yielded at least one selection consistent solution
#   2. The proportion of runs in the simulation where AIC yielded a selection consistent solution
#   3. The proportion of runs in the simulation where AIC yielded a solution that contains a selection consistent solution
#   4. The proportion of runs in the simulation where BIC yielded a selection consistent solution
#   5. The proportion of runs in the simulation where BIC yielded a solution that contains a selection consistent solution
summariseSim = function(distribution, file_name, basis_func, tol){
  # distribution: (str) string specifying the distribution, must be "normal" or "gamma"
  # file_name: (str) the name of the file where the data is stored
  # basis_func: (int) the basis function used. Currently, only 12 is supported
  # tol: (double) tolerance used to check consistent solutions, values less than this tolerance will be considered 0
  
  # Check this is being run in the data folder
  if(sub("^.*/", "", getwd()) != "Data"){
    print("sumariseSim must be run from the Data folder!")
    return(-1)
  }
  
  # Ensure distribution string is valid
  if(!(distribution %in% c("normal", "gamma", "lognormal"))){
    print('Invalid distribution string passed. Currently only "normal", "lognormal and "gamma" are supported.')
    return(-1)
  }
  
  # Ensure the basis function is valid
  if(basis_func != 12){
    print('Invalid basis function passed. Currently 12 is supported.')
    return(-1)
  }
  
  # Import data, throwing an error if import fails
  simData = try(read.csv(file_name, header = TRUE), silent = TRUE)
  if(inherits(simData, "try-error")){
    print("File not found. Double check the specified file name.")
    return(-1)
  }
  
  # Get the number of runs
  numRuns = max(simData$Run)
  
  # Conditionally create string to select columns for basis function terms, depending on distribution
  # This may need to be changed conditionally if other distributions or other basis functions are used
  if(distribution == "gamma"){
    condStr = "1.$|4.$"
  } else if (distribution == "lognormal"){
    condStr = "1.$|2.$"
  } else{ # distribution is normal
    condStr = "4.$|5.$"
  }
  
  # Compute proportion of paths with at least one consitent solution
  consistSols <- simData %>%
    filter(
      if_all(
        select(., matches(condStr) & starts_with("beta")) %>% names(),
        ~ abs(.) > tol
      ),
      if_all(
        select(., !matches(condStr) & starts_with("beta")) %>% names(),
        ~ abs(.) <= tol
      )
    )
  
  consistSolsProp = nrow(distinct(consistSols, Run))/numRuns
  
  # Compute proportions for AIC
  # Get rows where AIC is the minimum
  minAIC = simData %>%
    group_by(Run) %>%
    slice(which.min(AIC)) %>%
    ungroup()
  
  # Compute proportion that are consistent
  minAICSols = minAIC %>%
    filter(
      if_all(
        select(., matches(condStr) & starts_with("beta")) %>% names(),
        ~ abs(.) > tol
      ),
      if_all(
        select(., !matches(condStr) & starts_with("beta")) %>% names(),
        ~ abs(.) <= tol
      )
    )
  
  consistMinAIC = nrow(minAICSols)/numRuns
  
  # Compute porportion that contain true basis function
  minAICSubs = minAIC %>%
    filter(
      if_all(
        select(., matches(condStr) & starts_with("beta")) %>% names(),
        ~ abs(.) > tol
      )
    )
  
  subMinAIC = nrow(minAICSubs)/numRuns
  
  # Compute proportions for BIC
  # Get rows where BIC is the minimum
  minBIC = simData %>%
    group_by(Run) %>%
    slice(which.min(BIC)) %>%
    ungroup()
  
  # Compute proportion that are consistent
  minBICSols = minBIC %>%
    filter(
      if_all(
        select(., matches(condStr) & starts_with("beta")) %>% names(),
        ~ abs(.) > tol
      ),
      if_all(
        select(., !matches(condStr) & starts_with("beta")) %>% names(),
        ~ abs(.) <= tol
      )
    )
  
  consistMinBIC = nrow(minBICSols)/numRuns
  
  # Compute porportion that contain true basis function
  minBICSubs = minBIC %>%
    filter(
      if_all(
        select(., matches(condStr) & starts_with("beta")) %>% names(),
        ~ abs(.) > tol
      )
    )
  
  subMinBIC = nrow(minBICSubs)/numRuns
  
  # return vector showing proportions
  return(c(runs = consistSolsProp, 
           AIC = consistMinAIC,
           AIC_sub = subMinAIC,
           BIC = consistMinBIC,
           BIC_sub = subMinBIC))
}

# Define a simple wrapper function that, for data from an AICBIC simulation, computes:
#   1. The proportion of runs in the simulation where AIC yielded a selection consistent solution
#   2. The proportion of runs in the simulation where AIC yielded a solution that contains a selection consistent solution
#   3. The proportion of runs in the simulation where BIC yielded a selection consistent solution
#   4. The proportion of runs in the simulation where BIC yielded a solution that contains a selection consistent solution
#   NOTE: This function currently assumes you are passing a simulation that considered all sub-basis functions of model 12
summariseAICBICSim = function(distribution, file_name){
  # distribution: (str) string specifying the distribution, must be "normal" or "gamma"
  # file_name: (str) the name of the file where the data is stored
  
  # Check this is being run in the data folder
  if(sub("^.*/", "", getwd()) != "Data"){
    print("sumariseSim must be run from the Data folder!")
    return(-1)
  }
  
  # Ensure distribution string is valid
  if(!(distribution %in% c("normal", "gamma", "lognormal"))){
    print('Invalid distribution string passed. Currently only "normal", "lognormal" and "gamma" are supported.')
    return(-1)
  }
  
  # Import data, throwing an error if import fails
  simData = try(read.csv(file_name, header = TRUE), silent = TRUE)
  if(inherits(simData, "try-error")){
    print("File not found. Double check the specified file name.")
    return(-1)
  }
  
  # Conditionally assign IDs for consistent solutions based on distribution
  # This may need to be changed conditionally if other distributions or other basis functions are used
  if(distribution == "gamma"){
    consistID = 9
    consistSubIDs = c(9, 11, 13, 15, 25, 27, 29, 31)
  } else if (distribution == "lognormal"){
    consistID = 3
    consistSubIDs = c(3, 7, 11, 15, 19, 23, 27, 31)
  } else{ # distribution is normal
    consistID = 24
    consistSubIDs = 24:31
  }
  
  # Get the number of runs
  numRuns = max(simData$Run)
  
  # Compute proportions for AIC
  # Get rows where AIC is the minimum
  minAIC = simData %>%
    group_by(Run) %>%
    slice(which.min(AIC)) %>%
    ungroup()
  
  # Number of consistent solutions
  minAICConsistRows = minAIC %>%
    filter(subFuncID == !!consistID)
  
  consistMinAIC = nrow(minAICConsistRows)/numRuns
  
  # Number of sub-consistent solutions
  minAICSubRows = minAIC %>%
    filter(subFuncID %in% !!consistSubIDs)
  
  subMinAIC = nrow(minAICSubRows)/numRuns
  
  # Compute proportions for BIC
  # Get rows where BIC is the minimum
  minBIC = simData %>%
    group_by(Run) %>%
    slice(which.min(BIC)) %>%
    ungroup()
  
  # Number of consistent solutions
  minBICConsistRows = minBIC %>%
    filter(subFuncID == !!consistID)
  
  consistMinBIC = nrow(minBICConsistRows)/numRuns
  
  # Number of sub-consistent solutions
  minBICSubRows = minBIC %>%
    filter(subFuncID %in% !!consistSubIDs)
  
  subMinBIC = nrow(minBICSubRows)/numRuns
  
  # Return vector showing proportions
  return(c(AIC = consistMinAIC,
           AIC_sub = subMinAIC,
           BIC = consistMinBIC,
           BIC_sub = subMinBIC))
}


