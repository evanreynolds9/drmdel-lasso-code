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
runSimulation = function(distribution, n, d, model, lambdaVals, adaptive = FALSE, runs = 1000){
  # x ,n_total, n_samples, m, d, model, lambda_vals, adaptive = FALSE
  
  # Ensure distribution is normal or gamma
  if(!(distribution %in% c("normal", "gamma"))){
    stop('Invalid distribution string passed. Currently only "normal" and "gamma" are supported.')
  }
  
  # Define parameters for the distributions
  if(distribution == "normal"){
    mu_0 = 2
    mu_1 = 1.8
    mu_2 = 1.9
    sigma_0 = 2
    sigma_1 = 0.89
    sigma_2 = 0.875
  } else if (distribution == "gamma"){
    rate_0 = 2
    rate_1 = 1.4
    rate_2 = 1.2
    shape_0 = 1.8
    shape_1 = 1.2
    shape_2 = 1
  }
  
  # Set the seed based on the sample size
  if(n == 250){
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
    } else{ # The distribution is gamma
      x0_test = rgamma(n,shape=shape_0,rate=rate_0)
      x1_test = rgamma(n,shape=shape_1,rate=rate_1)
      x2_test = rgamma(n,shape=shape_2,rate=rate_2)
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
  if(!(distribution %in% c("normal", "gamma"))){
    print('Invalid distribution string passed. Currently only "normal" and "gamma" are supported.')
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



