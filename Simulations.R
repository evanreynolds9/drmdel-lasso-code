# Run simulations

# Call init.R to setup shared library and all functions
# This should be run with the working directory as the project root
# It will not change the working directory
source("R\\init.R")

# Set working directory to Data folder - all simulation results will be saved here
# Also, other functions will pull data from here to compute summary statistics and must be called from this folder
setwd("..")
setwd("Data")

# Set parameter values
distribution = "gamma"
paramSetup = 3
n = 1000
model = 12
lambdaVals = c(0, ((n*3)^(1/3))*c(0.1,0.25,0.5,1,2,5,10,25,50))
adaptive = TRUE
runs = 50

# Compute d based on model
if(model %in% 1:4){
  d = 1
}else if(model %in% 5:6){
  d = 2
}else if(model %in% 7:10){
  d = 3
}else if(model == 11){
  d = 4
}else if(model == 12){
  d = 5
}else{
  stop("Invalid model provided - adjust the program header.")
}

# Create string for estimator type
if(adaptive){
  adapStr = "adap"
} else{
  adapStr = "reg"
}

### lasso-drm simulation ###
# Run simulation
simData = runSimulation(distribution, paramSetup, n, d, model, lambdaVals, adaptive = adaptive, runs = runs)

# Write data
simDataDf = as.data.frame(simData)

fileName = paste0(paste("sim_results", distribution, paramSetup, adapStr, n, sep = "_"), ".csv")
write.csv(simDataDf, fileName, row.names = FALSE)

# Get proportions from simulations
simSum = summariseSim(distribution, fileName, basis_func = model, tol = 1e-12)
print(simSum)

### AIC/BIC simulation ###
# Run the same simulation using the AIC/BIC method
AICBICSimData = simAICBIC(distribution, paramSetup, n, runs = runs)
AICBICSimDataDf = as.data.frame(AICBICSimData)

fileNameAICBIC = paste0(paste("sim_results_AICBIC", distribution, paramSetup, n, sep = "_"), ".csv")
write.csv(AICBICSimDataDf, fileNameAICBIC, row.names = FALSE)

# Get porportions for AIC/BIC simulation
AICBICSum = summariseAICBICSim(distribution, fileNameAICBIC)
print(AICBICSum)


# Check results for multiple simulations, if desired:

checkAllSims(3, "normal", 1000)
checkAllSims(3,"gamma", 1000)
