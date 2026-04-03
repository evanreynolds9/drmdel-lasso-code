# Run small simulations locally to decide which parameter setup to run on server for each distn

# Call init.R to setup shared library and all functions
# This should be run with the working directory as the project root
# It will not change the working directory
source("R\\init.R")

# Set working directory to Data folder - all simulation results will be saved here
# Also, other functions will pull data from here to compute summary statistics and must be called from this folder
setwd("Data")

# Set parameter values
n = 50
model = 12
lambdaVals = c(0, ((n*3)^(1/3))*c(0.01,0.05,0.1,0.25,0.5,1,2,5,10))
runs = 100

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

distributions = c("normal", "gamma", "lognormal")
estTypes = c("reg", "adap")

# Run through distributions and parameter setup, computing a regular and adaptive simulation for each
for(distribution in distributions){
  for(paramSetup in 1:3){
    for(estType in estTypes){
      if(estType == "reg"){
        adaptive = FALSE
      }else{
        adaptive = TRUE
      }
      simData = runSimulation(distribution=distribution, paramSetup=paramSetup, n=n, model=model, d=d, 
                              lambdaVals=lambdaVals, adaptive=adaptive, runs=runs)
      simDataDf = as.data.frame(simData)
      fileName = paste0(paste("local_sim_results", distribution, paramSetup, estType, n, sep = "_"), ".csv")
      write.csv(simDataDf, fileName, row.names = FALSE)
      
      # Summarise
      simSum = summariseSim(distribution=distribution, file_name=fileName, basis_func=model, tol=1e-12)
      print("Results for:")
      print(paste("   Distribution:", distribution))
      print(paste("   ParameterSet:", paramSetup))
      print(paste("   EstimateType:", estType))
      print(simSum)
    }
  }
}