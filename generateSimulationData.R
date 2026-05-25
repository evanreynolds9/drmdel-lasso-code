# Create simulation data for simulation study
source("R\\init.R")

target_folder = "simulation_data"

# create function to return a matrix with the simulated data
generateSimulationData = function(runs, distribution, paramSetup, n){
  simMatrix = matrix(0, nrow=runs, ncol=n*3) # 3 is hard coded because all simulations use m=2
  generator = createRandomGenerator(distribution=distribution, paramSetup=paramSetup, n=n) # will set seed
  for(i in 1:runs){
    x = generator()
    simMatrix[i,] = x
  }
  return(simMatrix)
}

# Use second paramsetup for every distribution
paramSetup = 2
runs = 1000

# Simulate data for each distribution
for(distribution in c("normal", "gamma", "lognormal")){
  for(n in c(250, 500, 1000, 2500, 5000)){
    fileName = paste0(target_folder, "\\", distribution, "_", n, "_data.csv")
    print(fileName)
    simMatrixN = generateSimulationData(runs=runs, distribution=distribution, paramSetup=paramSetup, n=n)
    simDataDf = as.data.frame(simMatrixN)
    write.csv(simDataDf, file=fileName, row.names = FALSE)
  }
}
