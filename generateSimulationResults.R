# Get distribution, n arguments from command line
args = commandArgs(trailingOnly = TRUE)

if(length(args) != 2){
  stop("Usage: Rscript generateSimulationResults n distribution")
}else if(!(as.integer(args[1]) %in% c(250, 500, 1000, 2500, 5000))){
  stop("n must be one of the following: 250, 500, 1000, 2500, 5000.")
}else if(!(args[2] %in% c("normal", "gamma", "lognormal", "test"))){
  stop("distribution must be one of the following: ")
}

# Run simulation study
source("R\\init.R")

source_data_folder = "simulation_data"
target_folder = "simulation_results"

# Hard coded parameters
d_sim = 5
model_sim = 12

# Parameters passed from command line
n_sim = as.integer(args[1])
distribution_sim = args[2]

# Read sim data
x_file_name = paste0(source_data_folder, "\\", distribution_sim, "_", n_sim, "_data.csv")
x_sim = as.matrix(read.csv(x_file_name))
# Remove column names
colnames(x_sim) = NULL

# Function to run simulation 
run_full_simulation = function(x, n, distribution){
  # Generate lambda values - multiply n by 3 because it represents size of only one sample
  lambda_vals_sim = ((n*3)^(1/3))*c(0, 0.05, 0.1, 0.5, 1, 2, 5, 10, 25, 50)
  
  print("Running regular simulation.")
  reg_results = runSimulation(x = x, n = n, model = model_sim, d = d_sim, 
                              lambdaVals = lambda_vals_sim)
  
  reg_file = paste0(target_folder, "\\", distribution, "_", n, "_reg_data.csv")
  reg_df = as.data.frame(reg_results)
  write.csv(reg_df, file=reg_file, row.names = FALSE)
  print("Regular file written to csv")
  
  print("Running adaptive simulation.")
  adap_results = runSimulation(x = x, n = n, model = model_sim, d = d_sim, 
                              lambdaVals = lambda_vals_sim, adaptive = TRUE)
  
  adap_file = paste0(target_folder, "\\", distribution, "_", n, "_adap_data.csv")
  adap_df = as.data.frame(adap_results)
  write.csv(adap_df, file=adap_file, row.names = FALSE)
  print("Adaptive file written to csv")
  
  print("Running aicbic simulation.")
  aicbic_results = simAICBIC(x = x, n = n)
  
  aicbic_file = paste0(target_folder, "\\", distribution, "_", n, "_aicbic_data.csv")
  aicbic_df = as.data.frame(aicbic_results)
  write.csv(aicbic_df, file=aicbic_file, row.names = FALSE)
  print("AICBIC file written to csv")
}

run_full_simulation(x = x_sim, n = n_sim, distribution = distribution_sim)