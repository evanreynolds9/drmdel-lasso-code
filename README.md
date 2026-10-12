# drmdel-lasso-code
This repository is for code and data accompanying the (currently unwritten) paper on basis function selection in DRMs with Lasso.

Much of the original code for the project was taken from the `drmdel` R package, authored by Dr. Song Cai. Lasso functionality was added by myself in 2024-2025.

## Getting Started

The code is currently not suitable for a production environment, but can be used cases by following the below instructions:

1. Clone the repository
1. Ensure you have a C compiler installed on your machine that is discoverable by R.
1. Confirm that the unit tests all pass on your machine by running the `testthat` tests in `unitTests.R`. The easiest way to do this is through RStudio.

## Running the simulations

To generate the simulation data, run the `generateSimulationData.R` script. This will populate the `simulation_data` folder with all necessary sample data.

To generate the simulation results by running the BCGD optimizer on both the regular and adaptive estimators, as well as running the comparable AIC/BIC method, use the 
```bash
RScript generateSimulationResults.R n distribution
```
and set n and distribution as arguments (they must be acceptable by the program).

These jobs may take a long time to run. If you would like to order multiple jobs on a server, you can do something like this:
```bash
for n in 250 500 1000 2500 5000
do
  for dist in normal gamma lognormal
  do
    nohup Rscript generateSimulationResults.R "$n" "$dist" &> "out_${n}_${dist}"
  done
done
```
