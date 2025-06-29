# drmdel-lasso-code
This repository is for code and data accompanying the (currently unwritten) paper on basis function selection in DRMs with Lasso.

Much of the original code for the project was taken from the `drmdel` R package, authored by Dr. Song Cai. Lasso functionality was added by myself in 2024-2025.

## Getting Started

The code is currently not suitable for a production environment, but can be used for suitable cases by following the below instructions:

1. Ensure you have a C compiler installed on your machine.
2. Create a .Renviron file in root directory, and set the name you choose for your shared library to the following variable: `SHARED_LIB`
3. In your R session, from the project root, run the `init.R` file to build the shared library and load all functions. This can be done with `source("R\\init.R")`.

`The `Simulations.R` in the simulation gives an example of the usage.