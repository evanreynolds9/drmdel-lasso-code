# Quick script to give the non-zero coefficient values of the DRM parameters

#normal/lognormal
mu_0 = 2
mu_1 = 1.8
mu_2 = 1.9
sigma_0 = 2
sigma_1 = 0.89
sigma_2 = 0.875
#gamma
rate_0 = 1
rate_1 = 1.25
rate_2 = 1.5
shape_0 = 4
shape_1 = 3
shape_2 = 2

#normal parameters
b15_norm = mu_1/sigma_1^2 - mu_0/sigma_0^2
b17_norm = 1/(2*sigma_0^2) - 1/(2*sigma_1^2)
b25_norm = mu_2/sigma_2^2 - mu_0/sigma_0^2
b27_norm = 1/(2*sigma_0^2) - 1/(2*sigma_2^2)
round(c(b15_norm,b17_norm),2)
round(c(b25_norm,b27_norm),2)

#gamma
b11_gamma = shape_1 - shape_0
b15_gamma = rate_0 - rate_1
b21_gamma = shape_2 - shape_0
b25_gamma = rate_0 - rate_2
round(c(b11_gamma,b15_gamma),2)
round(c(b21_gamma,b25_gamma),2)

#lognormal
b11_lognorm = b15_norm
b13_lognorm = b17_norm
b21_lognorm = b25_norm
b23_lognorm = b27_norm
round(c(b11_lognorm,b13_lognorm), 2)
round(c(b21_lognorm,b23_lognorm), 2)


