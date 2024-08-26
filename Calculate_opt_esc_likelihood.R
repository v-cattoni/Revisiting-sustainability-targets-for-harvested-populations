#calculating optimal escapement likelihood

# ----- reading in files -----
bh_r <- readRDS("./data_r_k/bh_r_likelihood.rds")
rk_r <- readRDS("./data_r_k/rk_r_likelihood.rds")
hs_r <- readRDS("./data_r_k/hs_r_likelihood.rds")

bh_k <- readRDS("./data_r_k/bh_k_likelihood.rds")
rk_k <- readRDS("./data_r_k/rk_k_likelihood.rds")
hs_k <- readRDS("./data_r_k/hs_k_likelihood.rds")

k_max <- readRDS("./data_r_k/Ri_max.rds")

# ----- constants -----
n_datasets = 284

# ----- intialising -----
bh_opt_esc_scaled_max_pop <- rk_opt_esc_scaled_max_pop <- hs_opt_esc_scaled_max_pop <- rep(NA, n_datasets)
bh_opt_esc_scaled_k1 <- rk_opt_esc_scaled_k1 <- hs_opt_esc_scaled_k1 <- rep(NA, n_datasets)

# ----- functions -----
#Beverton-Holt optimal escapement
bh_opt_esc <- function(r, k) {
  return(k * (r ** 0.5 - 1) / (r - 1))
}

#Hockey-Stick optimal escpaement
hs_opt_esc <- function(r, k) {
  return(k / r)
}

#Ricker optimal escapement
# Define the function to find the root of
derivative_rk <- function(X) {
  # Calculate the derivative of ri with respect to X
  derivative <-
    (p ^ (1 - (X / k))) * (k - X * log(p)) / k - 1

  # Return the derivative
  return(derivative)
}


for(i in 1:n_datasets){
  bh_opt_esc_scaled_max_pop[i] <- bh_opt_esc(bh_r[i], bh_k[i])
  hs_opt_esc_scaled_max_pop[i] <- hs_opt_esc(hs_r[i], min(1, hs_k[i]))
  p = rk_r[i]
  k = rk_k[i]
  result = uniroot(derivative_rk, lower = 0, upper = 10*k)
  rk_opt_esc_scaled_max_pop[i] <- result$root
  
  
  bh_opt_esc_scaled_k1[i] <- bh_opt_esc(bh_r[i], 1)
  hs_opt_esc_scaled_k1[i] <- hs_opt_esc(hs_r[i], 1)
  p = rk_r[i]
  k = 1
  result = uniroot(derivative_rk, lower = 0, upper = 10*k)
  rk_opt_esc_scaled_k1[i] <- result$root
}

saveRDS(bh_opt_esc_scaled_max_pop, file = "./data_opt_esc/bh_opt_esc_prop_max_pop_likelihood.rds")
saveRDS(rk_opt_esc_scaled_max_pop, file = "./data_opt_esc/rk_opt_esc_prop_max_pop_likelihood.rds")
saveRDS(hs_opt_esc_scaled_max_pop, file = "./data_opt_esc/hs_opt_esc_prop_max_pop_likelihood.rds")

saveRDS(bh_opt_esc_scaled_k1, file = "./data_opt_esc/bh_opt_esc_prop_fitted_k_likelihood.rds")
saveRDS(rk_opt_esc_scaled_k1, file = "./data_opt_esc/rk_opt_esc_prop_fitted_k_likelihood.rds")
saveRDS(hs_opt_esc_scaled_k1, file = "./data_opt_esc/hs_opt_esc_prop_fitted_k_likelihood.rds")