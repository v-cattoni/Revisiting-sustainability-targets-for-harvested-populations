# Calculate Optimal Escapement Bayes

# Reading in files ####
Ri <- readRDS("./data_files/Ri.rds")
Ei <- readRDS("./data_files/Ei.rds")

# Constants ####
years <- seq(1950, 2020, 1) # year of each row of data
n_rows <- length(years) # num rows to iterate through all data
n_datasets <- length(Ri) # number of dataets = number of cols 
n_iter = 10000 # number of iterations used in mcmc model
n_branch = 4 # number of branches used in mcmc model

# Optimal Escapement Functions ####
#Beverton-Holt optimal escapement
bh_opt_esc_f <- function(r, k) {
  return(k * (r ** 0.5 - 1) / (r - 1))
}

#Hockey-Stick optimal escpaement
hs_opt_esc_f <- function(r, k) {
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


# Initialising ####
bio_mat <- matrix(NA, n_rows, n_datasets) # biomass (MT)
catch_mat <- matrix(NA, n_rows, n_datasets) # catch (MT)
Ri_max <- rep(NaN,n_datasets) # max Ri (MT)
bh_r_bayes <- rk_r_bayes  <- hs_r_bayes <- rep(NaN,n_datasets*n_iter*n_branch) #growth rate
bh_k_bayes <- rk_k_bayes  <- hs_k_bayes <- rep(NaN,n_datasets*n_iter*n_branch) #carrying capacity
bh_opt_esc_prop_fitted_k_bayes <- rk_opt_esc_prop_fitted_k_bayes  <- hs_opt_esc_prop_fitted_k_bayes <- rep(NaN,n_datasets*n_iter*n_branch) #carrying capacity prop of fitted k
bh_opt_esc_prop_max_pop_bayes <- rk_opt_esc_prop_max_pop_bayes  <- hs_opt_esc_prop_max_pop_bayes <- rep(NaN,n_datasets*n_iter*n_branch) #carrying capacity prop of max pop


# Storing max biomass for each species ####
for (i in 1:length(Ri)) {
  Ri_max[i] <- max(Ri[[i]], Ei[[i]])
}

# Main ####
# storing growth rate and carrying capacity by dataset, model, branch, iteration
for (dataset_id in 1:n_datasets){
  print(dataset_id)
  load(sprintf("./mcmc_data/results_dataset%d.RData",dataset_id))
  for (branch in 1:n_branch){
    for (iter in 1:n_iter){
      bh_r_bayes[(dataset_id-1)*n_iter*n_branch + (branch - 1)*n_iter + iter] <- bh_data@sim$samples[[branch]]$r[iter]
      bh_k_bayes[(dataset_id-1)*n_iter*n_branch + (branch - 1)*n_iter + iter] <- bh_data@sim$samples[[branch]]$k[iter]/Ri_max[dataset_id]
      hs_r_bayes[(dataset_id-1)*n_iter*n_branch + (branch - 1)*n_iter + iter] <- hs_data@sim$samples[[branch]]$r[iter]
      hs_k_bayes[(dataset_id-1)*n_iter*n_branch + (branch - 1)*n_iter + iter] <- hs_data@sim$samples[[branch]]$k[iter]/Ri_max[dataset_id]
      rk_r_bayes[(dataset_id-1)*n_iter*n_branch + (branch - 1)*n_iter + iter] <- rk_data@sim$samples[[branch]]$r[iter]
      rk_k_bayes[(dataset_id-1)*n_iter*n_branch + (branch - 1)*n_iter + iter] <- rk_data@sim$samples[[branch]]$k[iter]/Ri_max[dataset_id]
    }
  }
}

# # Storing optimal escapement as a proportion of the fitted K ####
for (i in 1:(n_datasets*n_branch*n_iter)){
  #print(i)
  bh_opt_esc_prop_fitted_k_bayes[i] <- bh_opt_esc_f(bh_r_bayes[i], 1)
  hs_opt_esc_prop_fitted_k_bayes[i] <- hs_opt_esc_f(hs_r_bayes[i], 1)
  p <- rk_r_bayes[i]
  k <- 1
  result <- uniroot(derivative_rk, lower = 0, upper = 10 * k)
  rk_opt_esc_prop_fitted_k_bayes[i] <- result$root
}

# Storing optimal escapement as a proportion of the max population ####
for (i in 1:(n_datasets*n_branch*n_iter)){
  bh_opt_esc_prop_max_pop_bayes[i] <- bh_opt_esc_f(bh_r_bayes[i], bh_k_bayes[i])
  hs_opt_esc_prop_max_pop_bayes[i] <- hs_opt_esc_f(hs_r_bayes[i], hs_k_bayes[i])
  p <- rk_r_bayes[i]
  k <- rk_k_bayes[i]
  result <- uniroot(derivative_rk, lower = 0, upper = 10 * k)
  rk_opt_esc_prop_max_pop_bayes[i] <- result$root
}

# Saving files ####
#saving opt esc data
saveRDS(bh_r_bayes, file = "./data_r_k/bh_r_bayes.rds")
saveRDS(hs_r_bayes, file = "./data_r_k/hs_r_bayes.rds")
saveRDS(rk_r_bayes, file = "./data_r_k/rk_r_bayes.rds")

saveRDS(bh_k_bayes, file = "./data_r_k/bh_k_bayes.rds")
saveRDS(hs_k_bayes, file = "./data_r_k/hs_k_bayes.rds")
saveRDS(rk_k_bayes, file = "./data_r_k/rk_k_bayes.rds")


saveRDS(bh_opt_esc_prop_fitted_k_bayes, file = "./opt_esc_data/bh_opt_esc_prop_fitted_k_bayes.rds")
saveRDS(hs_opt_esc_prop_fitted_k_bayes, file = "./opt_esc_data/hs_opt_esc_prop_fitted_k_bayes.rds")
saveRDS(rk_opt_esc_prop_fitted_k_bayes, file = "./opt_esc_data/rk_opt_esc_prop_fitted_k_bayes.rds")

saveRDS(bh_opt_esc_prop_max_pop_bayes, file = "./opt_esc_data/bh_opt_esc_prop_max_pop_bayes.rds")
saveRDS(hs_opt_esc_prop_max_pop_bayes, file = "./opt_esc_data/hs_opt_esc_prop_max_pop_bayes.rds")
saveRDS(rk_opt_esc_prop_max_pop_bayes, file = "./opt_esc_data/rk_opt_esc_prop_max_pop_bayes.rds")

