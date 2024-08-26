#calculating weighted aic for all datasets

n_datasets = 284
bh_ml <- rk_ml <- hs_ml <- rep(NA, n_datasets) #maximum likelihood
bh_r <- rk_r <- hs_r <- rep(NA, n_datasets) #growth rate
bh_k <- rk_k <- hs_k <- rep(NA, n_datasets) # carrying capacity
bh_k_scaled <- rk_k_scaled <- hs_k_scaled <- rep(NA, n_datasets) # carrying capacity scaled
n_datapoints <- rep(NA, n_datasets)
bh_sigma_y <- rk_sigma_y <- hs_sigma_y <- rep(NA, n_datasets) #sigma_y
bh_r_mcmc_mean <- rk_r_mcmc_mean <- hs_r_mcmc_mean <- rep(NA, n_datasets)
bh_r_mcmc_median <- rk_r_mcmc_median <- hs_r_mcmc_median <- rep(NA, n_datasets)
k_max <- rep(NA, n_datasets) # stores max pop
Ri_data <- readRDS("./data_files/Ri.rds") #recruitment
Ei_data <- readRDS("./data_files/Ei.rds") #escapement

loglike_bh <- function(params){ 
  
  # Make sure the ordering you're using matches what you give for the initial values
  r <- params[1]
  k <- params[2]
  sigma_y <- params[3]
  
  # The log likelihood
  val <- sum(dlnorm(Ri, log(Ei) + log(r) - log(1 + (r - 1) * Ei / k), sigma_y, log = TRUE))
  
  # Return negative log likelihood because we're minimising
  return(-1*val) 
}

loglike_rk <- function(params){ 
  
  # Make sure the ordering you're using matches what you give for the initial values
  r <- params[1]
  k <- params[2]
  sigma_y <- params[3]
  
  # The log likelihood
  val <- sum(dlnorm(Ri, log(Ei) + (1- Ei / k) * log(r), sigma_y, log = TRUE))
  
  # Return negative log likelihood because we're minimising
  return(-1*val) 
}

loglike_hs <- function(params){ 
  
  # Make sure the ordering you're using matches what you give for the initial values
  r <- params[1]
  k <- params[2]
  sigma_y <- params[3]
  
  # The log likelihood
  val <- sum(dlnorm(Ri, log(pmin(r * Ei, k)), sigma_y, log = TRUE))
  
  # Return negative log likelihood because we're minimising
  return(-1*val) 
}

for (i in 1:n_datasets){
  Ri <- Ri_data[[i]]
  Ei <- Ei_data[[i]]
  k_max[i] <- max(Ri)
  print(i)
  load(sprintf("./mcmc_data/results_dataset%d.RData",i))
  
  bh_r_mcmc <- c(bh_data@sim$samples[[1]]$r, bh_data@sim$samples[[2]]$r, bh_data@sim$samples[[3]]$r, bh_data@sim$samples[[4]]$r)
  bh_k_mcmc <- c(bh_data@sim$samples[[1]]$k, bh_data@sim$samples[[2]]$k, bh_data@sim$samples[[3]]$k, bh_data@sim$samples[[4]]$k)
  bh_sigma_y_mcmc <- c(bh_data@sim$samples[[1]]$sigma_y, bh_data@sim$samples[[2]]$sigma_y, bh_data@sim$samples[[3]]$sigma_y, bh_data@sim$samples[[4]]$sigma_y)
  rk_r_mcmc <- c(rk_data@sim$samples[[1]]$r, rk_data@sim$samples[[2]]$r, rk_data@sim$samples[[3]]$r, rk_data@sim$samples[[4]]$r)
  rk_k_mcmc <- c(rk_data@sim$samples[[1]]$k, rk_data@sim$samples[[2]]$k, rk_data@sim$samples[[3]]$k, rk_data@sim$samples[[4]]$k)
  rk_sigma_y_mcmc <- c(rk_data@sim$samples[[1]]$sigma_y, rk_data@sim$samples[[2]]$sigma_y, rk_data@sim$samples[[3]]$sigma_y, rk_data@sim$samples[[4]]$sigma_y)
  hs_r_mcmc <- c(hs_data@sim$samples[[1]]$r, hs_data@sim$samples[[2]]$r, hs_data@sim$samples[[3]]$r, hs_data@sim$samples[[4]]$r)
  hs_k_mcmc <- c(hs_data@sim$samples[[1]]$k, hs_data@sim$samples[[2]]$k, hs_data@sim$samples[[3]]$k, hs_data@sim$samples[[4]]$k)
  hs_sigma_y_mcmc <- c(hs_data@sim$samples[[1]]$sigma_y, hs_data@sim$samples[[2]]$sigma_y, hs_data@sim$samples[[3]]$sigma_y, hs_data@sim$samples[[4]]$sigma_y)
  
  bh_r_mcmc_mean[i] <- mean(bh_r_mcmc)
  rk_r_mcmc_mean[i] <- mean(rk_r_mcmc)
  hs_r_mcmc_mean[i] <- mean(hs_r_mcmc)
  
  bh_r_mcmc_median[i] <- median(bh_r_mcmc)
  rk_r_mcmc_median[i] <- median(rk_r_mcmc)
  hs_r_mcmc_median[i] <- median(hs_r_mcmc)
  
  bh_r_init <- bh_r_mcmc[1]
  bh_k_init <- bh_k_mcmc[1]
  bh_sigma_y_init <- bh_sigma_y_mcmc[1]
  
  rk_r_init <- bh_r_mcmc[1]
  rk_k_init <- bh_k_mcmc[1]
  rk_sigma_y_init <- bh_sigma_y_mcmc[1]
  hs_r_init <- bh_r_mcmc[1]
  hs_k_init <- bh_k_mcmc[1]
  hs_sigma_y_init <- bh_sigma_y_mcmc[1]
  
  bh_ml_init <- Inf
  rk_ml_init <- Inf
  hs_ml_init <- Inf
  
  
  for (j in 1:length(bh_r_mcmc)){
    params <- c(bh_r_mcmc[j], bh_k_mcmc[j], bh_sigma_y_mcmc[j])
    #if (!is.na(loglike_bh(params)) && loglike_bh(params) < bh_ml_init && bh_r_mcmc[j] > 1.00000000001 && bh_k_mcmc[j] > 0.1*max(Ri)) {
    if (loglike_bh(params) < bh_ml_init){
      bh_r_init <- bh_r_mcmc[j]
      bh_k_init <- bh_k_mcmc[j]
      bh_sigma_y_init <- bh_sigma_y_mcmc[j]
      bh_ml_init <- loglike_bh(params)
    }
    
    params <- c(rk_r_mcmc[j], rk_k_mcmc[j], rk_sigma_y_mcmc[j])
    #if (!is.na(loglike_rk(params)) && loglike_rk(params) < rk_ml_init && rk_r_mcmc[j] > 1.00000000001 && rk_k_mcmc[j] > 0.1*max(Ri)) {
    if (loglike_rk(params) < rk_ml_init){  
      rk_ml_init <- loglike_rk(params)
      rk_r_init <- rk_r_mcmc[j]
      rk_k_init <- rk_k_mcmc[j]
      rk_sigma_y_init <- rk_sigma_y_mcmc[j]
    }
    
    params <- c(hs_r_mcmc[j], hs_k_mcmc[j], hs_sigma_y_mcmc[j])
    #if (!is.na(loglike_hs(params)) && loglike_hs(params) < hs_ml_init && hs_r_mcmc[j] > 1.00000000001 && hs_k_mcmc[j] > 0.1*max(Ri)) {
    if (loglike_hs(params) < hs_ml_init){ 
      hs_ml_init <- loglike_hs(params)
      hs_r_init <- hs_r_mcmc[j]
      hs_k_init <- hs_k_mcmc[j]
      hs_sigma_y_init <- hs_sigma_y_mcmc[j]
    }
    
  }
  
  
  # calculating maximum likelihood
  params_bh <- c(bh_r_init, bh_k_init, bh_sigma_y_init)
  params_rk <- c(rk_r_init, rk_k_init, rk_sigma_y_init)
  params_hs <- c(hs_r_init, hs_k_init, hs_sigma_y_init)
  
  
  #storing results
  #this is where i am trying to bound the growth rate and carrying capacity in the optimisation
  bh_opt <- hjkb(fn = loglike_bh, par = params_bh, lower = c(1.00000000001, 0.1*max(Ri), 0), upper = c(10, 10*max(Ri), Inf))
  bh_ml[i] <- bh_opt$value
  bh_r[i] <- bh_opt$par[1]
  bh_k[i] <- bh_opt$par[2]
  bh_k_scaled[i] <- bh_opt$par[2]/max(Ri)
  bh_sigma_y[i] <- bh_opt$par[3]
  
  rk_opt <- hjkb(fn = loglike_rk, par = params_rk, lower = c(1.00000000001, 0.1*max(Ri), 0), upper = c(10, 10*max(Ri), Inf))
  rk_ml[i] <- rk_opt$value
  rk_r[i] <- rk_opt$par[1]
  rk_k[i] <- rk_opt$par[2]
  rk_k_scaled[i] <- rk_opt$par[2]/max(Ri)
  rk_sigma_y[i] <- rk_opt$par[3]
  
  hs_opt <- hjkb(fn = loglike_hs, par = params_hs, lower = c(1.00000000001, 0.1*max(Ri), 0), upper = c(10, 10*max(Ri), Inf))
  hs_ml[i] <- hs_opt$value
  hs_r[i] <- hs_opt$par[1]
  hs_k[i] <- hs_opt$par[2]
  hs_k_scaled[i] <- hs_opt$par[2]/max(Ri)
  hs_sigma_y[i] <- hs_opt$par[3]
  
  n_datapoints[i] <- length(Ri)
}

# calculating aic

bh_aic <- 2*bh_ml
rk_aic <- 2*rk_ml
hs_aic <- 2*hs_ml

min_aic <- rep(NA, n_datasets)

for(i in 1:n_datasets){
  min_aic[i] <- min(c(bh_aic[i], rk_aic[i], hs_aic[i]))
}

bh_delta_aic <- bh_aic - min_aic
rk_delta_aic <- rk_aic - min_aic
hs_delta_aic <- hs_aic - min_aic

sum_neg_half_delta_aic <- exp((-1/2)*bh_delta_aic) + exp((-1/2)*rk_delta_aic) + exp((-1/2)*hs_delta_aic)

bh_waic <- exp((-1/2)*bh_delta_aic)/sum_neg_half_delta_aic
rk_waic <- exp((-1/2)*rk_delta_aic)/sum_neg_half_delta_aic
hs_waic <- exp((-1/2)*hs_delta_aic)/sum_neg_half_delta_aic

saveRDS(bh_waic, file = "./data_waic/bh_waic.rds")
saveRDS(rk_waic, file = "./data_waic/rk_waic.rds")
saveRDS(hs_waic, file = "./data_waic/hs_waic.rds")

saveRDS(n_datapoints, file = "./data_n_datapoints/n_datapoints.rds")
saveRDS(order(n_datapoints, decreasing = FALSE), file = "./data_n_datapoints/n_datapoints_increasing.rds")

saveRDS(bh_r, file = "./data_r_k/bh_r_likelihood.rds")
saveRDS(rk_r, file = "./data_r_k/rk_r_likelihood.rds")
saveRDS(hs_r, file = "./data_r_k/hs_r_likelihood.rds")

saveRDS(bh_k_scaled, file = "./data_r_k/bh_k_likelihood.rds")
saveRDS(rk_k_scaled, file = "./data_r_k/rk_k_likelihood.rds")
saveRDS(hs_k_scaled, file = "./data_r_k/hs_k_likelihood.rds")

saveRDS(k_max, file = "./data_r_k/Ri_max.rds")
