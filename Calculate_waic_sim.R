#Calculate weighted AIC for all datasets for simulation

n_sims <- 20
bh_ml_bh <- rk_ml_bh <- hs_ml_bh <- rep(NA, n_sims) #maximum likelihood when data is generated from bh
bh_ml_rk <- rk_ml_rk <- hs_ml_rk <- rep(NA, n_sims) #maximum likelihood when data is generated from rk
bh_ml_hs <- rk_ml_hs <- hs_ml_hs <- rep(NA, n_sims) #maximum likelihood when data is generated from hs
bh_r_bh <- rk_r_bh <- hs_r_bh <- rep(NA, n_sims) #growth rate when data is generated from bh
bh_r_rk <- rk_r_rk <- hs_r_rk <- rep(NA, n_sims) #growth rate when data is generated from rk
bh_r_hs <- rk_r_hs <- hs_r_hs <- rep(NA, n_sims) #growth rate when data is generated from hs
bh_k_bh <- rk_k_bh <- hs_k_bh <- rep(NA, n_sims) # carrying capacity when data is generated from bh
bh_k_rk <- rk_k_rk <- hs_k_rk <- rep(NA, n_sims) # carrying capacity when data is generated from rk
bh_k_hs <- rk_k_hs <- hs_k_hs <- rep(NA, n_sims) # carrying capacity when data is generated from hs
bh_k_scaled_bh <- rk_k_scaled_bh <- hs_k_scaled_bh <- rep(NA, n_sims) # carrying capacity scaled when data is generated from bh
bh_k_scaled_rk <- rk_k_scaled_rk <- hs_k_scaled_rk <- rep(NA, n_sims) # carrying capacity scaled when data is generated from rk
bh_k_scaled_hs <- rk_k_scaled_hs <- hs_k_scaled_hs <- rep(NA, n_sims) # carrying capacity scaled when data is generated from hs
bh_sigma_y_bh <- rk_sigma_y_bh <- hs_sigma_y_bh <- rep(NA, n_sims) #sigma_y when data is generated from bh
bh_sigma_y_rk <- rk_sigma_y_rk <- hs_sigma_y_rk <- rep(NA, n_sims) #sigma_y when data is generated from rk 
bh_sigma_y_hs <- rk_sigma_y_hs <- hs_sigma_y_hs <- rep(NA, n_sims) #sigma_y when data is generated from hs
Ri_data_bh <- Ri_data_rk <- Ri_data_hs <- rep(NA, n_sims)
Ei_data_bh <- Ei_data_rk <- Ei_data_hs <- rep(NA, n_sims)
bh_n_datapoints <- rk_n_datapoints <- hs_n_datapoints <- rep(NA, n_sims)

for(i in 1:20){
  load(sprintf("./data_files_sim/sim_data_Matt_repeat%g.RData", i))
  Ri_data_bh[i] <- list(sim_bh$Ri)
  Ei_data_bh[i] <- list(sim_bh$Ei)
  Ri_data_rk[i] <- list(sim_rk$Ri)
  Ei_data_rk[i] <- list(sim_rk$Ei)
  Ri_data_hs[i] <- list(sim_hs$Ri)
  Ei_data_hs[i] <- list(sim_hs$Ei)
}

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

for (i in 1:20){
  print(i)
  for(model in c("bh", "rk", "hs")){
    print(model)
    
    if(model == "bh"){
      Ri <- Ri_data_bh[[i]]
      Ei <- Ei_data_bh[[i]]
      bh_n_datapoints[i] <- length(Ei)
    }
    
    if(model == "rk"){
      Ri <- Ri_data_rk[[i]]
      Ei <- Ei_data_rk[[i]]
      rk_n_datapoints[i] <- length(Ei)
    }
    
    if(model == "hs"){
      Ri <- Ri_data_hs[[i]]
      Ei <- Ei_data_hs[[i]]
      hs_n_datapoints[i] <- length(Ei)
    }
    
  
    load(sprintf("./mcmc_data_sim/results_sim_Matt_model%s_repeat%g.RData", model, i))

    bh_r_mcmc <- c(bh_data@sim$samples[[1]]$r, bh_data@sim$samples[[2]]$r, bh_data@sim$samples[[3]]$r, bh_data@sim$samples[[4]]$r)
    bh_k_mcmc <- c(bh_data@sim$samples[[1]]$k, bh_data@sim$samples[[2]]$k, bh_data@sim$samples[[3]]$k, bh_data@sim$samples[[4]]$k)
    bh_sigma_y_mcmc <- c(bh_data@sim$samples[[1]]$sigma_y, bh_data@sim$samples[[2]]$sigma_y, bh_data@sim$samples[[3]]$sigma_y, bh_data@sim$samples[[4]]$sigma_y)
    rk_r_mcmc <- c(rk_data@sim$samples[[1]]$r, rk_data@sim$samples[[2]]$r, rk_data@sim$samples[[3]]$r, rk_data@sim$samples[[4]]$r)
    rk_k_mcmc <- c(rk_data@sim$samples[[1]]$k, rk_data@sim$samples[[2]]$k, rk_data@sim$samples[[3]]$k, rk_data@sim$samples[[4]]$k)
    rk_sigma_y_mcmc <- c(rk_data@sim$samples[[1]]$sigma_y, rk_data@sim$samples[[2]]$sigma_y, rk_data@sim$samples[[3]]$sigma_y, rk_data@sim$samples[[4]]$sigma_y)
    hs_r_mcmc <- c(hs_data@sim$samples[[1]]$r, hs_data@sim$samples[[2]]$r, hs_data@sim$samples[[3]]$r, hs_data@sim$samples[[4]]$r)
    hs_k_mcmc <- c(hs_data@sim$samples[[1]]$k, hs_data@sim$samples[[2]]$k, hs_data@sim$samples[[3]]$k, hs_data@sim$samples[[4]]$k)
    hs_sigma_y_mcmc <- c(hs_data@sim$samples[[1]]$sigma_y, hs_data@sim$samples[[2]]$sigma_y, hs_data@sim$samples[[3]]$sigma_y, hs_data@sim$samples[[4]]$sigma_y)
    
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
      if (loglike_bh(params) < bh_ml_init){
        bh_r_init <- bh_r_mcmc[j]
        bh_k_init <- bh_k_mcmc[j]
        bh_sigma_y_init <- bh_sigma_y_mcmc[j]
        bh_ml_init <- loglike_bh(params)
      }

      params <- c(rk_r_mcmc[j], rk_k_mcmc[j], rk_sigma_y_mcmc[j])
      if (loglike_rk(params) < rk_ml_init){
        rk_ml_init <- loglike_rk(params)
        rk_r_init <- rk_r_mcmc[j]
        rk_k_init <- rk_k_mcmc[j]
        rk_sigma_y_init <- rk_sigma_y_mcmc[j]
      }

      params <- c(hs_r_mcmc[j], hs_k_mcmc[j], hs_sigma_y_mcmc[j])
      if (loglike_hs(params) < hs_ml_init){
        hs_ml_init <- loglike_hs(params)
        hs_r_init <- hs_r_mcmc[j]
        hs_k_init <- hs_k_mcmc[j]
        hs_sigma_y_init <- hs_sigma_y_mcmc[j]
      }
    }
      params_bh <- c(bh_r_init, bh_k_init, bh_sigma_y_init)
      params_rk <- c(rk_r_init, rk_k_init, rk_sigma_y_init)
      params_hs <- c(hs_r_init, hs_k_init, hs_sigma_y_init)

      if(model == "bh"){
        bh_opt <- hjkb(fn = loglike_bh, par = params_bh, lower = c(1.00000000001, 0.1*max(Ri), 0), upper = c(10, 10*max(Ri), Inf))
        bh_ml_bh[i] <- bh_opt$value
        bh_r_bh[i] <- bh_opt$par[1]
        bh_k_bh[i] <- bh_opt$par[2]
        bh_k_scaled_bh[i] <- bh_opt$par[2]/max(Ri)
        bh_sigma_y_bh[i] <- bh_opt$par[3]

        rk_opt <- hjkb(fn = loglike_rk, par = params_rk, lower = c(1.00000000001, 0.1*max(Ri), 0), upper = c(10, 10*max(Ri), Inf))
        rk_ml_bh[i] <- rk_opt$value
        rk_r_bh[i] <- rk_opt$par[1]
        rk_k_bh[i] <- rk_opt$par[2]
        rk_k_scaled_bh[i] <- rk_opt$par[2]/max(Ri)
        rk_sigma_y_bh[i] <- rk_opt$par[3]

        hs_opt <- hjkb(fn = loglike_hs, par = params_hs, lower = c(1.00000000001, 0.01*max(Ri), 0), upper = c(10, 10*max(Ri), Inf))
        hs_ml_bh[i] <- hs_opt$value
        hs_r_bh[i] <- hs_opt$par[1]
        hs_k_bh[i] <- hs_opt$par[2]
        hs_k_scaled_bh[i] <- hs_opt$par[2]/max(Ri)
        hs_sigma_y_bh[i] <- hs_opt$par[3]
      }

      if(model == "rk"){
        bh_opt <- hjkb(fn = loglike_bh, par = params_bh, lower = c(1.00000000001, 0.1*max(Ri), 0), upper = c(10, 10*max(Ri), Inf))
        bh_ml_rk[i] <- bh_opt$value
        bh_r_rk[i] <- bh_opt$par[1]
        bh_k_rk[i] <- bh_opt$par[2]
        bh_k_scaled_rk[i] <- bh_opt$par[2]/max(Ri)
        bh_sigma_y_rk[i] <- bh_opt$par[3]

        rk_opt <- hjkb(fn = loglike_rk, par = params_rk, lower = c(1.00000000001, 0.1*max(Ri), 0), upper = c(10, 10*max(Ri), Inf))
        rk_ml_rk[i] <- rk_opt$value
        rk_r_rk[i] <- rk_opt$par[1]
        rk_k_rk[i] <- rk_opt$par[2]
        rk_k_scaled_rk[i] <- rk_opt$par[2]/max(Ri)
        rk_sigma_y_rk[i] <- rk_opt$par[3]

        hs_opt <- hjkb(fn = loglike_hs, par = params_hs, lower = c(1.00000000001, 0.1*max(Ri), 0), upper = c(10, 10*max(Ri), Inf))
        hs_ml_rk[i] <- hs_opt$value
        hs_r_rk[i] <- hs_opt$par[1]
        hs_k_rk[i] <- hs_opt$par[2]
        hs_k_scaled_rk[i] <- hs_opt$par[2]/max(Ri)
        hs_sigma_y_rk[i] <- hs_opt$par[3]
      }

      if(model == "hs"){
        bh_opt <- hjkb(fn = loglike_bh, par = params_bh, lower = c(1.00000000001, 0.1*max(Ri), 0), upper = c(10, 10*max(Ri), Inf))
        bh_ml_hs[i] <- bh_opt$value
        bh_r_hs[i] <- bh_opt$par[1]
        bh_k_hs[i] <- bh_opt$par[2]
        bh_k_scaled_hs[i] <- bh_opt$par[2]/max(Ri)
        bh_sigma_y_hs[i] <- bh_opt$par[3]

        rk_opt <- hjkb(fn = loglike_rk, par = params_rk, lower = c(1.00000000001, 0.1*max(Ri), 0), upper = c(10, 10*max(Ri), Inf))
        rk_ml_hs[i] <- rk_opt$value
        rk_r_hs[i] <- rk_opt$par[1]
        rk_k_hs[i] <- rk_opt$par[2]
        rk_k_scaled_hs[i] <- rk_opt$par[2]/max(Ri)
        rk_sigma_y_hs[i] <- rk_opt$par[3]

        hs_opt <- hjkb(fn = loglike_hs, par = params_hs, lower = c(1.00000000001, 0.1*max(Ri), 0), upper = c(10, 10*max(Ri), Inf))
        hs_ml_hs[i] <- hs_opt$value
        hs_r_hs[i] <- hs_opt$par[1]
        hs_k_hs[i] <- hs_opt$par[2]
        hs_k_scaled_hs[i] <- hs_opt$par[2]/max(Ri)
        hs_sigma_y_hs[i] <- hs_opt$par[3]
      }

  }
}

# calculating AIC

bh_AIC_bh <- 2*bh_ml_bh
bh_AIC_rk <- 2*bh_ml_rk
bh_AIC_hs <- 2*bh_ml_hs
rk_AIC_bh <- 2*rk_ml_bh
rk_AIC_rk <- 2*rk_ml_rk
rk_AIC_hs <- 2*rk_ml_hs
hs_AIC_bh <- 2*hs_ml_bh
hs_AIC_rk <- 2*hs_ml_rk
hs_AIC_hs <- 2*hs_ml_hs

min_AIC_bh <- min_AIC_rk <- min_AIC_hs <-rep(NA, n_sims)

for(i in 1:n_sims){
  min_AIC_bh[i] <- min(c(bh_AIC_bh[i], rk_AIC_bh[i], hs_AIC_bh[i]))
  min_AIC_rk[i] <- min(c(bh_AIC_rk[i], rk_AIC_rk[i], hs_AIC_rk[i]))
  min_AIC_hs[i] <- min(c(bh_AIC_hs[i], rk_AIC_hs[i], hs_AIC_hs[i]))
}

bh_delta_AIC_bh <- bh_AIC_bh - min_AIC_bh
bh_delta_AIC_rk <- bh_AIC_rk - min_AIC_rk
bh_delta_AIC_hs <- bh_AIC_hs - min_AIC_hs
rk_delta_AIC_bh <- rk_AIC_bh - min_AIC_bh
rk_delta_AIC_rk <- rk_AIC_rk - min_AIC_rk
rk_delta_AIC_hs <- rk_AIC_hs - min_AIC_hs
hs_delta_AIC_bh <- hs_AIC_bh - min_AIC_bh
hs_delta_AIC_rk <- hs_AIC_rk - min_AIC_rk
hs_delta_AIC_hs <- hs_AIC_hs - min_AIC_hs

sum_neg_half_delta_AIC_bh <- exp((-1/2)*bh_delta_AIC_bh) + exp((-1/2)*rk_delta_AIC_bh) + exp((-1/2)*hs_delta_AIC_bh)
sum_neg_half_delta_AIC_rk <- exp((-1/2)*bh_delta_AIC_rk) + exp((-1/2)*rk_delta_AIC_rk) + exp((-1/2)*hs_delta_AIC_rk)
sum_neg_half_delta_AIC_hs <- exp((-1/2)*bh_delta_AIC_hs) + exp((-1/2)*rk_delta_AIC_hs) + exp((-1/2)*hs_delta_AIC_hs)

bh_waic_bh <- exp((-1/2)*bh_delta_AIC_bh)/sum_neg_half_delta_AIC_bh
bh_waic_rk <- exp((-1/2)*bh_delta_AIC_rk)/sum_neg_half_delta_AIC_rk
bh_waic_hs <- exp((-1/2)*bh_delta_AIC_hs)/sum_neg_half_delta_AIC_hs
rk_waic_bh <- exp((-1/2)*rk_delta_AIC_bh)/sum_neg_half_delta_AIC_bh
rk_waic_rk <- exp((-1/2)*rk_delta_AIC_rk)/sum_neg_half_delta_AIC_rk
rk_waic_hs <- exp((-1/2)*rk_delta_AIC_hs)/sum_neg_half_delta_AIC_hs
hs_waic_bh <- exp((-1/2)*hs_delta_AIC_bh)/sum_neg_half_delta_AIC_bh
hs_waic_rk <- exp((-1/2)*hs_delta_AIC_rk)/sum_neg_half_delta_AIC_rk
hs_waic_hs <- exp((-1/2)*hs_delta_AIC_hs)/sum_neg_half_delta_AIC_hs

saveRDS(bh_waic_bh, file = "./data_waic/bh_waic_sim_bh_true.rds")
saveRDS(rk_waic_bh, file = "./data_waic/rk_waic_sim_bh_true.rds")
saveRDS(hs_waic_bh, file = "./data_waic/hs_waic_sim_bh_true.rds")
saveRDS(bh_waic_rk, file = "./data_waic/bh_waic_sim_rk_true.rds")
saveRDS(rk_waic_rk, file = "./data_waic/rk_waic_sim_rk_true.rds")
saveRDS(hs_waic_rk, file = "./data_waic/hs_waic_sim_rk_true.rds")
saveRDS(bh_waic_hs, file = "./data_waic/bh_waic_sim_hs_true.rds")
saveRDS(rk_waic_hs, file = "./data_waic/rk_waic_sim_hs_true.rds")
saveRDS(hs_waic_hs, file = "./data_waic/hs_waic_sim_hs_true.rds")

saveRDS(bh_n_datapoints, file = "./data_n_datapoints/n_datapoints_sim.rds")

saveRDS(order(bh_n_datapoints, decreasing = FALSE), file = "./data_n_datapoints/n_datapoints_increasing_sim.rds")
