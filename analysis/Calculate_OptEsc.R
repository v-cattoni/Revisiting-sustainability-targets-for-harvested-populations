# plot opt esc 

# reading in files ####
load("./data_files/RAMCore[asmt][v4.495].rdata") #RAM DATABASE
#this was visually identified manually recorded in the below csv file
cols <- read.csv("./data_files/non_det_cols.csv", header <- FALSE)[, 1]

hs_bf <- readRDS("./opt_esc_data/hs_bf.rds")
bh_bf <- readRDS("./opt_esc_data/bh_bf.rds")
rk_bf <- readRDS("./opt_esc_data/rk_bf.rds")

hs_sig_prob <- readRDS("./opt_esc_data/hs_sig_prob.rds")
bh_sig_prob <- readRDS("./opt_esc_data/bh_sig_prob.rds")
rk_sig_prob <- readRDS("./opt_esc_data/rk_sig_prob.rds")

# constants ####
years <- seq(1950, 2020, 1) # year of each row of data
n_rows <- length(years) # num rows to iterate through all data
phi <- 0.188 # constant in pella-tomlinson model
n_datasets <- length(cols) # number of dataets = number of cols 
n_iter = 10000 # number of iterations used in mcmc model
n_branch = 4 # number of branches used in mcmc model
lw <- 2
ca <- 1.2
cl <- 1.3

# functions ####
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

#Growth rates corresponding to equal optimal escapements as a proportion of fitted carrying capacity
bh_vs_hs_r <- function(r) {
  return((r - 1) / ((r ** 0.5) - 1))
}

# initialising ####
bio_mat <-
  matrix(NA, n_rows, n_datasets) # biomass (MT)
catch_mat <-
  matrix(NA, n_rows, n_datasets) # catch (MT)
Ri_max <-
  rep(NaN,n_datasets) # max Ri (MT)
bh_r <- rk_r  <- hs_r <- rep(NaN,n_datasets*n_iter*n_branch) #growth rate
bh_k <- rk_k  <- hs_k <- rep(NaN,n_datasets*n_iter*n_branch) #carrying capacity
bh_opt_esc <- rk_opt_esc  <- hs_opt_esc <- rep(NaN,n_datasets*n_iter*n_branch) #carrying capacity prop of fitted k
bh_opt_esc_scaled <- rk_opt_esc_scaled  <- hs_opt_esc_scaled <- rep(NaN,n_datasets*n_iter*n_branch) #carrying capacity prop of max pop
bh_opt_esc_bf <- rk_opt_esc_bf  <- hs_opt_esc_bf <- c() #carrying capacity prop of fitted k when model is bf
bh_opt_esc_scaled_bf <- rk_opt_esc_scaled_bf  <- hs_opt_esc_scaled_bf <- c() #carrying capacity prop of max pop when models is bf
bh_opt_esc_sig_prob <- rk_opt_esc_sig_prob  <- hs_opt_esc_sig_prob <- c() #carrying capacity prop of fitted k when model is bf
bh_opt_esc_scaled_sig_prob <- rk_opt_esc_scaled_sig_prob  <- hs_opt_esc_scaled_sig_prob <- c() #carrying capacity prop of max pop when models is bf
bh_pi <- hs_pi <- rk_pi  <- rep(NA, n_datasets) # best parameter index per model per dataset (40 000)
bh_opt_esc_pi <- hs_opt_esc_pi <- rk_opt_esc_pi  <- rep(NA, n_datasets) # opt esc for best parameter index
bh_opt_esc_scaled_pi <- hs_opt_esc_scaled_pi <- rk_opt_esc_scaled_pi  <- rep(NA, n_datasets) # as above as prop of max pop
bh_opt_esc_pi_bf <- hs_opt_esc_pi_bf <- rk_opt_esc_pi_bf <- pt_opt_esc_pi_bf <- c() # best fit
bh_opt_esc_scaled_pi_bf <- hs_opt_esc_scaled_pi_bf <- rk_opt_esc_scaled_pi_bf <-  c() # best fit of max pop
bh_opt_esc_pi_sig_prob <- hs_opt_esc_pi_sig_prob <- rk_opt_esc_pi_sig_prob <-  c() # model prob > 0.25
bh_opt_esc_scaled_pi_sig_prob <- hs_opt_esc_scaled_pi_sig_prob <- rk_opt_esc_scaled_pi_sig_prob <-  c() # prop of max pop when model prob > 0.25
hs_k_kmax_flag <- rep(NA, n_datasets*n_iter*n_branch)
hs_opt_esc_scaled_kmax <- hs_opt_esc_scaled_kmin <- rep(NA, n_datasets*n_branch*n_iter)
hs_opt_esc_scaled_pi_kmin <- hs_opt_esc_scaled_pi_kmax <- rep(NA, n_datasets)
hs_opt_esc_scaled_pi_bf_kmin <- hs_opt_esc_scaled_pi_bf_kmax <- c()

# sorting data ####
# creates a matrix of Biomass data and a matrix of Catch data with only
# overlapping entries.
for (i in cols) {
  for (j in 1:n_rows) {
    if (FALSE %in% is.na(Bio[, i][j]) &
        FALSE %in% is.na(MCatch[, i][j])) {
      bio_mat[j, match(i, cols)] <- c(Bio[, i][j])
      catch_mat[j, match(i, cols)] <-
        c(MCatch[, i][j])
    }
  }
}

# storing max biomass for each species ####
for (i in 1:length(cols)) {
  bios <- bio_mat[, i][!is.na(bio_mat[, i])]
  catchs <- catch_mat[, i][!is.na(catch_mat[, i])]
  X <-
    bios[1:(length(bios) - 1)] - catchs[1:(length(bios) - 1)]
  Y <-
    bios[2:length(bios)]
  Ri_max[i] <- 
    max(X, Y)
}


# storing growth rate and carrying capacity by dataset, model, branch, iteration
for (dataset_id in 1:n_datasets){
  print(dataset_id)
  load(sprintf("./mcmc_data/results_dataset%d.RData",dataset_id))
  hs_pi[dataset_id] <- which.max(c(hs_data@sim$samples[[1]]$lp__, hs_data@sim$samples[[2]]$lp__, hs_data@sim$samples[[3]]$lp__, hs_data@sim$samples[[4]]$lp__))
  bh_pi[dataset_id] <- which.max(c(bh_data@sim$samples[[1]]$lp__, bh_data@sim$samples[[2]]$lp__, bh_data@sim$samples[[3]]$lp__, bh_data@sim$samples[[4]]$lp__))
  rk_pi[dataset_id] <- which.max(c(rk_data@sim$samples[[1]]$lp__, rk_data@sim$samples[[2]]$lp__, rk_data@sim$samples[[3]]$lp__, rk_data@sim$samples[[4]]$lp__))
  for (branch in 1:n_branch){
    for (iter in 1:n_iter){
      bh_r[(dataset_id-1)*n_iter*n_branch + (branch - 1)*n_iter + iter] <- bh_data@sim$samples[[branch]]$r[iter]
      bh_k[(dataset_id-1)*n_iter*n_branch + (branch - 1)*n_iter + iter] <- bh_data@sim$samples[[branch]]$k[iter]/Ri_max[dataset_id]
      hs_r[(dataset_id-1)*n_iter*n_branch + (branch - 1)*n_iter + iter] <- hs_data@sim$samples[[branch]]$r[iter]
      hs_k[(dataset_id-1)*n_iter*n_branch + (branch - 1)*n_iter + iter] <- hs_data@sim$samples[[branch]]$k[iter]/Ri_max[dataset_id]
      rk_r[(dataset_id-1)*n_iter*n_branch + (branch - 1)*n_iter + iter] <- rk_data@sim$samples[[branch]]$r[iter]
      rk_k[(dataset_id-1)*n_iter*n_branch + (branch - 1)*n_iter + iter] <- rk_data@sim$samples[[branch]]$k[iter]/Ri_max[dataset_id]
      if (min(hs_data@sim$samples[[branch]]$k[iter]/Ri_max[dataset_id], 1) == 1){
        hs_k_kmax_flag[(dataset_id-1)*n_iter*n_branch + (branch - 1)*n_iter + iter] <- 1
      }
      else{
        hs_k_kmax_flag[(dataset_id-1)*n_iter*n_branch + (branch - 1)*n_iter + iter] <- 0
      }
    }
  }
}

print("opt esc fitted k")
# storing optimal escapement as a proportion of the fitted K ####
for (i in 1:(n_datasets*n_branch*n_iter)){
  print(i)
  bh_opt_esc[i] <- bh_opt_esc_f(bh_r[i], 1)
  hs_opt_esc[i] <- hs_opt_esc_f(hs_r[i], 1)
  p <- rk_r[i]
  k <- 1
  result <- uniroot(derivative_rk, lower = 0, upper = 10 * k)
  rk_opt_esc[i] <- result$root
  }

print("opt esc max pop")
# storing optimal escapement as a proportion of the max population ####
for (i in 1:(n_datasets*n_branch*n_iter)){
  print(i)
  bh_opt_esc_scaled[i] <- bh_opt_esc_f(bh_r[i], bh_k[i])
  hs_opt_esc_scaled[i] <- hs_opt_esc_f(hs_r[i], hs_k[i])
  p <- rk_r[i]
  k <- rk_k[i]
  result <- uniroot(derivative_rk, lower = 0, upper = 10 * k)
  rk_opt_esc_scaled[i] <- result$root
  if(hs_k[i] == 1){
    hs_opt_esc_scaled_kmax[i] <- hs_opt_esc_f(hs_r[i], hs_k[i])
  }
  else{
    hs_opt_esc_scaled_kmin[i] <- hs_opt_esc_f(hs_r[i], hs_k[i])
  }
  }

# storing optimal escapement for best parameters for every model dataset pair
for (i in 1:n_datasets){
  bh_opt_esc_pi[i] <- bh_opt_esc[(i-1)*n_iter*n_branch + bh_pi[i]]
  hs_opt_esc_pi[i] <- hs_opt_esc[(i-1)*n_iter*n_branch + hs_pi[i]]
  rk_opt_esc_pi[i] <- rk_opt_esc[(i-1)*n_iter*n_branch + rk_pi[i]]
  bh_opt_esc_scaled_pi[i] <- bh_opt_esc_scaled[(i-1)*n_iter*n_branch + bh_pi[i]]
  hs_opt_esc_scaled_pi[i] <- hs_opt_esc_scaled[(i-1)*n_iter*n_branch + hs_pi[i]]
  rk_opt_esc_scaled_pi[i] <- rk_opt_esc_scaled[(i-1)*n_iter*n_branch + rk_pi[i]]
  if(hs_k[(i-1)*n_iter*n_branch + hs_pi[i]] == 1){
    hs_opt_esc_scaled_pi_kmax[i] <- hs_opt_esc_scaled[(i-1)*n_iter*n_branch + hs_pi[i]]
  }
  else{
    hs_opt_esc_scaled_pi_kmin[i] <- hs_opt_esc_scaled[(i-1)*n_iter*n_branch + hs_pi[i]]
  }
  if(hs_k[(i-1)*n_iter*n_branch + hs_pi[i]] >= 1 & i %in% hs_bf){
    hs_opt_esc_scaled_pi_bf_kmax <- append(hs_opt_esc_scaled_pi_bf_kmax, hs_opt_esc_scaled[(i-1)*n_iter*n_branch + hs_pi[i]])
  }
  if(hs_k[(i-1)*n_iter*n_branch + hs_pi[i]] < 1 & i %in% hs_bf){
    hs_opt_esc_scaled_pi_bf_kmin <- append(hs_opt_esc_scaled_pi_bf_kmax, hs_opt_esc_scaled[(i-1)*n_iter*n_branch + hs_pi[i]])
  }
}

# storing optimal escapement for best fit models
for (i in hs_bf){
  hs_opt_esc_bf <- append(hs_opt_esc_bf, hs_opt_esc[((i - 1) * n_iter * n_branch + 1):(i * n_iter * n_branch)])
  hs_opt_esc_scaled_bf <- append(hs_opt_esc_scaled_bf, hs_opt_esc_scaled[((i - 1) * n_iter * n_branch + 1):(i * n_iter * n_branch)])
  hs_opt_esc_pi_bf <- append(hs_opt_esc_pi_bf, hs_opt_esc_pi[i])
  hs_opt_esc_scaled_pi_bf <- append(hs_opt_esc_scaled_pi_bf, hs_opt_esc_scaled_pi[i])
}

for (i in bh_bf){
  bh_opt_esc_bf <- append(bh_opt_esc_bf, bh_opt_esc[((i - 1) * n_iter * n_branch + 1):(i * n_iter * n_branch)])
  bh_opt_esc_scaled_bf <- append(bh_opt_esc_scaled_bf, bh_opt_esc_scaled[((i - 1) * n_iter * n_branch + 1):(i * n_iter * n_branch)])
  bh_opt_esc_pi_bf <- append(bh_opt_esc_pi_bf, bh_opt_esc_pi[i])
  bh_opt_esc_scaled_pi_bf <- append(bh_opt_esc_scaled_pi_bf, bh_opt_esc_scaled_pi[i])
}

for (i in rk_bf){
  rk_opt_esc_bf <- append(rk_opt_esc_bf, rk_opt_esc[((i - 1) * n_iter * n_branch + 1):(i * n_iter * n_branch)])
  rk_opt_esc_scaled_bf <- append(rk_opt_esc_scaled_bf, rk_opt_esc_scaled[((i - 1) * n_iter * n_branch + 1):(i * n_iter * n_branch)])
  rk_opt_esc_pi_bf <- append(rk_opt_esc_pi_bf, rk_opt_esc_pi[i])
  rk_opt_esc_scaled_pi_bf <- append(rk_opt_esc_scaled_pi_bf, rk_opt_esc_scaled_pi[i])
}

# storing optimal escapement for models with prob > 0.05
for (i in hs_sig_prob){
  hs_opt_esc_sig_prob <- append(hs_opt_esc_sig_prob, hs_opt_esc[((i - 1) * n_iter * n_branch + 1):(i * n_iter * n_branch)])
  hs_opt_esc_scaled_sig_prob <- append(hs_opt_esc_scaled_sig_prob, hs_opt_esc_scaled[((i - 1) * n_iter * n_branch + 1):(i * n_iter * n_branch)])
  hs_opt_esc_pi_sig_prob <- append(hs_opt_esc_pi_sig_prob, hs_opt_esc_pi[i])
  hs_opt_esc_scaled_pi_sig_prob <- append(hs_opt_esc_scaled_pi_sig_prob, hs_opt_esc_scaled_pi[i])
}

for (i in bh_sig_prob){
  bh_opt_esc_sig_prob <- append(bh_opt_esc_sig_prob, bh_opt_esc[((i - 1) * n_iter * n_branch + 1):(i * n_iter * n_branch)])
  bh_opt_esc_scaled_sig_prob <- append(bh_opt_esc_scaled_sig_prob, bh_opt_esc_scaled[((i - 1) * n_iter * n_branch + 1):(i * n_iter * n_branch)])
  bh_opt_esc_pi_sig_prob <- append(bh_opt_esc_pi_sig_prob, bh_opt_esc_pi[i])
  bh_opt_esc_scaled_pi_sig_prob <- append(bh_opt_esc_scaled_pi_sig_prob, bh_opt_esc_scaled_pi[i])
}

for (i in rk_sig_prob){
  rk_opt_esc_sig_prob <- append(rk_opt_esc_sig_prob, rk_opt_esc[((i - 1) * n_iter * n_branch + 1):(i * n_iter * n_branch)])
  rk_opt_esc_scaled_sig_prob <- append(rk_opt_esc_scaled_sig_prob, rk_opt_esc_scaled[((i - 1) * n_iter * n_branch + 1):(i * n_iter * n_branch)])
  rk_opt_esc_pi_sig_prob <- append(rk_opt_esc_pi_sig_prob, rk_opt_esc_pi[i])
  rk_opt_esc_scaled_pi_sig_prob <- append(rk_opt_esc_scaled_pi_sig_prob, rk_opt_esc_scaled_pi[i])
}




#saving opt esc data
saveRDS(bh_r, file = "./opt_esc_data/bh_r.rds")
saveRDS(hs_r, file = "./opt_esc_data/hs_r.rds")
saveRDS(rk_r, file = "./opt_esc_data/rk_r.rds")

saveRDS(bh_k, file = "./opt_esc_data/bh_k.rds")
saveRDS(hs_k, file = "./opt_esc_data/hs_k.rds")
saveRDS(rk_k, file = "./opt_esc_data/rk_k.rds")

saveRDS(bh_pi, file = "./opt_esc_data/bh_pi.rds")
saveRDS(hs_pi, file = "./opt_esc_data/hs_pi.rds")
saveRDS(rk_pi, file = "./opt_esc_data/rk_pi.rds")

saveRDS(bh_opt_esc, file = "./opt_esc_data/bh_opt_esc.rds")
saveRDS(hs_opt_esc, file = "./opt_esc_data/hs_opt_esc.rds")
saveRDS(rk_opt_esc, file = "./opt_esc_data/rk_opt_esc.rds")

saveRDS(bh_opt_esc_scaled, file = "./opt_esc_data/bh_opt_esc_scaled.rds")
saveRDS(hs_opt_esc_scaled, file = "./opt_esc_data/hs_opt_esc_scaled.rds")
saveRDS(rk_opt_esc_scaled, file = "./opt_esc_data/rk_opt_esc_scaled.rds")

saveRDS(bh_opt_esc_pi, file = "./opt_esc_data/bh_opt_esc_pi.rds")
saveRDS(hs_opt_esc_pi, file = "./opt_esc_data/hs_opt_esc_pi.rds")
saveRDS(rk_opt_esc_pi, file = "./opt_esc_data/rk_opt_esc_pi.rds")

saveRDS(bh_opt_esc_scaled_pi, file = "./opt_esc_data/bh_opt_esc_scaled_pi.rds")
saveRDS(hs_opt_esc_scaled_pi, file = "./opt_esc_data/hs_opt_esc_scaled_pi.rds")
saveRDS(rk_opt_esc_scaled_pi, file = "./opt_esc_data/rk_opt_esc_scaled_pi.rds")

saveRDS(bh_opt_esc_bf, file = "./opt_esc_data/bh_opt_esc_bf.rds")
saveRDS(hs_opt_esc_bf, file = "./opt_esc_data/hs_opt_esc_bf.rds")
saveRDS(rk_opt_esc_bf, file = "./opt_esc_data/rk_opt_esc_bf.rds")

saveRDS(bh_opt_esc_scaled_bf, file = "./opt_esc_data/bh_opt_esc_scaled_bf.rds")
saveRDS(hs_opt_esc_scaled_bf, file = "./opt_esc_data/hs_opt_esc_scaled_bf.rds")
saveRDS(rk_opt_esc_scaled_bf, file = "./opt_esc_data/rk_opt_esc_scaled_bf.rds")

saveRDS(bh_opt_esc_sig_prob, file = "./opt_esc_data/bh_opt_esc_sig_prob.rds")
saveRDS(hs_opt_esc_sig_prob, file = "./opt_esc_data/hs_opt_esc_sig_prob.rds")
saveRDS(rk_opt_esc_sig_prob, file = "./opt_esc_data/rk_opt_esc_sig_prob.rds")

saveRDS(bh_opt_esc_scaled_sig_prob, file = "./opt_esc_data/bh_opt_esc_scaled_sig_prob.rds")
saveRDS(hs_opt_esc_scaled_sig_prob, file = "./opt_esc_data/hs_opt_esc_scaled_sig_prob.rds")
saveRDS(rk_opt_esc_scaled_sig_prob, file = "./opt_esc_data/rk_opt_esc_scaled_sig_prob.rds")

saveRDS(bh_opt_esc_pi_bf, file = "./opt_esc_data/bh_opt_esc_pi_bf.rds")
saveRDS(hs_opt_esc_pi_bf, file = "./opt_esc_data/hs_opt_esc_pi_bf.rds")
saveRDS(rk_opt_esc_pi_bf, file = "./opt_esc_data/rk_opt_esc_pi_bf.rds")

saveRDS(bh_opt_esc_scaled_pi_bf, file = "./opt_esc_data/bh_opt_esc_scaled_pi_bf.rds")
saveRDS(hs_opt_esc_scaled_pi_bf, file = "./opt_esc_data/hs_opt_esc_scaled_pi_bf.rds")
saveRDS(rk_opt_esc_scaled_pi_bf, file = "./opt_esc_data/rk_opt_esc_scaled_pi_bf.rds")

saveRDS(bh_opt_esc_pi_sig_prob, file = "./opt_esc_data/bh_opt_esc_pi_sig_prob.rds")
saveRDS(hs_opt_esc_pi_sig_prob, file = "./opt_esc_data/hs_opt_esc_pi_sig_prob.rds")
saveRDS(rk_opt_esc_pi_sig_prob, file = "./opt_esc_data/rk_opt_esc_pi_sig_prob.rds")

saveRDS(bh_opt_esc_scaled_pi_sig_prob, file = "./opt_esc_data/bh_opt_esc_scaled_pi_sig_prob.rds")
saveRDS(hs_opt_esc_scaled_pi_sig_prob, file = "./opt_esc_data/hs_opt_esc_scaled_pi_sig_prob.rds")
saveRDS(rk_opt_esc_scaled_pi_sig_prob, file = "./opt_esc_data/rk_opt_esc_scaled_pi_sig_prob.rds")
