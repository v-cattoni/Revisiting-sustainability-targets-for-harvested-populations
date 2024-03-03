#### Getting the dataset index ####
# The following two lines get the input from the csv file
args<-commandArgs(TRUE)
my_repeat <- eval( parse(text=args[1]) )

setwd("/home/southl/BayesFish_NoLoglike")

#### READING IN FILES ####
load("./data_files/RAMCore[asmt][v4.495].rdata") #ram database
cols <- read.csv("./data_files/non_det_cols.csv", header <- FALSE)[, 1] #non-deterministic data sets

#### CONSTANTS ####
num_rows <- 
  length(Bio[,1]) #number of years
num_cols <- 
  length(cols) #number of non-deterministic data sets

#### INITIALISING ####
bio_mat <-
  matrix(NA, num_rows, num_cols) #biomass (year, species)
catch_mat <-
  matrix(NA, num_rows, num_cols) #catch(MT) (year, species)

#### SORTING DATA ####
#Creates a matrix of Biomass data and a matrix of Catch data with only overlapping entries.
for (i in cols) {
  for (j in 1:num_rows) {
    if (FALSE %in% is.na(Bio[, i][j]) &
        FALSE %in% is.na(MCatch[, i][j])) {
      bio_mat[j, match(i, cols)] <- c(Bio[, i][j])
      catch_mat[j, match(i, cols)] <-
        c(MCatch[, i][j])
    }
  }
}

#### MAIN ####
set.seed(my_repeat+100)
dataset_id <- sample(1:284,size=1)
load(file = sprintf("./mcmc_data/results_dataset%d.RData",dataset_id))

print(dataset_id)
bios_real <- bio_mat[, dataset_id][!is.na(bio_mat[, dataset_id])] 
catchs_real <- catch_mat[, dataset_id][!is.na(catch_mat[, dataset_id])]
props_real <- catchs_real/bios_real
n <- length(bios_real)-1


##################
# MODEL: Beverton-Holt
set.seed(my_repeat)

# Parameter value approximately drawn from posterior
r <- bh_data@sim$samples[[1]]$r[10^4]
k <- bh_data@sim$samples[[1]]$k[10^4]
sigma_y <- bh_data@sim$samples[[1]]$sigma_y[10^4]

# Simulating data from the model (Matt's approach)
bios <- rep(NaN,n+1)
catchs <- Ei <- rep(NaN,n)
bios[1] <- bios_real[1]
props <- sample(props_real,n+1,replace=TRUE)
for (i in 1:n){
  catchs[i] <- bios[i]*props[i]
  Ei[i] <- bios[i] - catchs[i] #Escapement = biomass - catch in year i
  bios[i+1] <- rlnorm(1,log(Ei[i]) + log(r) - log(1 + (r - 1) * Ei[i] / k), sigma_y)
}
Ri <- bios[2:length(bios)] #Recruitment = biomass in year t+1

# Saving the data
sim_bh <- list(n = n,
             Ri = Ri,
             Ei = Ei,
             max_r = max(Ri)
)


##################
# MODEL: Ricker
set.seed(my_repeat)

# Parameter value approximately drawn from posterior
r <- rk_data@sim$samples[[1]]$r[10^4]
k <- rk_data@sim$samples[[1]]$k[10^4]
sigma_y <- rk_data@sim$samples[[1]]$sigma_y[10^4]

# Simulating data from the model (Matt's approach)
bios <- rep(NaN,n+1)
catchs <- Ei <- rep(NaN,n)
bios[1] <- bios_real[1]
props <- sample(props_real,n+1,replace=TRUE)
for (i in 1:n){
  catchs[i] <- bios[i]*props[i]
  Ei[i] <- bios[i] - catchs[i] #Escapement = biomass - catch in year i
  bios[i+1] <- rlnorm(1,log(Ei[i]) + (1 - Ei[i] / k) * log(r), sigma_y)
}
Ri <- bios[2:length(bios)] #Recruitment = biomass in year t+1

# Saving the data
sim_rk <- list(n = n,
                Ri = Ri,
                Ei = Ei,
                max_r = max(Ri)
)

##################
# MODEL: Pella-Tomlinson
set.seed(my_repeat)

# Parameter value approximately drawn from posterior
r <- pt_data@sim$samples[[1]]$r[10^4]
k <- pt_data@sim$samples[[1]]$k[10^4]
sigma_y <- pt_data@sim$samples[[1]]$sigma_y[10^4]

# Simulating data from the model (Matt's approach)
bios <- rep(NaN,n+1)
catchs <- Ei <- rep(NaN,n)
bios[1] <- bios_real[1]
props <- sample(props_real,n+1,replace=TRUE)
for (i in 1:n){
  catchs[i] <- bios[i]*props[i]
  Ei[i] <- bios[i] - catchs[i] #Escapement = biomass - catch in year i
  bios[i+1] <- rlnorm(1,log(Ei[i] + (r - 1) * Ei[i] * (1 - (Ei[i] / k) ^ 0.188)), sigma_y)
}
Ri <- bios[2:length(bios)] #Recruitment = biomass in year t+1

# Saving the data
sim_pt <- list(n = n,
               Ri = Ri,
               Ei = Ei,
               max_r = max(Ri)
)

##################
# MODEL: Hockey-Stick
set.seed(my_repeat)

# Parameter value approximately drawn from posterior
r <- hs_data@sim$samples[[1]]$r[10^4]
k <- hs_data@sim$samples[[1]]$k[10^4]
sigma_y <- hs_data@sim$samples[[1]]$sigma_y[10^4]

# Simulating data from the model (Matt's approach)
bios <- rep(NaN,n+1)
catchs <- Ei <- rep(NaN,n)
bios[1] <- bios_real[1]
props <- sample(props_real,n+1,replace=TRUE)
for (i in 1:n){
  catchs[i] <- bios[i]*props[i]
  Ei[i] <- bios[i] - catchs[i] #Escapement = biomass - catch in year i
  bios[i+1] <- rlnorm(1,log(min(r * Ei[i], k)), sigma_y)
}
Ri <- bios[2:length(bios)] #Recruitment = biomass in year t+1

# Saving the data
sim_hs <- list(n = n,
               Ri = Ri,
               Ei = Ei,
               max_r = max(Ri)
)

save(sim_bh, sim_rk, sim_pt, sim_hs, dataset_id, file = sprintf("./data_files_sim/sim_data_Matt_repeat%d.RData",my_repeat))

