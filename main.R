#### Getting the dataset index ####
# The following two lines get the input from the csv file
args<-commandArgs(TRUE)
dataset_id <- eval( parse(text=args[1]) )

#### INSTALLING PACKAGES ####
# install.packages(c("bridgesampling", "coda", "rstan", "lbfgsb3c", "loo"))

# Loading packages
mylib.loc <- "/home/southl/BayesFish_NoLoglike/pkg"
library(bridgesampling)
library(coda)
library(rstan)
library(lbfgsb3c, lib.loc=mylib.loc)
library(loo)

set.seed(123)
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
i <- dataset_id

print(i)
bios <- bio_mat[, i][!is.na(bio_mat[, i])] 
catchs <- catch_mat[, i][!is.na(catch_mat[, i])]

#data
Ei <- bios[1:(length(bios) - 1)] - catchs[1:(length(bios) - 1)] #Escapement = biomass - catch in year i
Ri <- bios[2:length(bios)] #Recruitment = biomass in year t+1

warmups <- 5*10^4

total_iterations <- 10^5

max_treedepth <-  30

n_chains <-  4

adapt_delta <- 0.999

data <- list(n = length(Ei),
             Ri = Ri,
             Ei = Ei,
             max_r = max(Ri)
)


#### BEVERTON-HOLT ####
set.seed(123)
bh_fit <- stan(
  file = "./models/bh_model.stan",
  data = data,
  chains = n_chains,
  warmup = warmups,
  thin = 10,
  iter = total_iterations,
  control = list(max_treedepth = max_treedepth,
                 adapt_delta = adapt_delta)
)
bh_data <- bh_fit
bh_logml <- bridge_sampler(bh_fit)$logml


#### RICKER ####
set.seed(123)
rk_fit <- stan(
  file = "./models/rk_model.stan",
  data = data,
  chains = n_chains,
  warmup = warmups,
  thin = 10,
  iter = total_iterations,
  control = list(max_treedepth = max_treedepth,
                 adapt_delta = adapt_delta)
)
rk_data <- rk_fit
rk_logml <- bridge_sampler(rk_fit)$logml

#### PELLA-TOMLINSON ####
set.seed(123)
pt_fit <- stan(
  file = "./models/pt_model.stan",
  data = data,
  chains = n_chains,
  warmup = warmups,
  thin = 10,
  iter = total_iterations,
  control = list(max_treedepth = max_treedepth,
                 adapt_delta = adapt_delta)
)
pt_data <- pt_fit
pt_logml <- bridge_sampler(pt_fit)$logml

#### HOCKEY-STICK ####
set.seed(123)
hs_fit <- stan(
  file = "./models/hs_model.stan",
  data = data,
  chains = n_chains,
  warmup = warmups,
  thin = 10,
  iter = total_iterations,
  control = list(max_treedepth = max_treedepth,
                 adapt_delta = adapt_delta)
)
hs_data <- hs_fit
hs_logml <- bridge_sampler(hs_fit)$logml

#### VANILLA MODEL ####
set.seed(123)
vanilla_fit <- stan(
  file = "./models/vanilla_model.stan",
  data = data,
  chains = n_chains,
  warmup = warmups,
  thin = 10,
  iter = total_iterations,
  control = list(max_treedepth = max_treedepth,
                 adapt_delta = adapt_delta)
)
vanilla_data <- vanilla_fit
vanilla_logml <- bridge_sampler(vanilla_fit)$logml

save(bh_data, rk_data, pt_data, hs_data, vanilla_data, bh_logml, rk_logml, pt_logml, hs_logml, vanilla_logml, file = sprintf("./mcmc_data/results_dataset%d.RData",dataset_id))


