#### Getting the dataset index ####
# The following two lines get the input from the csv file
args<-commandArgs(TRUE)
args <- eval( parse(text=args[1]) )
my_repeat <- args%%20 # Gets the index out of 20
model_index <- (args-my_repeat)/20+1
if (my_repeat==0){
  my_repeat <- 20
  model_index <- model_index - 1
}

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
load(sprintf("./data_files_sim/sim_data_Matt_repeat%d.RData",my_repeat)) # simulated dataset
warmups <- 5*10^4

total_iterations <- 10^5

max_treedepth <-  30

n_chains <-  4

adapt_delta <- 0.999

if (model_index==1){
  true_model <- "bh"
  data <- sim_bh
} else if (model_index==2){
  true_model <- "rk"
  data <- sim_rk
} else if (model_index==3){
  true_model <- "pt"
  data <- sim_pt
} else if (model_index==4){
  true_model <- "hs"
  data <- sim_hs
}

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

save(bh_data, rk_data, pt_data, hs_data, vanilla_data, bh_logml, rk_logml, pt_logml, hs_logml, vanilla_logml, file = sprintf("./mcmc_data_sim/results_sim_Matt_model%s_repeat%d.RData",true_model,my_repeat))


