# Load and combine the datasets

n_datasets <- 20
for(x in c('Matt', 'David')){
  for(model in c('hs', 'bh', 'rk', 'pt')){
    bh <- rk <- pt <- hs <- vanilla <- rep(NaN,n_datasets)
    for (dataset_id in 1:n_datasets){
      load(sprintf("./mcmc_data_sim/results_sim_%s_model%s_repeat%d.RData",x, model, dataset_id))
      bh[dataset_id] <- bh_logml
      rk[dataset_id] <- rk_logml
      pt[dataset_id] <- pt_logml
      hs[dataset_id] <- hs_logml
      vanilla[dataset_id] <- vanilla_logml
      print(dataset_id)
    }
    
    # Stable calculation of a vector dividing by its sum, but
    # where everything (including the response) is on the log
    # scale for stability
    logsumexp <- function(x){
      myMax <- max(x)
      x <- x - myMax
      return (log(sum(exp(x))) + myMax)
    }
    
    # Combining all log marginal likelihood estimates
    logZ <- cbind(bh, rk, hs, pt, vanilla)
    
    # Stable calculations of model probability. 
    
    # Including vanilla model
    p_van <- t(apply(logZ,1,function(x) exp(x - logsumexp(x))))
    
    # Excluding vanilla model
    p <- t(apply(logZ[,1:3],1,function(x) exp(x - logsumexp(x))))
    p <- cbind(p,rep(NaN,n_datasets), rep(NaN,n_datasets))
    
    # Putting together the results in a data frame
    modelchoice <- data.frame(method = c(rep("bh",n_datasets),rep("rk",n_datasets),rep("hs",n_datasets),rep("pt",n_datasets),rep("vanilla",n_datasets)),
                              dataset = rep(1:n_datasets,5),
                              logZ = c(bh, rk, hs, pt, vanilla),
                              p = c(p),
                              p_van = c(p_van))
    
    # Checking that it's combined correctly
    # test <- subset(modelchoice,dataset==11)
    # sum(test$p_van)
    
    save(modelchoice,file = sprintf("./mcmc_data_sim/results_modelchoice_sim_%s_model%s.RData", x, model))
  }
}

for(model in c('hs', 'bh', 'rk')){
  load(sprintf("./mcmc_data_sim/results_modelchoice_sim_Matt_model%s.RData", model))
  if(model == "bh"){
    saveRDS(modelchoice$p[modelchoice$method == "bh"], file = "./data_model_prob/bh_model_prob_sim_bh_true.rds")
    saveRDS(modelchoice$p[modelchoice$method == "rk"], file = "./data_model_prob/rk_model_prob_sim_bh_true.rds")
    saveRDS(modelchoice$p[modelchoice$method == "hs"], file = "./data_model_prob/hs_model_prob_sim_bh_true.rds")
  }
  
  if(model == "rk"){
    saveRDS(modelchoice$p[modelchoice$method == "bh"], file = "./data_model_prob/bh_model_prob_sim_rk_true.rds")
    saveRDS(modelchoice$p[modelchoice$method == "rk"], file = "./data_model_prob/rk_model_prob_sim_rk_true.rds")
    saveRDS(modelchoice$p[modelchoice$method == "hs"], file = "./data_model_prob/hs_model_prob_sim_rk_true.rds")
  }
  
  if(model == "hs"){
    saveRDS(modelchoice$p[modelchoice$method == "bh"], file = "./data_model_prob/bh_model_prob_sim_hs_true.rds")
    saveRDS(modelchoice$p[modelchoice$method == "rk"], file = "./data_model_prob/rk_model_prob_sim_hs_true.rds")
    saveRDS(modelchoice$p[modelchoice$method == "hs"], file = "./data_model_prob/hs_model_prob_sim_hs_true.rds")
  }
}