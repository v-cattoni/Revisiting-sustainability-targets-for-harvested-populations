#### Reading in files ####
load("./data_files/RAMCore[asmt][v4.495].rdata") 
cols <- read.csv("./data_files/non_det_cols.csv", header <- FALSE)[, 1]
load("./mcmc_data/results_modelchoice.RData")

#initialising
years <- seq(1950, 2020, 1) #year of each row of data
num_rows <- length(years) #num rows to iterate through all data
num_cols <- length(cols) #num cols to iterate through all data
bio_mat <- matrix(NA, num_rows, num_cols) #biomass(MT)
catch_mat <- matrix(NA, num_rows, num_cols) #catch(MT)
n_years <- rep(NA, num_cols)

#### Sorting data ####
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

for (i in 1:length(cols)) {
  bios <- bio_mat[, i][!is.na(bio_mat[, i])]
  catchs <- catch_mat[, i][!is.na(catch_mat[, i])]
  X <-
    bios[1:(length(bios) - 1)] - catchs[1:(length(bios) - 1)]
  Y <- bios[2:length(bios)]
  n_years[i] <- length(X)
}


# colours ####
bh_col <- 'darkblue'
bh_col_tr <- '#4281C388'
rk_col <- 'red'
rk_col_tr <- '#E4231388'
hs_col <- 'darkgreen'
hs_col_tr <- '#00843788'


#ordering x axis
n_years_index_acsending = order(n_years)
hs_model_prob_data_ascending <- modelchoice$p[modelchoice$method == "hs"][n_years_index_acsending]
bh_model_prob_data_ascending <- modelchoice$p[modelchoice$method == "bh"][n_years_index_acsending]
rk_model_prob_data_ascending <- modelchoice$p[modelchoice$method == "rk"][n_years_index_acsending]



pdf(file = './plots_final/3_way_model_prob_data_increasing_1.pdf',
    width = 15,
    height = 7.5)
par(mar=c(5.1, 3.1, 2.1, 2.1), xpd=TRUE)
barplot(rk_model_prob_data_ascending + bh_model_prob_data_ascending + hs_model_prob_data_ascending, 
        col = bh_col, 
        border = bh_col,
        cex.lab = 2,
        xaxt = 'n',
        yaxt = 'n',
        width = rep(0.138, 284)
)
title(ylab = "Model Probability", cex.lab = 2, line = 1)
title(xlab = "Number of datapoints", cex.lab = 2, line = 3)
axis(1, at=0:47, labels=sort(n_years[seq(1, length(n_years), 6)]))
axis(2, at=seq(0, 1, 0.2), labels=c(0, 0.2, 0.4, 0.6, 0.8, 1), line = -2.25)

barplot(hs_model_prob_data_ascending + rk_model_prob_data_ascending,
        col = rk_col, 
        border = rk_col, 
        xaxt = 'n',
        yaxt = 'n',
        add = TRUE,
        width = rep(0.138, 284)
)
barplot(hs_model_prob_data_ascending,
        col = hs_col,
        border = hs_col, 
        xaxt = 'n',
        yaxt = 'n',
        add = TRUE,
        width = rep(0.138, 284)
)
dev.off()

pdf(file = './plots_final/3_way_model_prob_data_increasing_2.pdf',
    width = 15,
    height = 7.5)
par(mar=c(5.1, 3.1, 5.1, 2.1), xpd=TRUE)
barplot(rk_model_prob_data_ascending + bh_model_prob_data_ascending + hs_model_prob_data_ascending, 
        col = hs_col, 
        border = hs_col,
        cex.lab = 2,
        xaxt = 'n',
        yaxt = 'n',
        width = rep(0.138, 284)
)
title(main = "Model Probability Comparison by Dataset", cex.main = 2)
title(ylab = "Model Probability", cex.lab = 2, line = 1)
title(xlab = "Number of datapoints", cex.lab = 2, line = 3)
axis(1, at=0:47, labels=sort(n_years[seq(1, length(n_years), 6)]))
axis(2, at=seq(0, 1, 0.2), labels=c(0, 0.2, 0.4, 0.6, 0.8, 1), line = -2.25)

barplot(rk_model_prob_data_ascending + bh_model_prob_data_ascending,
        col = rk_col, 
        border = rk_col, 
        xaxt = 'n',
        yaxt = 'n',
        add = TRUE,
        width = rep(0.138, 284)
)
barplot(bh_model_prob_data_ascending,
        col = bh_col,
        border = bh_col, 
        xaxt = 'n',
        yaxt = 'n',
        add = TRUE,
        width = rep(0.138, 284)
)
dev.off()

pdf(file = './plots_final/3_way_model_prob_data_increasing_3.pdf',
    width = 15,
    height = 7.5)
par(mar=c(5.1, 3.1, 5.1, 2.1), xpd=TRUE)
barplot(rk_model_prob_data_ascending + bh_model_prob_data_ascending + hs_model_prob_data_ascending, 
        col = rk_col, 
        border = rk_col,
        cex.lab = 2,
        xaxt = 'n',
        yaxt = 'n',
        width = rep(0.138, 284)
)
title(main = "Model Probability Comparison by Dataset", cex.main = 2)
title(ylab = "Model Probability", cex.lab = 2, line = 1)
title(xlab = "Number of datapoints", cex.lab = 2, line = 3)
axis(1, at=0:47, labels=sort(n_years[seq(1, length(n_years), 6)]))
axis(2, at=seq(0, 1, 0.2), labels=c(0, 0.2, 0.4, 0.6, 0.8, 1), line = -2.25)

barplot(hs_model_prob_data_ascending + bh_model_prob_data_ascending,
        col = bh_col, 
        border = bh_col, 
        xaxt = 'n',
        yaxt = 'n',
        add = TRUE,
        width = rep(0.138, 284)
)
barplot(hs_model_prob_data_ascending,
        col = hs_col,
        border = hs_col, 
        xaxt = 'n',
        yaxt = 'n',
        add = TRUE,
        width = rep(0.138, 284)
)
dev.off()


# count number of times each model has highest probability ####

hs_bf <- bh_bf <- rk_bf <- c() # best fit model per dataset
hs_bf_par <- bh_bf_par <- rk_bf_par <- c() # best fit parameter values per model dataset
hs_sig_prob <- bh_sig_prob <- rk_sig_prob <- c() # best fit parameter values per model dataset
for(i in 1:length(modelchoice$p[modelchoice$method == "hs"])){

  # best fit models
  if(modelchoice$p[modelchoice$method == "hs"][i] == max(modelchoice$p[modelchoice$dataset == i][1:3])){
    hs_bf <- append(hs_bf, i)
  }
  if(modelchoice$p[modelchoice$method == "bh"][i] == max(modelchoice$p[modelchoice$dataset == i][1:3])){
    bh_bf <- append(bh_bf, i)
  }
  if(modelchoice$p[modelchoice$method == "rk"][i] == max(modelchoice$p[modelchoice$dataset == i][1:3])){
    rk_bf <- append(rk_bf, i)
  }

  # models with prob > 0.25
  if(modelchoice$p[modelchoice$method == "hs"][i] > 1/3){
    hs_sig_prob <- append(hs_sig_prob, i)
  }
  if(modelchoice$p[modelchoice$method == "bh"][i] > 1/3){
    bh_sig_prob <- append(bh_sig_prob, i)
  }
  if(modelchoice$p[modelchoice$method == "rk"][i] > 1/3){
    rk_sig_prob <- append(rk_sig_prob, i)
  }
}

saveRDS(hs_bf, file = "./opt_esc_data/hs_bf.rds")
saveRDS(bh_bf, file = "./opt_esc_data/bh_bf.rds")
saveRDS(rk_bf, file = "./opt_esc_data/rk_bf.rds")

saveRDS(hs_sig_prob, file = "./opt_esc_data/hs_sig_prob.rds")
saveRDS(bh_sig_prob, file = "./opt_esc_data/bh_sig_prob.rds")
saveRDS(rk_sig_prob, file = "./opt_esc_data/rk_sig_prob.rds")