# colours ####
bh_col <- 'darkblue'
bh_col_tr <- '#4281C388'
rk_col <- 'red'
rk_col_tr <- '#E4231388'
hs_col <- 'darkgreen'
hs_col_tr <- '#00843788'
pt_col <- '#DB7093'
pt_col_tr <- '#DB709388'

n_datasets <- 20
# counting fits
hs_hs_bf_Matt <- bh_hs_bf_Matt <- rk_hs_bf_Matt <- c() # best fit model per dataset when data is generated from hs
hs_bh_bf_Matt <- bh_bh_bf_Matt <- rk_bh_bf_Matt <- c() # best fit model per dataset when data is generated from bh
hs_rk_bf_Matt <- bh_rk_bf_Matt <- rk_rk_bf_Matt <- c() # best fit model per dataset when data is generated from rk

hs_hs_bf_David <- bh_hs_bf_David <- rk_hs_bf_David <- c() # best fit model per dataset when data is generated from hs
hs_bh_bf_David <- bh_bh_bf_David <- rk_bh_bf_David <- c() # best fit model per dataset when data is generated from bh
hs_rk_bf_David <- bh_rk_bf_David <- rk_rk_bf_David <- c() # best fit model per dataset when data is generated from rk

for(x in c('Matt', 'David')){
  for(model in c('hs', 'bh', 'rk', 'pt')){
    bh <- rk <- pt <- hs <- vanilla <- rep(NaN,n_datasets)
    load(sprintf("./mcmc_data_sim/results_modelchoice_sim_%s_model%s.RData", x, model))
    if(model == "bh"){
      if(x == "David"){
      for(i in 1:length(modelchoice$p[modelchoice$method == "bh"])){
        # best fit models
        if(modelchoice$p[modelchoice$method == "hs"][i] == max(modelchoice$p[modelchoice$dataset == i][1:3])){
          hs_bh_bf_David <- append(hs_bh_bf_David, i)
        }
        if(modelchoice$p[modelchoice$method == "bh"][i] == max(modelchoice$p[modelchoice$dataset == i][1:3])){
          bh_bh_bf_David <- append(bh_bh_bf_David, i)
        }
        if(modelchoice$p[modelchoice$method == "rk"][i] == max(modelchoice$p[modelchoice$dataset == i][1:3])){
          rk_bh_bf_David <- append(rk_bh_bf_David, i)
        }
      }}
      if(x == "Matt"){
        for(i in 1:length(modelchoice$p[modelchoice$method == "bh"])){
          # best fit models
          if(modelchoice$p[modelchoice$method == "hs"][i] == max(modelchoice$p[modelchoice$dataset == i][1:3])){
            hs_bh_bf_Matt <- append(hs_bh_bf_Matt, i)
          }
          if(modelchoice$p[modelchoice$method == "bh"][i] == max(modelchoice$p[modelchoice$dataset == i][1:3])){
            bh_bh_bf_Matt <- append(bh_bh_bf_Matt, i)
          }
          if(modelchoice$p[modelchoice$method == "rk"][i] == max(modelchoice$p[modelchoice$dataset == i][1:3])){
            rk_bh_bf_Matt <- append(rk_bh_bf_Matt, i)
          }
        }}
  
      bh_model_prob_descending_index <- order(modelchoice$p[modelchoice$method == "bh"], decreasing = TRUE)
      hs_bh_model_prob_descending <- modelchoice$p[modelchoice$method == "hs"][bh_model_prob_descending_index]
      bh_bh_model_prob_descending <- modelchoice$p[modelchoice$method == "bh"][bh_model_prob_descending_index]
      rk_bh_model_prob_descending <- modelchoice$p[modelchoice$method == "rk"][bh_model_prob_descending_index]
      
      pdf(file = sprintf('./plots_final/3_way_model_prob_bh_decreasing_sim_%s.pdf', x),
          width = 15,
          height = 7.5)
      par(mar=c(4.1, 5.1, 5.1, 2.1), xpd=TRUE)
      barplot(bh_bh_model_prob_descending + rk_bh_model_prob_descending + hs_bh_model_prob_descending,
              col = hs_col,
              border = hs_col,
              cex.lab = 2,
              axis.lty = 0
      )
      title(main = "Model Probability Comparison by Dataset", cex.main = 2)
      title(ylab = "Model Probability", cex.lab = 2)
      title(xlab = "Dataset", cex.lab = 2, line = 1.5)

      barplot(bh_bh_model_prob_descending + rk_bh_model_prob_descending,
              col = rk_col,
              border = rk_col,
              add = TRUE
      )
      barplot(bh_bh_model_prob_descending,
              col = bh_col,
              border = bh_col,
              add = TRUE
      )
      dev.off()
    }

    if(model == "hs"){
      if(x == "David"){
        for(i in 1:length(modelchoice$p[modelchoice$method == "hs"])){
          # best fit models
          if(modelchoice$p[modelchoice$method == "hs"][i] == max(modelchoice$p[modelchoice$dataset == i][1:3])){
            hs_hs_bf_David <- append(hs_hs_bf_David, i)
          }
          if(modelchoice$p[modelchoice$method == "bh"][i] == max(modelchoice$p[modelchoice$dataset == i][1:3])){
            bh_hs_bf_David <- append(bh_hs_bf_David, i)
          }
          if(modelchoice$p[modelchoice$method == "rk"][i] == max(modelchoice$p[modelchoice$dataset == i][1:3])){
            rk_hs_bf_David <- append(rk_hs_bf_David, i)
          }
        }}
      if(x == "Matt"){
        for(i in 1:length(modelchoice$p[modelchoice$method == "hs"])){
          # best fit models
          if(modelchoice$p[modelchoice$method == "hs"][i] == max(modelchoice$p[modelchoice$dataset == i][1:3])){
            hs_hs_bf_Matt <- append(hs_hs_bf_Matt, i)
          }
          if(modelchoice$p[modelchoice$method == "bh"][i] == max(modelchoice$p[modelchoice$dataset == i][1:3])){
            bh_hs_bf_Matt <- append(bh_hs_bf_Matt, i)
          }
          if(modelchoice$p[modelchoice$method == "rk"][i] == max(modelchoice$p[modelchoice$dataset == i][1:3])){
            rk_hs_bf_Matt <- append(rk_hs_bf_Matt, i)
          }
        }}
      hs_model_prob_descending_index <- order(modelchoice$p[modelchoice$method == "hs"], decreasing = TRUE)
      hs_hs_model_prob_descending <- modelchoice$p[modelchoice$method == "hs"][hs_model_prob_descending_index]
      bh_hs_model_prob_descending <- modelchoice$p[modelchoice$method == "bh"][hs_model_prob_descending_index]
      rk_hs_model_prob_descending <- modelchoice$p[modelchoice$method == "rk"][hs_model_prob_descending_index]

      pdf(file = sprintf('./plots_final/3_way_model_prob_hs_decreasing_sim_%s.pdf', x),
          width = 15,
          height = 7.5)
      par(mar=c(4.1, 5.1, 5.1, 2.1), xpd=TRUE)
      barplot(hs_hs_model_prob_descending + bh_hs_model_prob_descending + rk_hs_model_prob_descending,
              col = rk_col,
              border = rk_col,
              cex.lab = 2,
              axis.lty = 0
              )
      title(main = "Model Probability Comparison by Dataset", cex.main = 2)
      title(ylab = "Model Probability", cex.lab = 2)
      title(xlab = "Dataset", cex.lab = 2, line = 1.5)

      barplot(hs_hs_model_prob_descending + bh_hs_model_prob_descending,
              col = bh_col,
              border = bh_col,
              add = TRUE
              )
      barplot(hs_hs_model_prob_descending,
              col = hs_col,
              border = hs_col,
              add = TRUE
              )
      dev.off()
    }

    if(model == "rk"){
      if(x == "David"){
        for(i in 1:length(modelchoice$p[modelchoice$method == "rk"])){
          # best fit models
          if(modelchoice$p[modelchoice$method == "hs"][i] == max(modelchoice$p[modelchoice$dataset == i][1:3])){
            hs_rk_bf_David <- append(hs_rk_bf_David, i)
          }
          if(modelchoice$p[modelchoice$method == "bh"][i] == max(modelchoice$p[modelchoice$dataset == i][1:3])){
            bh_rk_bf_David <- append(bh_rk_bf_David, i)
          }
          if(modelchoice$p[modelchoice$method == "rk"][i] == max(modelchoice$p[modelchoice$dataset == i][1:3])){
            rk_rk_bf_David <- append(rk_rk_bf_David, i)
          }
        }}
      
      if(x == "Matt"){
        for(i in 1:length(modelchoice$p[modelchoice$method == "rk"])){
          # best fit models
          if(modelchoice$p[modelchoice$method == "hs"][i] == max(modelchoice$p[modelchoice$dataset == i][1:3])){
            hs_rk_bf_Matt <- append(hs_rk_bf_Matt, i)
          }
          if(modelchoice$p[modelchoice$method == "bh"][i] == max(modelchoice$p[modelchoice$dataset == i][1:3])){
            bh_rk_bf_Matt <- append(bh_rk_bf_Matt, i)
          }
          if(modelchoice$p[modelchoice$method == "rk"][i] == max(modelchoice$p[modelchoice$dataset == i][1:3])){
            rk_rk_bf_Matt <- append(rk_rk_bf_Matt, i)
          }
        }}
      rk_model_prob_descending_index <- order(modelchoice$p[modelchoice$method == "rk"], decreasing = TRUE)
      hs_rk_model_prob_descending <- modelchoice$p[modelchoice$method == "hs"][rk_model_prob_descending_index]
      bh_rk_model_prob_descending <- modelchoice$p[modelchoice$method == "bh"][rk_model_prob_descending_index]
      rk_rk_model_prob_descending <- modelchoice$p[modelchoice$method == "rk"][rk_model_prob_descending_index]

      pdf(file = sprintf('./plots_final/3_way_model_prob_rk_decreasing_sim_%s.pdf', x),
          width = 15,
          height = 7.5)
      par(mar=c(4.1, 5.1, 5.1, 3.1), xpd=TRUE)
      barplot(rk_rk_model_prob_descending + bh_rk_model_prob_descending + hs_rk_model_prob_descending,
              col = hs_col,
              border = hs_col,
              cex.lab = 2,
              axis.lty = 0
      )
      title(main = "Model Probability Comparison by Dataset", cex.main = 2)
      title(ylab = "Model Probability", cex.lab = 2)
      title(xlab = "Dataset", cex.lab = 2, line = 1.5)

      barplot(rk_rk_model_prob_descending + bh_rk_model_prob_descending,
              col = bh_col,
              border = bh_col,
              add = TRUE
      )
      barplot(rk_rk_model_prob_descending,
              col = rk_col,
              border = rk_col,
              add = TRUE
      )
      dev.off()
    }
  }
}
