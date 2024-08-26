#reading in files

#model probailities
bh_model_prob <- readRDS("./data_model_prob/bh_model_prob.rds")
rk_model_prob <- readRDS("./data_model_prob/rk_model_prob.rds")
hs_model_prob <- readRDS("./data_model_prob/hs_model_prob.rds")

#weighted AIC's
bh_waic <- readRDS("./data_waic/bh_waic.rds")
rk_waic <- readRDS("./data_waic/rk_waic.rds")
hs_waic <- readRDS("./data_waic/hs_waic.rds")

n_datapoints <- readRDS("./data_n_datapoints/n_datapoints.rds")
n_datapoints_increasing <- readRDS("./data_n_datapoints/n_datapoints_increasing.rds")

#counts
bh_highest_model_prob <- rk_highest_model_prob <- hs_highest_model_prob <- c()
bh_highest_waic <- rk_highest_waic <- hs_highest_waic <- c()
bh_model_prob_bigger_than_50_percent <- rk_model_prob_bigger_than_50_percent <- hs_model_prob_bigger_than_50_percent <- c()

for(i in 1:length(n_datapoints)){
  #model prob
  if(bh_model_prob[i] == max(c(bh_model_prob[i], rk_model_prob[i], hs_model_prob[i]))){
    bh_highest_model_prob <- append(bh_highest_model_prob, i)
  }
  if(rk_model_prob[i] == max(c(bh_model_prob[i], rk_model_prob[i], hs_model_prob[i]))){
    rk_highest_model_prob <- append(rk_highest_model_prob, i)
  }
  if(hs_model_prob[i] == max(c(bh_model_prob[i], rk_model_prob[i], hs_model_prob[i]))){
    hs_highest_model_prob <- append(hs_highest_model_prob, i)
  }
  #waic
  if(bh_waic[i] == max(c(bh_waic[i], rk_waic[i], hs_waic[i]))){
    bh_highest_waic <- append(bh_highest_waic, i)
  }
  if(rk_waic[i] == max(c(bh_waic[i], rk_waic[i], hs_waic[i]))){
    rk_highest_waic <- append(rk_highest_waic, i)
  }
  if(hs_waic[i] == max(c(bh_waic[i], rk_waic[i], hs_waic[i]))){
    hs_highest_waic <- append(hs_highest_waic, i)
  }
}

n_bh_highest_model_prob <- length(bh_highest_model_prob)
n_rk_highest_model_prob <- length(rk_highest_model_prob)
n_hs_highest_model_prob <- length(hs_highest_model_prob)

n_bh_highest_waic <- length(bh_highest_waic)
n_rk_highest_waic <- length(rk_highest_waic)
n_hs_highest_waic <- length(hs_highest_waic)


#plots 
bh_col <- 'darkblue'
bh_col_tr <- 'lightblue'
rk_col <- 'red'
rk_col_tr <- 'pink'
hs_col <- 'darkgreen'
hs_col_tr <- 'lightgreen'


# waic ####
pdf(file = './plots_final/waic.pdf',
      width = 15,
      height = 7.5)
par(mar=c(5.1, 5.1, 1.1, 1.1), xpd=TRUE)
barplot(bh_waic[n_datapoints_increasing] + rk_waic[n_datapoints_increasing] + hs_waic[n_datapoints_increasing],
        col = rk_col,
        border = rk_col,
        cex.lab = 2,
        xaxt = 'n',
        yaxt = 'n',
        width = rep(0.138, 284))
title(ylab = "Akaike Weight", cex.lab = 3, line = 2)
title(xlab = "Number of datapoints", cex.lab = 3, line = 4)
axis(1, at=0:47, labels=sort(n_datapoints[seq(1, length(n_datapoints), 6)]), cex.axis = 2)
axis(2, at=seq(0, 1, 0.2), labels=c(0, 0.2, 0.4, 0.6, 0.8, 1), line = -2.25, cex.axis = 2)
barplot(bh_waic[n_datapoints_increasing] + hs_waic[n_datapoints_increasing],
        col = bh_col,
        border = bh_col,
        xaxt = 'n',
        yaxt = 'n',
        add = TRUE,
        width = rep(0.138, 284))
barplot(hs_waic[n_datapoints_increasing],
        col = hs_col,
        border = hs_col,
        xaxt = 'n',
        yaxt = 'n',
        add = TRUE,
        width = rep(0.138, 284))
dev.off()
# model prob ####
pdf(file = './plots_final/model_prob.pdf',
    width = 15,
    height = 7.5)
par(mar=c(5.1, 5.1, 1.1, 1.1), xpd=TRUE)
barplot(bh_model_prob[n_datapoints_increasing] + rk_model_prob[n_datapoints_increasing] + hs_model_prob[n_datapoints_increasing],
        col = rk_col,
        border = rk_col,
        cex.lab = 2,
        xaxt = 'n',
        yaxt = 'n',
        width = rep(0.138, 284))
title(ylab = "Model Probability", cex.lab = 3, line = 2)
title(xlab = "Number of datapoints", cex.lab = 3, line = 4)
axis(1, at=0:47, labels=sort(n_datapoints[seq(1, length(n_datapoints), 6)]), cex.axis = 2)
axis(2, at=seq(0, 1, 0.2), labels=c(0, 0.2, 0.4, 0.6, 0.8, 1), line = -2.25, cex.axis = 2)
barplot(bh_model_prob[n_datapoints_increasing] + hs_model_prob[n_datapoints_increasing],
        col = bh_col,
        border = bh_col,
        xaxt = 'n',
        yaxt = 'n',
        add = TRUE,
        width = rep(0.138, 284))
barplot(hs_model_prob[n_datapoints_increasing],
        col = hs_col,
        border = hs_col,
        xaxt = 'n',
        yaxt = 'n',
        add = TRUE,
        width = rep(0.138, 284))
dev.off()

# model prob and waic combined ####
pdf(file = './plots_final/model_prob_waic_combined.pdf',
    width = 15,
    height = 10)
par(mar=c(3, 5, 3, 0), xpd=TRUE, mfrow=c(2, 1), oma =c(3, 0, 0, 0))
barplot(bh_model_prob[n_datapoints_increasing] + rk_model_prob[n_datapoints_increasing] + hs_model_prob[n_datapoints_increasing],
        col = rk_col,
        border = rk_col,
        cex.lab = 2,
        xaxt = 'n',
        yaxt = 'n',
        width = rep(0.138, 284))
title(ylab = "Model Probability", cex.lab = 3, line = 2)
text(0.2, 1.1, "(a)", xpd = NA, cex = 2, adj = 0)

axis(2, at=seq(0, 1, 0.2), labels=c(0, 0.2, 0.4, 0.6, 0.8, 1), line = -2.25, cex.axis = 2)
barplot(bh_model_prob[n_datapoints_increasing] + hs_model_prob[n_datapoints_increasing],
        col = bh_col,
        border = bh_col,
        xaxt = 'n',
        yaxt = 'n',
        add = TRUE,
        width = rep(0.138, 284))
barplot(hs_model_prob[n_datapoints_increasing],
        col = hs_col,
        border = hs_col,
        xaxt = 'n',
        yaxt = 'n',
        add = TRUE,
        width = rep(0.138, 284))

par(mar=c(3, 5, 2, 0), xpd=TRUE)
barplot(bh_waic[n_datapoints_increasing] + rk_waic[n_datapoints_increasing] + hs_waic[n_datapoints_increasing],
        col = rk_col,
        border = rk_col,
        cex.lab = 2,
        xaxt = 'n',
        yaxt = 'n',
        width = rep(0.138, 284))
title(ylab = "Akaike Weight", cex.lab = 3, line = 2)
mtext("Number of Datapoints", side = 1, cex = 3, line = 4)
text(0.2, 1.1, "(b)", xpd = NA, cex = 2, adj = 0)
axis(1, at=0:47, labels=sort(n_datapoints[seq(1, length(n_datapoints), 6)]), cex.axis = 2)
axis(2, at=seq(0, 1, 0.2), labels=c(0, 0.2, 0.4, 0.6, 0.8, 1), line = -2.25, cex.axis = 2)
barplot(bh_waic[n_datapoints_increasing] + hs_waic[n_datapoints_increasing],
        col = bh_col,
        border = bh_col,
        xaxt = 'n',
        yaxt = 'n',
        add = TRUE,
        width = rep(0.138, 284))
barplot(hs_waic[n_datapoints_increasing],
        col = hs_col,
        border = hs_col,
        xaxt = 'n',
        yaxt = 'n',
        add = TRUE,
        width = rep(0.138, 284))
dev.off()

# scat model prob vs waic ####
pdf(file = './plots_final/scat_model_prob_vs_waic.pdf',
    width = 15,
    height = 15)
par(mar=c(7.1, 7.1, 1, 1), xpd=TRUE)
plot(c(0, 1), c(0, 1), col='white', pch=16, cex=0.000001,
     xlab = NA,
     ylab = NA,
     xaxt='n', 
     yaxt='n')
title(xlab = "Model Probability", cex.lab = 3, line = 4)
title(ylab = "Akaike Weight", cex.lab = 3, line = 4)
axis(1, at=seq(0, 1, 0.2), labels=c(0, 0.2, 0.4, 0.6, 0.8, 1))
axis(2, at=seq(0, 1, 0.2), labels=c(0, 0.2, 0.4, 0.6, 0.8, 1))
points(bh_model_prob, bh_waic, pch=15, cex=1.5, col=bh_col)
points(rk_model_prob, rk_waic, pch=16, cex=1.5, col=rk_col)
points(hs_model_prob, hs_waic, pch=17, cex=1.5, col=hs_col)
lines(c(-0.04, 1.04), c(-0.04, 1.04), lwd=2)
dev.off()

# fav_model_hist_waic ####
pdf(file = './plots_final/hist_fav_model_waic.pdf',
    width = 7.5,
    height = 7.5)
par(mar=c(6, 6, 1, 1), xpd=TRUE)
hist(c(bh_waic[bh_highest_waic], 
       rk_waic[rk_highest_waic], 
       hs_waic[hs_highest_waic]),
     xlim = c(0, 1),
     ylim = c(0, 50),
     cex.axis = 1.8,
     xlab = '',
     ylab = '',
     main = '',
     xaxt = 'n')
axis(1, at=seq(0, 1, 0.2), labels=seq(0, 1, 0.2), line = -1.15, cex.axis = 1.8)
mtext("Akaike Weight", side = 1, cex = 2.5, line = 4)
mtext("Frequency", side = 2, cex = 2.5, line = 4)
dev.off()
# fav_model_hist_model_prob_waic_combined ####
pdf(file = './plots_final/hist_fav_model_model_prob_waic_combined.pdf',
    width = 15,
    height = 7.5)
par(mar=c(0, 0, 0, 2), xpd=TRUE, mfrow=c(1, 2), oma =c(6, 6, 1, 1))
hist(c(bh_model_prob[bh_highest_model_prob], 
       rk_model_prob[rk_highest_model_prob], 
       hs_model_prob[hs_highest_model_prob]),
     xlim = c(0, 1),
     ylim = c(0, 50),
     cex.axis = 1.8,
     xlab = '',
     ylab = '',
     main = '',
     xaxt = 'n')
axis(1, at=seq(0, 1, 0.2), labels=seq(0, 1, 0.2), line = -1.15, cex.axis = 1.8)
mtext("Frequency", side = 2, cex = 2.5, line = 3.5)
mtext("Model Probability", side = 1, cex = 2.5, line = 4)
par(mar=c(0, 2, 0, 0))
hist(c(bh_waic[bh_highest_waic], 
       rk_waic[rk_highest_waic], 
       hs_waic[hs_highest_waic]),
     xlim = c(0, 1),
     ylim = c(0, 50),
     cex.axis = 1.5,
     xlab = '',
     ylab = '',
     main = '',
     xaxt = 'n',
     yaxt = 'n')
axis(1, at=seq(0, 1, 0.2), labels=seq(0, 1, 0.2), line = -1.15, cex.axis = 1.8)
mtext("Akaike Weight", side = 1, cex = 2.5, line = 4)
dev.off()
