#reading in files

#model probailities
#bh true
bh_model_prob_sim_bh_true <- readRDS("./data_model_prob/bh_model_prob_sim_bh_true.rds")
rk_model_prob_sim_bh_true <- readRDS("./data_model_prob/rk_model_prob_sim_bh_true.rds")
hs_model_prob_sim_bh_true <- readRDS("./data_model_prob/hs_model_prob_sim_bh_true.rds")

#rk true
bh_model_prob_sim_rk_true <- readRDS("./data_model_prob/bh_model_prob_sim_rk_true.rds")
rk_model_prob_sim_rk_true <- readRDS("./data_model_prob/rk_model_prob_sim_rk_true.rds")
hs_model_prob_sim_rk_true <- readRDS("./data_model_prob/hs_model_prob_sim_rk_true.rds")

#hs true
bh_model_prob_sim_hs_true <- readRDS("./data_model_prob/bh_model_prob_sim_hs_true.rds")
rk_model_prob_sim_hs_true <- readRDS("./data_model_prob/rk_model_prob_sim_hs_true.rds")
hs_model_prob_sim_hs_true <- readRDS("./data_model_prob/hs_model_prob_sim_hs_true.rds")

#weighted AIC's
#bh true
bh_waic_sim_bh_true <- readRDS("./data_waic/bh_waic_sim_bh_true.rds")
rk_waic_sim_bh_true <- readRDS("./data_waic/rk_waic_sim_bh_true.rds")
hs_waic_sim_bh_true <- readRDS("./data_waic/hs_waic_sim_bh_true.rds")

#rk true
bh_waic_sim_rk_true <- readRDS("./data_waic/bh_waic_sim_rk_true.rds")
rk_waic_sim_rk_true <- readRDS("./data_waic/rk_waic_sim_rk_true.rds")
hs_waic_sim_rk_true <- readRDS("./data_waic/hs_waic_sim_rk_true.rds")

#hs true
bh_waic_sim_hs_true <- readRDS("./data_waic/bh_waic_sim_hs_true.rds")
rk_waic_sim_hs_true <- readRDS("./data_waic/rk_waic_sim_hs_true.rds")
hs_waic_sim_hs_true <- readRDS("./data_waic/hs_waic_sim_hs_true.rds")

#number of datapoints
n_datapoints_increasing_sim <- readRDS("./data_n_datapoints/n_datapoints_increasing_sim.rds")
n_datapoints_sim <- readRDS("./data_n_datapoints/n_datapoints_sim.rds")

#plots 

#plots 
bh_col <- 'darkblue'
bh_col_tr <- 'lightblue'
rk_col <- 'red'
rk_col_tr <- 'pink'
hs_col <- 'darkgreen'
hs_col_tr <- 'lightgreen'

# waic bh true ####

pdf(file = './plots_final/model_prob_waic_combined_sim.pdf',
    width = 15,
    height = 10)
par(mar=c(1, 1, 5, 1), xpd=TRUE, mfrow=c(2, 3), oma =c(5, 5, 0,0))
# model prob bh true
par(mar=c(1, 1, 5, 1), xpd=TRUE)
barplot(bh_model_prob_sim_bh_true[n_datapoints_increasing_sim] + rk_model_prob_sim_bh_true[n_datapoints_increasing_sim] + hs_model_prob_sim_bh_true[n_datapoints_increasing_sim],
        col = rk_col,
        border = rk_col,
        cex.lab = 2,
        xaxt = 'n',
        yaxt = 'n',
        width = rep(0.138, 284))
text(0.2, 1.1, "(a) Beverton-Holt True", xpd = NA, cex = 3, adj = 0)
mtext("Model Probability", side = 2, cex = 2.5, line = 3.5)
axis(2, at=seq(0, 1, 0.25), labels=c(0, 0.25, 0.5, 0.75, 1), line = -0.7, cex.axis = 2)
barplot(bh_model_prob_sim_bh_true[n_datapoints_increasing_sim] + hs_model_prob_sim_bh_true[n_datapoints_increasing_sim],
        col = bh_col,
        border = bh_col,
        xaxt = 'n',
        yaxt = 'n',
        add = TRUE,
        width = rep(0.138, 284))
barplot(hs_model_prob_sim_bh_true[n_datapoints_increasing_sim],
        col = hs_col,
        border = hs_col,
        xaxt = 'n',
        yaxt = 'n',
        add = TRUE,
        width = rep(0.138, 284))


# model prob rk true
par(mar=c(1, 1, 5, 1), xpd=TRUE)
barplot(bh_model_prob_sim_rk_true[n_datapoints_increasing_sim] + rk_model_prob_sim_rk_true[n_datapoints_increasing_sim] + hs_model_prob_sim_rk_true[n_datapoints_increasing_sim],
        col = rk_col,
        border = rk_col,
        cex.lab = 2,
        xaxt = 'n',
        yaxt = 'n',
        width = rep(0.138, 284))
text(0.2, 1.1, "(b) Ricker True", xpd = NA, cex = 3, adj = 0)
axis(2, at=seq(0, 1, 0.25), labels=c(0, 0.25, 0.5, 0.75, 1), line = -0.7, cex.axis = 2)
barplot(bh_model_prob_sim_rk_true[n_datapoints_increasing_sim] + hs_model_prob_sim_rk_true[n_datapoints_increasing_sim],
        col = bh_col,
        border = bh_col,
        xaxt = 'n',
        yaxt = 'n',
        add = TRUE,
        width = rep(0.138, 284))
barplot(hs_model_prob_sim_rk_true[n_datapoints_increasing_sim],
        col = hs_col,
        border = hs_col,
        xaxt = 'n',
        yaxt = 'n',
        add = TRUE,
        width = rep(0.138, 284))

# model prob hs true
par(mar=c(1, 1, 5, 1), xpd=TRUE)
barplot(bh_model_prob_sim_hs_true[n_datapoints_increasing_sim] + rk_model_prob_sim_hs_true[n_datapoints_increasing_sim] + hs_model_prob_sim_hs_true[n_datapoints_increasing_sim],
        col = rk_col,
        border = rk_col,
        cex.lab = 2,
        yaxt = 'n',
        width = rep(0.138, 284))
text(0.2, 1.1, "(c) Hockey-Stick True", xpd = NA, cex = 3, adj = 0)
axis(2, at=seq(0, 1, 0.25), labels=c(0, 0.25, 0.5, 0.75, 1), line = -0.7, cex.axis = 2)
barplot(bh_model_prob_sim_hs_true[n_datapoints_increasing_sim] + hs_model_prob_sim_hs_true[n_datapoints_increasing_sim],
        col = bh_col,
        border = bh_col,
        xaxt = 'n',
        yaxt = 'n',
        add = TRUE,
        width = rep(0.138, 284))
barplot(hs_model_prob_sim_hs_true[n_datapoints_increasing_sim],
        col = hs_col,
        border = hs_col,
        xaxt = 'n',
        yaxt = 'n',
        add = TRUE,
        width = rep(0.138, 284))

barplot(bh_waic_sim_bh_true[n_datapoints_increasing_sim] + rk_waic_sim_bh_true[n_datapoints_increasing_sim] + hs_waic_sim_bh_true[n_datapoints_increasing_sim],
        col = rk_col,
        border = rk_col,
        cex.lab = 2,
        xaxt = 'n',
        yaxt = 'n',
        width = rep(0.138, 284))
text(0.2, 1.1, "(d)", xpd = NA, cex = 3, adj = 0)
mtext("Akaike Weight", side = 2, cex = 2.5, line = 3.5)
axis(2, at=seq(0, 1, 0.25), labels=c(0, 0.25, 0.5, 0.75, 1), line = -0.7, cex.axis = 2)
axis(1, at=seq(0.1, 3.25, length.out = 20), labels=n_datapoints_sim[n_datapoints_increasing_sim], cex.axis = 2, line = 0.2)
barplot(bh_waic_sim_bh_true[n_datapoints_increasing_sim] + hs_waic_sim_bh_true[n_datapoints_increasing_sim],
        col = bh_col,
        border = bh_col,
        xaxt = 'n',
        yaxt = 'n',
        add = TRUE,
        width = rep(0.138, 284))
barplot(hs_waic_sim_bh_true[n_datapoints_increasing_sim],
        col = hs_col,
        border = hs_col,
        xaxt = 'n',
        yaxt = 'n',
        add = TRUE,
        width = rep(0.138, 284))


# waic rk true
par(mar=c(1, 1, 5, 1), xpd=TRUE)
barplot(bh_waic_sim_rk_true[n_datapoints_increasing_sim] + rk_waic_sim_rk_true[n_datapoints_increasing_sim] + hs_waic_sim_rk_true[n_datapoints_increasing_sim],
        col = rk_col,
        border = rk_col,
        cex.lab = 2,
        xaxt = 'n',
        yaxt = 'n',
        width = rep(0.138, 284))
text(0.2, 1.1, "(e)", xpd = NA, cex = 3, adj = 0)
axis(1, at=seq(0.1, 3.25, length.out = 20), labels=n_datapoints_sim[n_datapoints_increasing_sim], cex.axis = 2, line = 0.2)
axis(2, at=seq(0, 1, 0.25), labels=c(0, 0.25, 0.5, 0.75, 1), line = -0.7, cex.axis = 2)
barplot(bh_waic_sim_rk_true[n_datapoints_increasing_sim] + hs_waic_sim_rk_true[n_datapoints_increasing_sim],
        col = bh_col,
        border = bh_col,
        xaxt = 'n',
        yaxt = 'n',
        add = TRUE,
        width = rep(0.138, 284))
barplot(hs_waic_sim_rk_true[n_datapoints_increasing_sim],
        col = hs_col,
        border = hs_col,
        xaxt = 'n',
        yaxt = 'n',
        add = TRUE,
        width = rep(0.138, 284))

# waic hs true
par(mar=c(1, 1, 5, 1), xpd=TRUE)
barplot(bh_waic_sim_hs_true[n_datapoints_increasing_sim] + rk_waic_sim_hs_true[n_datapoints_increasing_sim] + hs_waic_sim_hs_true[n_datapoints_increasing_sim],
        col = rk_col,
        border = rk_col,
        cex.lab = 2,
        yaxt = 'n',
        width = rep(0.138, 284))
text(0.2, 1.1, "(f)", xpd = NA, cex = 3, adj = 0)
axis(1, at=seq(0.1, 3.25, length.out = 20), labels=n_datapoints_sim[n_datapoints_increasing_sim], cex.axis = 2, line = 0.2)
axis(2, at=seq(0, 1, 0.25), labels=c(0, 0.25, 0.5, 0.75, 1), line = -0.7, cex.axis = 2)
mtext("Number of datapoints", side = 1, cex = 2.5, line = 4.5, adj = -8.5)
barplot(bh_waic_sim_hs_true[n_datapoints_increasing_sim] + hs_waic_sim_hs_true[n_datapoints_increasing_sim],
        col = bh_col,
        border = bh_col,
        xaxt = 'n',
        yaxt = 'n',
        add = TRUE,
        width = rep(0.138, 284))
barplot(hs_waic_sim_hs_true[n_datapoints_increasing_sim],
        col = hs_col,
        border = hs_col,
        xaxt = 'n',
        yaxt = 'n',
        add = TRUE,
        width = rep(0.138, 284))
dev.off()


