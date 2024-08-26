# plotting optimal escapement bayes and likelihood

# reading in data ####

# optimal escapement as a proportion of fitted k
bh_opt_esc_prop_fitted_k_bayes <-      readRDS("./data_opt_esc/bh_opt_esc_prop_fitted_k_bayes.rds")
hs_opt_esc_prop_fitted_k_bayes <-      readRDS("./data_opt_esc/hs_opt_esc_prop_fitted_k_bayes.rds")
rk_opt_esc_prop_fitted_k_bayes <-      readRDS("./data_opt_esc/rk_opt_esc_prop_fitted_k_bayes.rds")
bh_opt_esc_prop_fitted_k_likelihood <- readRDS("./data_opt_esc/bh_opt_esc_prop_fitted_k_likelihood.rds")
rk_opt_esc_prop_fitted_k_likelihood <- readRDS("./data_opt_esc/rk_opt_esc_prop_fitted_k_likelihood.rds")
hs_opt_esc_prop_fitted_k_likelihood <- readRDS("./data_opt_esc/hs_opt_esc_prop_fitted_k_likelihood.rds")

# optimal escapement as a proportion of max pop
bh_opt_esc_prop_max_pop_bayes <-      readRDS("./data_opt_esc/bh_opt_esc_prop_max_pop_bayes.rds")
hs_opt_esc_prop_max_pop_bayes <-      readRDS("./data_opt_esc/hs_opt_esc_prop_max_pop_bayes.rds")
rk_opt_esc_prop_max_pop_bayes <-      readRDS("./data_opt_esc/rk_opt_esc_prop_max_pop_bayes.rds")
bh_opt_esc_prop_max_pop_likelihood <- readRDS("./data_opt_esc/bh_opt_esc_prop_max_pop_likelihood.rds")
rk_opt_esc_prop_max_pop_likelihood <- readRDS("./data_opt_esc/rk_opt_esc_prop_max_pop_likelihood.rds")
hs_opt_esc_prop_max_pop_likelihood <- readRDS("./data_opt_esc/hs_opt_esc_prop_max_pop_likelihood.rds")

# constants ####
bh_col <- 'darkblue'
ri_col <- 'red'
hs_col <- 'darkgreen'

# opt esc hist combined likelihood ####
pdf(file = './plots_final/opt_esc_hist_combined_likelihood.pdf',
    width = 7,
    height = 3.5)
par(
  mfrow = c(3, 2),
  mar = c(0, 0, 0.7, 2),
  oma = c(3.5, 3.7, 0, 2)
)
hist(
  bh_opt_esc_prop_fitted_k_likelihood,
  yaxs = "i",
  yaxt = 'n',
  xlim = c(0, 1),
  ylim = c(0, 100),
  breaks = seq(from = -100, to = 100, by = 0.025),
  main = '',
  xlab = '',
  ylab = '',
  cex.lab = 1.5,
  xaxt = 'n',
  yaxt = 'n',
  col = bh_col,
  cex.axis = 2
)
text(0.25, 90, 'a) Proportion of fitted k', cex = 1.2)
axis(
  side = 2,
  at = c(0),
  labels = c(0),
  cex.axis = 1.25,
  mgp = c(2, 0.5, 0)
)
axis(1,
     pos = 0,
     cex.axis = 1.25,
     mgp = c(2, 0.5, 0))
lines(c(-1, 1), c(0, 0))
hist(
  bh_opt_esc_prop_max_pop_likelihood,
  yaxs = "i",
  yaxt = 'n',
  xlim = c(0, 2),
  ylim = c(0, 75),
  breaks = seq(from = -100, to = 100, by = 0.05),
  main = '',
  xlab = '',
  ylab = '',
  cex.lab = 1.5,
  xaxt = 'n',
  col = bh_col,
  cex.axis = 2
)
axis(
  side = 4,
  at = c(0),
  labels = c(0),
  cex.axis = 1.25,
  mgp = c(2, 0.5, 0)
)
axis(1,
     pos = 0,
     cex.axis = 1.25,
     mgp = c(2, 0.5, 0))
lines(c(-1, 1), c(0, 0))
text(0.7, 68, 'b) Proportion of max biomass', cex = 1.2)
hist(
  rk_opt_esc_prop_fitted_k_likelihood,
  yaxs = "i",
  yaxt = 'n',
  xlim = c(0, 1),
  ylim = c(0, 150),
  breaks = seq(from = -100, to = 100, by = 0.025),
  main = '',
  xlab = '',
  ylab = '',
  cex.lab = 1.5,
  xaxt = 'n',
  col = ri_col,
  cex.axis = 2
)
axis(
  side = 2,
  at = c(0),
  labels = c(0),
  cex.axis = 1.25,
  mgp = c(2, 0.5, 0)
)
axis(1,
     pos = 0,
     cex.axis = 1.25,
     mgp = c(2, 0.5, 0))
hist(
  rk_opt_esc_prop_max_pop_likelihood,
  yaxs = "i",
  yaxt = 'n',
  xlim = c(0, 2),
  ylim = c(0, 75),
  breaks = seq(from = -100, to = 100, by = 0.05),
  main = '',
  xlab = '',
  ylab = '',
  cex.lab = 1.5,
  xaxt = 'n',
  col = ri_col,
  cex.axis = 2
)
axis(
  side = 4,
  at = c(0),
  labels = c(0),
  cex.axis = 1.25,
  mgp = c(2, 0.5, 0)
)
axis(1,
     pos = 0,
     cex.axis = 1.25,
     mgp = c(2, 0.5, 0))
hist(
  hs_opt_esc_prop_fitted_k_likelihood,
  yaxs = "i",
  yaxt = 'n',
  xlim = c(0, 1),
  ylim = c(0, 35),
  breaks = seq(from = -100, to = 100, by = 0.025),
  main = '',
  xlab = '',
  ylab = '',
  xaxt = 'n',
  col = hs_col,
  cex.axis = 2
)
axis(
  side = 2,
  at = c(0),
  labels = c(0),
  cex.axis = 1.25,
  mgp = c(2, 0.5, 0)
)
axis(1,
     pos = 0,
     cex.axis = 1.25,
     mgp = c(2, 0.5, 0))
hist(
  hs_opt_esc_prop_max_pop_likelihood,
  yaxs = "i",
  yaxt = 'n',
  xlim = c(0, 2),
  ylim = c(0, 50),
  breaks = seq(from = -100, to = 100, by = 0.05),
  main = '',
  xlab = '',
  ylab = '',
  xaxt = 'n',
  col = hs_col,
  cex.axis = 2
)
axis(
  side = 4,
  at = c(0),
  labels = c(0),
  cex.axis = 1.25,
  mgp = c(2, 0.5, 0)
)
axis(1,
     pos = 0,
     cex.axis = 1.25,
     mgp = c(2, 0.5, 0))
mtext(
  side = 1,
  "Optimal escapement",
  outer = TRUE,
  padj = 1.9,
  cex = 1.3
)
mtext(
  side = 2,
  "Frequency",
  outer = TRUE,
  padj = -1.9,
  cex = 1.3
)
mtext(
  side = 1,
  "Hockey-Stick",
  outer = TRUE,
  padj = -4,
  adj = -1.6,
  at = 0.5,
  col = hs_col
)
mtext(
  side = 1,
  "Ricker",
  outer = TRUE,
  padj = -12,
  adj = -4,
  at = 0.5,
  col = ri_col
)
mtext(
  side = 1,
  "Beverton-Holt",
  outer = TRUE,
  padj = -21,
  adj = -1.5,
  at = 0.5,
  col = bh_col
)

dev.off()

# opt esc hist combined likelihood capped ####
pdf(file = './plots_final/opt_esc_hist_combined_likelihood_capped.pdf',
    width = 7,
    height = 3.5)
par(
  mfrow = c(3, 2),
  mar = c(0, 0, 0.7, 2),
  oma = c(3.5, 3.7, 0, 2)
)
hist(
  bh_opt_esc_prop_fitted_k_likelihood,
  yaxs = "i",
  yaxt = 'n',
  xlim = c(0, 1),
  ylim = c(0, 100),
  breaks = seq(from = -100, to = 100, by = 0.025),
  main = '',
  xlab = '',
  ylab = '',
  cex.lab = 1.5,
  xaxt = 'n',
  yaxt = 'n',
  col = bh_col,
  cex.axis = 2
)
text(0.25, 90, 'a) Proportion of fitted k', cex = 1.2)
axis(
  side = 2,
  at = c(0),
  labels = c(0),
  cex.axis = 1.25,
  mgp = c(2, 0.5, 0)
)
axis(1,
     pos = 0,
     cex.axis = 1.25,
     mgp = c(2, 0.5, 0))
lines(c(-1, 1), c(0, 0))
hist(
  bh_opt_esc_prop_max_pop_likelihood,
  yaxs = "i",
  yaxt = 'n',
  xlim = c(0, 2),
  ylim = c(0, 75),
  breaks = seq(from = -100, to = 100, by = 0.05),
  main = '',
  xlab = '',
  ylab = '',
  cex.lab = 1.5,
  xaxt = 'n',
  col = bh_col,
  cex.axis = 2
)
axis(
  side = 4,
  at = c(0),
  labels = c(0),
  cex.axis = 1.25,
  mgp = c(2, 0.5, 0)
)
axis(1,
     pos = 0,
     cex.axis = 1.25,
     mgp = c(2, 0.5, 0))
lines(c(-1, 1), c(0, 0))
text(0.7, 68, 'b) Proportion of max biomass', cex = 1.2)
hist(
  rk_opt_esc_prop_fitted_k_likelihood,
  yaxs = "i",
  yaxt = 'n',
  xlim = c(0, 1),
  ylim = c(0, 150),
  breaks = seq(from = -100, to = 100, by = 0.025),
  main = '',
  xlab = '',
  ylab = '',
  cex.lab = 1.5,
  xaxt = 'n',
  col = ri_col,
  cex.axis = 2
)
axis(
  side = 2,
  at = c(0),
  labels = c(0),
  cex.axis = 1.25,
  mgp = c(2, 0.5, 0)
)
axis(1,
     pos = 0,
     cex.axis = 1.25,
     mgp = c(2, 0.5, 0))
hist(
  rk_opt_esc_prop_max_pop_likelihood,
  yaxs = "i",
  yaxt = 'n',
  xlim = c(0, 2),
  ylim = c(0, 75),
  breaks = seq(from = -100, to = 100, by = 0.05),
  main = '',
  xlab = '',
  ylab = '',
  cex.lab = 1.5,
  xaxt = 'n',
  col = ri_col,
  cex.axis = 2
)
axis(
  side = 4,
  at = c(0),
  labels = c(0),
  cex.axis = 1.25,
  mgp = c(2, 0.5, 0)
)
axis(1,
     pos = 0,
     cex.axis = 1.25,
     mgp = c(2, 0.5, 0))
hist(
  hs_opt_esc_prop_fitted_k_likelihood,
  yaxs = "i",
  yaxt = 'n',
  xlim = c(0, 1),
  ylim = c(0, 35),
  breaks = seq(from = -100, to = 100, by = 0.025),
  main = '',
  xlab = '',
  ylab = '',
  xaxt = 'n',
  col = hs_col,
  cex.axis = 2
)
axis(
  side = 2,
  at = c(0),
  labels = c(0),
  cex.axis = 1.25,
  mgp = c(2, 0.5, 0)
)
axis(1,
     pos = 0,
     cex.axis = 1.25,
     mgp = c(2, 0.5, 0))
hist(
  hs_opt_esc_prop_max_pop_likelihood[hs_opt_esc_prop_max_pop_likelihood <= 1],
  yaxs = "i",
  yaxt = 'n',
  xlim = c(0, 2),
  ylim = c(0, 50),
  breaks = seq(from = -100, to = 100, by = 0.05),
  main = '',
  xlab = '',
  ylab = '',
  xaxt = 'n',
  col = hs_col,
  cex.axis = 2
)
axis(
  side = 4,
  at = c(0),
  labels = c(0),
  cex.axis = 1.25,
  mgp = c(2, 0.5, 0)
)
axis(1,
     pos = 0,
     cex.axis = 1.25,
     mgp = c(2, 0.5, 0))
mtext(
  side = 1,
  "Optimal escapement",
  outer = TRUE,
  padj = 1.9,
  cex = 1.3
)
mtext(
  side = 2,
  "Frequency",
  outer = TRUE,
  padj = -1.9,
  cex = 1.3
)
mtext(
  side = 1,
  "Hockey-Stick",
  outer = TRUE,
  padj = -4,
  adj = -1.6,
  at = 0.5,
  col = hs_col
)
mtext(
  side = 1,
  "Ricker",
  outer = TRUE,
  padj = -12,
  adj = -4,
  at = 0.5,
  col = ri_col
)
mtext(
  side = 1,
  "Beverton-Holt",
  outer = TRUE,
  padj = -21,
  adj = -1.5,
  at = 0.5,
  col = bh_col
)

dev.off()
# opt esc hist combined bayes ####
pdf(file = './plots_final/opt_esc_hist_combined_bayes.pdf',
    width = 7,
    height = 3.5)
par(
  mfrow = c(3, 2),
  mar = c(0, 0, 0.7, 2),
  oma = c(3.5, 3.7, 0, 2)
)
hist(
  bh_opt_esc_prop_fitted_k_bayes,
  yaxs = "i",
  yaxt = 'n',
  xlim = c(0, 1),
  ylim = c(0, 5000000),
  breaks = seq(from = -100, to = 100, by = 0.025),
  main = '',
  xlab = '',
  ylab = '',
  cex.lab = 1.5,
  xaxt = 'n',
  col = bh_col,
  cex.axis = 2
)
text(0.25, 4000000, 'a) Proportion of fitted k', cex = 1.2)
axis(
  side = 2,
  at = c(0, 75),
  labels = c(0, 75),
  cex.axis = 1.25,
  mgp = c(2, 0.5, 0)
)
axis(1,
     pos = 0,
     cex.axis = 1.25,
     mgp = c(2, 0.5, 0))
lines(c(-1, 1), c(0, 0))
hist(
  bh_opt_esc_prop_max_pop_bayes,
  yaxs = "i",
  yaxt = 'n',
  xlim = c(0, 2),
  ylim = c(0, 2000000),
  breaks = seq(from = -100, to = 100, by = 0.05),
  main = '',
  xlab = '',
  ylab = '',
  cex.lab = 1.5,
  xaxt = 'n',
  col = bh_col,
  cex.axis = 2
)
axis(
  side = 4,
  at = c(0, 50),
  labels = c(0, 50),
  cex.axis = 1.25,
  mgp = c(2, 0.5, 0)
)
axis(1,
     pos = 0,
     cex.axis = 1.25,
     mgp = c(2, 0.5, 0))
lines(c(-1, 1), c(0, 0))
text(0.7, 1600000, 'b) Proportion of max biomass', cex = 1.2)
hist(
  rk_opt_esc_prop_fitted_k_bayes,
  yaxs = "i",
  yaxt = 'n',
  xlim = c(0, 1),
  ylim = c(0, 8000000),
  breaks = seq(from = -100, to = 100, by = 0.025),
  main = '',
  xlab = '',
  ylab = '',
  cex.lab = 1.5,
  xaxt = 'n',
  col = ri_col,
  cex.axis = 2
)
axis(
  side = 2,
  at = c(0, 150),
  labels = c(0, 150),
  cex.axis = 1.25,
  mgp = c(2, 0.5, 0)
)
axis(1,
     pos = 0,
     cex.axis = 1.25,
     mgp = c(2, 0.5, 0))
hist(
  rk_opt_esc_prop_max_pop_bayes,
  yaxs = "i",
  yaxt = 'n',
  xlim = c(0, 2),
  ylim = c(0, 2000000),
  breaks = seq(from = -100, to = 100, by = 0.05),
  main = '',
  xlab = '',
  ylab = '',
  cex.lab = 1.5,
  xaxt = 'n',
  col = ri_col,
  cex.axis = 2
)
axis(
  side = 4,
  at = c(0, 50),
  labels = c(0, 50),
  cex.axis = 1.25,
  mgp = c(2, 0.5, 0)
)
axis(1,
     pos = 0,
     cex.axis = 1.25,
     mgp = c(2, 0.5, 0))
hist(
  hs_opt_esc_prop_fitted_k_bayes,
  yaxs = "i",
  yaxt = 'n',
  xlim = c(0, 1),
  ylim = c(0, 1500000),
  breaks = seq(from = -100, to = 100, by = 0.025),
  main = '',
  xlab = '',
  ylab = '',
  xaxt = 'n',
  col = hs_col,
  cex.axis = 2
)
axis(
  side = 2,
  at = c(0, 25),
  labels = c(0, 25),
  cex.axis = 1.25,
  mgp = c(2, 0.5, 0)
)
axis(1,
     pos = 0,
     cex.axis = 1.25,
     mgp = c(2, 0.5, 0))
hist(
  hs_opt_esc_prop_max_pop_bayes,
  yaxs = "i",
  yaxt = 'n',
  xlim = c(0, 2),
  ylim = c(0, 2200000),
  breaks = seq(from = -100, to = 100, by = 0.05),
  main = '',
  xlab = '',
  ylab = '',
  xaxt = 'n',
  col = hs_col,
  cex.axis = 2
)
axis(
  side = 4,
  at = c(0, 50),
  labels = c(0, 50),
  cex.axis = 1.25,
  mgp = c(2, 0.5, 0)
)
axis(1,
     pos = 0,
     cex.axis = 1.25,
     mgp = c(2, 0.5, 0))
mtext(
  side = 1,
  "Optimal escapement",
  outer = TRUE,
  padj = 1.9,
  cex = 1.3
)
mtext(
  side = 2,
  "Frequency",
  outer = TRUE,
  padj = -1.9,
  cex = 1.3
)
mtext(
  side = 1,
  "Hockey-Stick",
  outer = TRUE,
  padj = -4,
  adj = -1.6,
  at = 0.5,
  col = hs_col
)
mtext(
  side = 1,
  "Ricker",
  outer = TRUE,
  padj = -12,
  adj = -4,
  at = 0.5,
  col = ri_col
)
mtext(
  side = 1,
  "Beverton-Holt",
  outer = TRUE,
  padj = -21,
  adj = -1.5,
  at = 0.5,
  col = bh_col
)

dev.off()

