# plotting growth rate and carrying capacity bayes and likelihood

# reading in data ####

# growth rate
bh_r_bayes <- readRDS("./data_r_k/bh_r_bayes.rds")
hs_r_bayes <- readRDS("./data_r_k/hs_r_bayes.rds")
rk_r_bayes <- readRDS("./data_r_k/rk_r_bayes.rds")
bh_r_likelihood <- readRDS("./data_r_k/bh_r_likelihood.rds")
hs_r_likelihood <- readRDS("./data_r_k/hs_r_likelihood.rds")
rk_r_likelihood <- readRDS("./data_r_k/rk_r_likelihood.rds")

# carrying capacity as a proportion of max pop
bh_k_bayes <- readRDS("./data_r_k/bh_k_bayes.rds")
hs_k_bayes <- readRDS("./data_r_k/hs_k_bayes.rds")
rk_k_bayes <- readRDS("./data_r_k/rk_k_bayes.rds")
bh_k_likelihood <- readRDS("./data_r_k/bh_k_likelihood.rds")
hs_k_likelihood <- readRDS("./data_r_k/hs_k_likelihood.rds")
rk_k_likelihood <- readRDS("./data_r_k/rk_k_likelihood.rds")

# constants ####
bh_col <- 'darkblue'
ri_col <- 'red'
hs_col <- 'darkgreen'

# r k hist combined bayes ####
pdf(file = './plots_final/r_k_hist_combined_bayes.pdf',
    width = 7,
    height = 3.5)
par(
  mfrow = c(3, 2),
  mar = c(0, 0, 0.7, 2),
  oma = c(3.5, 3.7, 0, 2)
)
hist(
  bh_r_bayes,
  yaxs = "i",
  yaxt = 'n',
  xlim = c(1, 4),
  ylim = c(0, 1500000),
  breaks = seq(from = -100, to = 100, by = 0.075),
  main = '',
  xlab = '',
  ylab = '',
  cex.lab = 1.5,
  xaxt = 'n',
  col = bh_col,
  cex.axis = 2
)
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
  bh_k_bayes,
  yaxs = "i",
  yaxt = 'n',
  xlim = c(0, 4),
  ylim = c(0, 1500000),
  breaks = seq(from = -100, to = 100, by = 0.1),
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
hist(
  rk_r_bayes,
  yaxs = "i",
  yaxt = 'n',
  xlim = c(1, 4),
  ylim = c(0, 1500000),
  breaks = seq(from = -100, to = 100, by = 0.075),
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
  rk_k_bayes,
  yaxs = "i",
  yaxt = 'n',
  xlim = c(0, 4),
  ylim = c(0, 1500000),
  breaks = seq(from = -100, to = 100, by = 0.1),
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
  hs_r_bayes,
  yaxs = "i",
  yaxt = 'n',
  xlim = c(1, 4),
  ylim = c(0, 2500000),
  breaks = seq(from = -100, to = 100, by = 0.075),
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
  hs_k_bayes,
  yaxs = "i",
  yaxt = 'n',
  xlim = c(0, 4),
  ylim = c(0, 1500000),
  breaks = seq(from = -100, to = 100, by = 0.1),
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
  "Growth Rate (r)",
  outer = TRUE,
  padj = 1.9,
  adj = 0.13,
  cex = 1.3
)
mtext(
  side = 1,
  "Carrying Capacity (k)",
  outer = TRUE,
  padj = 1.9,
  adj = 0.85,
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
# r k hist combined likelihood ####
pdf(file = './plots_final/r_k_hist_combined_likelihood.pdf',
    width = 7,
    height = 3.5)
par(
  mfrow = c(3, 2),
  mar = c(0, 0, 0.7, 2),
  oma = c(3.5, 3.7, 0, 2)
)
hist(
  bh_r_likelihood,
  yaxs = "i",
  yaxt = 'n',
  xlim = c(1, 4),
  ylim = c(0, 1500000/40000),
  breaks = seq(from = -100, to = 100, by = 0.075),
  main = '',
  xlab = '',
  ylab = '',
  cex.lab = 1.5,
  xaxt = 'n',
  col = bh_col,
  cex.axis = 2
)
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
  bh_k_likelihood,
  yaxs = "i",
  yaxt = 'n',
  xlim = c(0, 4),
  ylim = c(0, 1500000/30000),
  breaks = seq(from = -100, to = 100, by = 0.1),
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
hist(
  rk_r_likelihood,
  yaxs = "i",
  yaxt = 'n',
  xlim = c(1, 4),
  ylim = c(0, 1500000/40000),
  breaks = seq(from = -100, to = 100, by = 0.075),
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
  rk_k_likelihood,
  yaxs = "i",
  yaxt = 'n',
  xlim = c(0, 4),
  ylim = c(0, 1500000/30000),
  breaks = seq(from = -100, to = 100, by = 0.1),
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
  hs_r_likelihood,
  yaxs = "i",
  yaxt = 'n',
  xlim = c(1, 4),
  ylim = c(0, 2500000/40000),
  breaks = seq(from = -100, to = 100, by = 0.075),
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
  hs_k_likelihood,
  yaxs = "i",
  yaxt = 'n',
  xlim = c(0, 4),
  ylim = c(0, 1500000/10000),
  breaks = seq(from = -100, to = 100, by = 0.1),
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
  "Growth Rate (r)",
  outer = TRUE,
  padj = 1.9,
  adj = 0.13,
  cex = 1.3
)
mtext(
  side = 1,
  "Carrying Capacity (k)",
  outer = TRUE,
  padj = 1.9,
  adj = 0.85,
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