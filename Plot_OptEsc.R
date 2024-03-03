# plotting opt esc

# reading in data ####
bh_r <- readRDS("./opt_esc_data/bh_r.rds")
hs_r <- readRDS("./opt_esc_data/hs_r.rds")
rk_r <- readRDS("./opt_esc_data/rk_r.rds")

bh_k <- readRDS("./opt_esc_data/bh_k.rds")
hs_k <- readRDS("./opt_esc_data/hs_k.rds")
rk_k <- readRDS("./opt_esc_data/rk_k.rds")

bh_bf <- readRDS("./opt_esc_data/bh_k.rds")
hs_bf <- readRDS("./opt_esc_data/hs_k.rds")
rk_bf <- readRDS("./opt_esc_data/rk_k.rds")

bh_opt_esc <- readRDS("./opt_esc_data/bh_opt_esc.rds")
hs_opt_esc <- readRDS("./opt_esc_data/hs_opt_esc.rds")
rk_opt_esc <- readRDS("./opt_esc_data/rk_opt_esc.rds")

bh_opt_esc_scaled <- readRDS("./opt_esc_data/bh_opt_esc_scaled.rds")
hs_opt_esc_scaled <- readRDS("./opt_esc_data/hs_opt_esc_scaled.rds")
rk_opt_esc_scaled <- readRDS("./opt_esc_data/rk_opt_esc_scaled.rds")

bh_opt_esc_bf <- readRDS("./opt_esc_data/bh_opt_esc_bf.rds")
hs_opt_esc_bf <- readRDS("./opt_esc_data/hs_opt_esc_bf.rds")
rk_opt_esc_bf <- readRDS("./opt_esc_data/rk_opt_esc_bf.rds")

bh_opt_esc_scaled_bf <- readRDS("./opt_esc_data/bh_opt_esc_scaled_bf.rds")
hs_opt_esc_scaled_bf <- readRDS("./opt_esc_data/hs_opt_esc_scaled_bf.rds")
rk_opt_esc_scaled_bf <- readRDS("./opt_esc_data/rk_opt_esc_scaled_bf.rds")

bh_opt_esc_sig_prob <- readRDS("./opt_esc_data/bh_opt_esc_sig_prob.rds")
hs_opt_esc_sig_prob <- readRDS("./opt_esc_data/hs_opt_esc_sig_prob.rds")
rk_opt_esc_sig_prob <- readRDS("./opt_esc_data/rk_opt_esc_sig_prob.rds")

bh_opt_esc_scaled_sig_prob <- readRDS("./opt_esc_data/bh_opt_esc_scaled_sig_prob.rds")
hs_opt_esc_scaled_sig_prob <- readRDS("./opt_esc_data/hs_opt_esc_scaled_sig_prob.rds")
rk_opt_esc_scaled_sig_prob <- readRDS("./opt_esc_data/rk_opt_esc_scaled_sig_prob.rds")

bh_opt_esc_pi <- readRDS("./opt_esc_data/bh_opt_esc_pi.rds")
hs_opt_esc_pi <- readRDS("./opt_esc_data/hs_opt_esc_pi.rds")
rk_opt_esc_pi <- readRDS("./opt_esc_data/rk_opt_esc_pi.rds")

bh_opt_esc_scaled_pi <- readRDS("./opt_esc_data/bh_opt_esc_scaled_pi.rds")
hs_opt_esc_scaled_pi <- readRDS("./opt_esc_data/hs_opt_esc_scaled_pi.rds")
rk_opt_esc_scaled_pi <- readRDS("./opt_esc_data/rk_opt_esc_scaled_pi.rds")


####
hs_opt_esc_B60_count <- sum(hs_opt_esc > 0.6, na.rm = TRUE)  #percentage of the time hockey-stick optimal escapement is greater than 60% fitted k
hs_opt_esc_B60_count_pi <- sum(hs_opt_esc_pi > 0.6, na.rm = TRUE)  #percentage of the time hockey-stick optimal escapement is greater than 60% fitted k
bh_opt_esc_max_count <- hs_opt_esc_max_count <- rk_opt_esc_max_count <- 0
bh_opt_esc_scaled_max_count <- hs_opt_esc_scaled_max_count <- rk_opt_esc_scaled_max_count <- 0
opt_esc_max <- rep(NA, 284)
opt_esc_scaled_max <- rep(NA, 284)

for(i in 1:284){
  opt_esc_max[i] <- max(c(bh_opt_esc_pi[i], hs_opt_esc_pi[i], rk_opt_esc_pi[i]))
  opt_esc_scaled_max[i] <- max(c(bh_opt_esc_scaled_pi[i], hs_opt_esc_scaled_pi[i], rk_opt_esc_scaled_pi[i]))

  if(bh_opt_esc_pi[i] == opt_esc_max[i]){
    bh_opt_esc_max_count <- bh_opt_esc_max_count + 1
  }

  if(hs_opt_esc_pi[i] == opt_esc_max[i]){
    hs_opt_esc_max_count <- hs_opt_esc_max_count + 1
  }

  if(rk_opt_esc_pi[i] == opt_esc_max[i]){
    rk_opt_esc_max_count <- rk_opt_esc_max_count + 1
  }

  if(bh_opt_esc_scaled_pi[i] == opt_esc_scaled_max[i]){
    bh_opt_esc_scaled_max_count <- bh_opt_esc_scaled_max_count + 1
  }

  if(hs_opt_esc_scaled_pi[i] == opt_esc_scaled_max[i]){
    hs_opt_esc_scaled_max_count <- hs_opt_esc_scaled_max_count + 1
  }

  if(rk_opt_esc_scaled_pi[i] == opt_esc_scaled_max[i]){
    rk_opt_esc_scaled_max_count <- rk_opt_esc_scaled_max_count + 1
  }
}

# outputs 
#%model%_opt_esc_max_count = opt esc as prop of fitted k
#%model%_opt_esc_scaled = opt esc as prop of max pop


# constants ####
bh_col <- 'darkblue'
bh_col_tr <- '#4281C388'
ri_col <- 'red'
ri_col_tr <- '#E4231388'
hs_col <- 'darkgreen'
hs_col_tr <- '#00843788'
lw <- 2
ca <- 1.2
cl <- 1.3

# functions ####
# growth rates corresponding to equal optimal escapements as a proportion of fitted carrying capacity
bh_vs_hs_r <- function(r) {
  return((r - 1) / ((r ** 0.5) - 1))
}

# opt esc hist combined ####

pdf(file = './plots_final/opt_esc_hist_combined.pdf',
    width = 7,
    height = 3.5)
par(
  mfrow = c(3, 2),
  mar = c(0, 0, 0.7, 2),
  oma = c(3.5, 3.7, 0, 2)
)
hist(
  bh_opt_esc,
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
  bh_opt_esc_scaled,
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
  rk_opt_esc,
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
  rk_opt_esc_scaled,
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
  hs_opt_esc,
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
  hs_opt_esc_scaled,
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
  cex = cl
)
mtext(
  side = 2,
  "Frequency",
  outer = TRUE,
  padj = -1.9,
  cex = cl
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




# r k hist combined ####
pdf(file = './plots_final/r_k_hist_combined.pdf',
    width = 7,
    height = 3.5)
par(
  mfrow = c(3, 2),
  mar = c(0, 0, 0.7, 2),
  oma = c(3.5, 3.7, 0, 2)
)
hist(
  bh_r,
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
  bh_k,
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
  rk_r,
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
  rk_k,
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
  hs_r,
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
  hs_k,
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
  cex = cl
)
mtext(
  side = 1,
  "Carrying Capacity (k)",
  outer = TRUE,
  padj = 1.9,
  adj = 0.85,
  cex = cl
)
mtext(
  side = 2,
  "Frequency",
  outer = TRUE,
  padj = -1.9,
  cex = cl
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