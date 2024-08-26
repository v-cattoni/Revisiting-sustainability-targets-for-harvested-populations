# Load necessary libraries
library(ggplot2)
library(reshape2)
library(viridis)
library(gridExtra)

# # Load data
# bh_r_likelihood <- readRDS('./data_r_k/bh_r_likelihood.RDS')
# bh_k_likelihood <- readRDS('./data_r_k/bh_k_likelihood.RDS')
# 
# hs_r_likelihood <- readRDS('./data_r_k/hs_r_likelihood.RDS')
# hs_k_likelihood <- readRDS('./data_r_k/hs_k_likelihood.RDS')
# 
# bh_r_bayes <- readRDS("./data_r_k/bh_r_bayes.rds")
# rk_r_bayes <- readRDS("./data_r_k/rk_r_bayes.rds")
# hs_r_bayes <- readRDS("./data_r_k/hs_r_bayes.rds")
# 
# bh_k_bayes <- readRDS("./data_r_k/bh_k_bayes.rds")
# rk_k_bayes <- readRDS("./data_r_k/rk_k_bayes.rds")
# hs_k_bayes <- readRDS("./data_r_k/hs_k_bayes.rds")
# 
# bh_opt_esc_prop_max_pop_bayes <- readRDS("./data_opt_esc/bh_opt_esc_prop_max_pop_bayes.rds")
# rk_opt_esc_prop_max_pop_bayes <- readRDS("./data_opt_esc/rk_opt_esc_prop_max_pop_bayes.rds")
# hs_opt_esc_prop_max_pop_bayes <- readRDS("./data_opt_esc/hs_opt_esc_prop_max_pop_bayes.rds")
# 
# bh_opt_esc_bayes_median_index <- rk_opt_esc_bayes_median_index <- hs_opt_esc_bayes_median_index <- rep(NA, 284)
# for(i in 1:284){
#   bh_median <- median(bh_opt_esc_prop_max_pop_bayes[((i-1)*40000):(i*40000)])
#   print(bh_median)
#   print(which(bh_opt_esc_prop_max_pop_bayes == bh_median))
#   bh_opt_esc_bayes_median_index[i] <- which(bh_opt_esc_prop_max_pop_bayes == bh_median)
# }


# Define functions
bh_opt_esc <- function(r, k) {
  return(k * (r ** 0.5 - 1) / (r - 1))
}

hs_opt_esc <- function(r, k) {
  return(k / r)
}

# Set parameters
n_vals <- 1000

rmin <- 1.0001
rmax <- 10
r_vals <- seq(from = rmin, to = rmax, length.out = n_vals)

kmin <- 0.1
kmax <- 10
k_vals <- exp(seq(log(kmin), log(kmax), length.out = n_vals))

# Initialize matrices
bh_opt_esc_mat <- matrix(NA, n_vals, n_vals)
hs_opt_esc_mat <- matrix(NA, n_vals, n_vals)

# Fill matrices
for(r in 1:n_vals) {
  for(k in 1:n_vals) {
    bh_opt_esc_mat[k, r] <- bh_opt_esc(r_vals[r], k_vals[k])
    hs_opt_esc_mat[k, r] <- hs_opt_esc(r_vals[r], k_vals[k])
  }
}

# Convert matrices to data frames for ggplot
bh_opt_esc_df <- melt(bh_opt_esc_mat)
hs_opt_esc_df <- melt(hs_opt_esc_mat)

# Add columns for r and k values
bh_opt_esc_df$r <- rep(r_vals, each = n_vals)
bh_opt_esc_df$k <- rep(k_vals, times = n_vals)
colnames(bh_opt_esc_df) <- c("Var1", "Var2", "value", "r", "k")

hs_opt_esc_df$r <- rep(r_vals, each = n_vals)
hs_opt_esc_df$k <- rep(k_vals, times = n_vals)
colnames(hs_opt_esc_df) <- c("Var1", "Var2", "value", "r", "k")

# Calculate common range for color scale
combined_min <- min(bh_opt_esc_df$value, hs_opt_esc_df$value)
combined_max <- max(bh_opt_esc_df$value, hs_opt_esc_df$value)

# Create heatmaps with likelihood points and solid black line at value = 0.6
bh_heatmap <- ggplot(bh_opt_esc_df, aes(x = r, y = k, fill = value)) +
  geom_tile() +
  scale_fill_gradient(low = "orange", high = "darkred", limits = c(combined_min, combined_max)) + # Common color scale
  scale_y_log10() + # Logarithmic scale for y-axis
  geom_hline(yintercept = 1, linetype = "dotted", color = "black") + # Horizontal dotted line at k = 1
  geom_contour(aes(z = value), color = "black", breaks = 0.6) + # Solid black line at value = 0.6
  geom_point(data = data.frame(r = bh_r_likelihood, k = bh_k_likelihood), 
             aes(x = r, y = k), color = "blue", size = 1.5, inherit.aes = FALSE, pch = 18) + 
  annotate("text", x = 10, y = 10, label = "B60", hjust = -0.1, vjust = 21, size = 3, color = "black") + # Label for contour line
  labs(title = "Beverton-Holt Optimal Escapement", x = "Growth rate (r)", y = "Carrying capacity (k) (log scale)") +
  theme_minimal() +
  theme(legend.position = "none") # Remove legend

hs_heatmap <- ggplot(hs_opt_esc_df, aes(x = r, y = k, fill = value)) +
  geom_tile() +
  scale_fill_gradient(low = "orange", high = "darkred", limits = c(combined_min, combined_max)) + # Common color scale
  scale_y_log10() + # Logarithmic scale for y-axis
  geom_hline(yintercept = 1, linetype = "dotted", color = "black") + # Horizontal dotted line at k = 1
  geom_contour(aes(z = value), color = "black", breaks = 0.6) + # Solid black line at value = 0.6
  geom_point(data = data.frame(r = hs_r_likelihood, k = hs_k_likelihood), 
             aes(x = r, y = k), color = "blue", size = 1.5, inherit.aes = FALSE, pch = 18) + 
  annotate("text", x = 10, y = 10, label = "B60", hjust = -0.1, vjust = 8, size = 3, color = "black") + # Label for contour line
  labs(title = "Hockey-Stick Optimal Escapement", x = "Growth rate (r)", y = NULL) + # Remove y-axis label
  theme_minimal() +
  theme(axis.title.y = element_blank()) # Remove y-axis label

# Combine the two heatmaps side by side and save as a PDF
pdf(file = './plots_final/opt_esc_heat_map_bh_hs_combined.pdf',
    width = 15,
    height = 7.5)

grid.arrange(bh_heatmap, hs_heatmap, ncol = 2, widths = c(0.85, 0.95), 
             padding = unit(1, "lines")) # Adjust widths and padding for space around labels

dev.off()