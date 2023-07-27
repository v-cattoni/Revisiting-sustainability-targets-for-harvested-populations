#### Reading in files ####
load("./data_files/RAMCore[asmt][v4.495].rdata") #RAM DATABASE
#this was visually identified manually recorded in the below csv file
cols <- read.csv("./data_files/non_det_cols.csv", header <- FALSE)[, 1]

#### Constants ####
years <- seq(1950, 2020, 1) #year of each row of data
num_rows <- length(years) #num rows to iterate through all data
num_cols <- length(cols) #num cols to iterate through all data
k_max <- 5 #upper bound on carrying capacity (multiple of max observed population)
r_max <- 100 #upper bound on growth rate
phi <- 0.188
col_names <- colnames(Bio)
bh_col <- 'darkblue'
bh_col_tr <- '#4281C388'
ri_col <- 'red'
ri_col_tr <- '#E4231388'
hs_col <- 'darkgreen'
hs_col_tr <- '#00843788'
pt_col <- '#DB7093'
pt_col_tr <- '#DB709388'
lw <- 2
ca <- 1.2
cl <- 1.3
                
                
#### Initialising vectors and constants ####
bio_mat <-
  matrix(NA, num_rows, num_cols) #biomass(MT)
catch_mat <-
  matrix(NA, num_rows, num_cols) #catch(MT)

r_bh <-
  c() #estimated growth rates for Beverton-Holt
r_hs <- c() #estimated growth rates for Hockey-Stick
r_ri <- c() #estimated growth rates for Ricker
r_pt <-
  c() #estimated growth rates for Pella-Tomlinson

k_bh <-
  c() #estimated carrying capacity/max pop for Beverton-Holt
k_hs <-
  c() #estimated carrying capacity/max pop for Hockey-Stick
k_ri <-
  c() #estimated carrying capacity/max pop Ricker
k_pt <-
  c() #estimated carrying capacity for Pella-Tomlinson

r2_bh <- c() #r-sqaured for Beverton-Holt
r2_hs <- c() #r-squared for Hockey-Stick
r2_ri <- c() #r-squared for Ricker
r2_pt <- c() #r-squared for Pella-Tomlinson

l <- c() #number of data points in each data set

num_bh <-
  0 #number of data sets where Beverton-Holt is best fit
num_hs <-
  0 #number of data sets where Hockey-Stick is best fit
num_ri <-
  0 #number of data sets where Ricker is best fit
num_pt <-
  0 #number of data sets where Pella-Tomlinson is best fit

num_bh_rs <-
  c() #list containing number of times Beverton-Holt is best for with increasing upper bound on r (used for cumulative frequency plot)
num_hs_rs <-
  c() #list containing number of times Hockey-Stick is best for with increasing upper bound on r (used for cumulative frequency plot)
num_ri_rs <-
  c() #list containing number of times Ricker is best for with increasing upper bound on r (used for cumulative frequency plot)
num_pt_rs <-
  c() #list containing number of times Pella-Tomlinson is best for with increasing upper bound on r (used for cumulative frequency plot)

opt_esc_bh <-
  c() #estimated optimal escapement for Beverton-Holt (k <- 1)
opt_esc_hs <-
  c() #estimated optimal escapement for Hockey-Stick (k <- 1)
opt_esc_ri <-
  c() #estimated optimal escapement for Ricker (k <- 1)
opt_esc_pt <-
  c() #estimated optimal escapement for Pella-Tomlinson (k <- 1)

opt_esc_bh_bf <-
  c() #estimated optimal escapement for Beverton-Holt (k <- 1) (best fit)
opt_esc_hs_bf <-
  c() #estimated optimal escapement for Hockey-Stick (k <- 1) (best fit)
opt_esc_ri_bf <-
  c() #estimated optimal escapement for Ricker (k <- 1) (best fit)
opt_esc_pt_bf <-
  c() #estimated optimal escapement for Pella-Tomlinson (k <- 1) (best fit)

opt_esc_bh_scaled <-
  c() #estimated optimal escapement for Beverton-Holt (proportion of max pop)
opt_esc_hs_scaled <-
  c() #estimated optimal escapement for Hockey-Stick (proportion of max pop)
opt_esc_ri_scaled <-
  c() #estimated optimal escapement for Ricker (proportion of max pop)
opt_esc_pt_scaled <-
  c() #estimated optimal escapement for Pella-Tomlinson (proportion of max pop)

opt_esc_bh_scaled_bf <-
  c() #estimated optimal escapement for Beverton-Holt (proportion of max pop) (best fit)
opt_esc_hs_scaled_bf <-
  c() #estimated optimal escapement for Hockey-Stick (proportion of max pop) (best fit)
opt_esc_ri_scaled_bf <-
  c() #estimated optimal escapement for Ricker (proportion of max pop) (best fit)
opt_esc_pt_scaled_bf <-
  c() #estimated optimal escapement for Pella-Tomlinson (proportion of max pop) (best fit)
opt_esc_data <-
  c() #records optimal escapement of the data

num_hs_lower_bh_r <- 0
num_hs_lower_bh_k <- 0
num_hs_higher_bh_opt_esc_prop_k <- 0
num_hs_higher_bh_opt_esc_prop_max_pop <- 0
num_hs_higher_b60_opt_esc_prop_k <- 0
num_bh_k_higher_1 <- 0
num_hs_higher_bh_opt_esc_replacing_1s <- 0
num_hs_1 <- 0

#### Functions ####

#Beverton-Holt
BH <- function(r, k = 1, x) {
  return(r * x / (1 + (r - 1) * x / k))
}

#Hockey-Stick
HS <- function(r, k = 1, x) {
  return(apply(cbind(r * x, k), 1, min))
}

#Ricker
RI <- function(r, k = 1, x) {
  return(x * r ^ (1 - x / k))
}

#Pella-Tomlinson
PT <- function(r, k = 1, x , phi = 0.188) {
  return(x + x * (r - 1) * (1 - (x / k) ^ phi))
}


#Beverton-Holt optimal escapement
bh_opt_esc <- function(r, k) {
  return(k * (r ** 0.5 - 1) / (r - 1))
}

#Hockey-Stick optimal escpaement
hs_opt_esc <- function(r, k) {
  return(k / r)
}

#Ricker optimal escapement
# Define the function to find the root of
derivative_ri <- function(X) {
  # Calculate the derivative of ri with respect to X
  derivative <-
    (p ^ (1 - (X / k))) * (k - X * log(p)) / k - 1
  
  # Return the derivative
  return(derivative)
}

#Pella-Tomlinson optimal escapement
pt_opt_esc <- ((1 / (1 + phi)) ^ (1 / phi))


#Sum of Squares Error
sse <- function(X, Y) {
  s <- 0
  for (i in 1:length(X)) {
    s <- s + (X[i] - Y[i]) ** 2
  }
  return(s)
}

#Sum of Squares Total
sst <- function(Y) {
  s <- 0
  ym <- mean(Y)
  for (i in 1:length(Y)) {
    s <- s + (Y[i] - ym) ** 2
  }
  return(s)
}

#R squared
r_squared <- function(X, Y) {
  return(1 - sse(X, Y) / sst(Y))
}

#Sum of Squares Error (Beverton-Holt)
#This is the function which will optimised in terms of r and k
sse_bh <- function(p, X, Y) {
  r <- p[1]
  k <- p[2]
  return(sse(Y, BH(r, k, X)))
}

#Sum of Squares Error (Hockey-Stick)
#This is the function which will optimised in terms of r and k
sse_hs <- function(p, X, Y) {
  r <- p[1]
  k <- p[2]
  return(sse(Y, HS(r, k, X)))
}

#Sum of Squares Error (Ricker)
#This is the function which will optimised in terms of r and k
sse_ri <- function(p, X, Y) {
  r <- p[1]
  k <- p[2]
  return(sse(Y, RI(r, k, X)))
}

#Sum of Squares Error (Pella-Tomlinson)
#This is the function which will optimised in terms of r and k
sse_pt <- function(p, X, Y) {
  r <- p[1]
  k <- p[2]
  return(sse(Y, PT(r, k, X)))
}

#Growth rates corresponding to equal optimal escapements as a proportion of fitted carrying capacity
bh_vs_hs_r <- function(r) {
  return((r - 1) / ((r ** 0.5) - 1))
}

#### Sorting data ####
#Creates a matrix of Biomass data and a matrix of Catch data with only
#overlapping entries.
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


#### Main ####
for (i in 1:length(cols)) {
  bios <- bio_mat[, i][!is.na(bio_mat[, i])]
  catchs <- catch_mat[, i][!is.na(catch_mat[, i])]
  X <-
    bios[1:(length(bios) - 1)] - catchs[1:(length(bios) - 1)]
  Y <- bios[2:length(bios)]
  
  #intial fact_guess before optimisation
  init_guess <- c(mean(Y / X), max(X, Y))
  fact_guess <- c(0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75)
  num_guess <- length(fact_guess)
  
  #initialising matrices to store data for each guess
  bh_r2_guess_mat <- matrix(NA, num_guess, num_guess)
  hs_r2_guess_mat <- matrix(NA, num_guess, num_guess)
  ri_r2_guess_mat <- matrix(NA, num_guess, num_guess)
  pt_r2_guess_mat <- matrix(NA, num_guess, num_guess)
  bh_r_guess_mat <- matrix(NA, num_guess, num_guess)
  hs_r_guess_mat <- matrix(NA, num_guess, num_guess)
  ri_r_guess_mat <- matrix(NA, num_guess, num_guess)
  pt_r_guess_mat <- matrix(NA, num_guess, num_guess)
  bh_k_guess_mat <- matrix(NA, num_guess, num_guess)
  hs_k_guess_mat <- matrix(NA, num_guess, num_guess)
  ri_k_guess_mat <- matrix(NA, num_guess, num_guess)
  pt_k_guess_mat <- matrix(NA, num_guess, num_guess)
  
  
  for (j in 1:num_guess) {
    for (k in 1:num_guess) {
      jk_guess <-
        c(init_guess[1] * fact_guess[j], init_guess[2] * fact_guess[k])
      #parameters estimated
      bh_opt <-
        optim(
          jk_guess,
          sse_bh,
          X = X,
          Y = Y,
          method = "L-BFGS-B",
          lower = c(1.00001, mean(Y)),
          upper = c(min(r_max, max(Y) / min(X)), max(Y) * k_max),
          control = c(factr <-
                        1e-15, factr <-
                        1e-15 * max(Y) / mean(Y / X))
        )
      hs_opt <-
        optim(
          jk_guess,
          sse_hs,
          X = X,
          Y = Y,
          method = "L-BFGS-B",
          lower = c(1.00001, mean(Y)),
          upper = c(min(r_max, max(Y) / min(X)), max(Y) * k_max),
          control = c(factr <-
                        1e-15, factr <-
                        1e-15 * max(Y) / mean(Y / X))
        )
      ri_opt <-
        optim(
          jk_guess,
          sse_ri,
          X = X,
          Y = Y,
          method = "L-BFGS-B",
          lower = c(1.00001, mean(Y)),
          upper = c(min(r_max, max(Y) / min(X)), max(Y) * k_max),
          control = c(factr <-
                        1e-15, factr <-
                        1e-15 * max(Y) / mean(Y / X))
        )
      pt_opt <-
        optim(
          jk_guess,
          sse_pt,
          X = X,
          Y = Y,
          method = "L-BFGS-B",
          lower = c(1.00001, mean(Y)),
          upper = c(min(r_max, max(Y) / min(X)), max(Y) * k_max),
          control = c(factr <-
                        1e-15, factr <-
                        1e-15 * max(Y) / mean(Y / X))
        )
      
      #estimated parameters stored
      if (hs_opt[[1]][1] > hs_opt[[1]][2] / min(X)) {
        hs_opt[[1]][1] <- hs_opt[[1]][2] / min(X)
      }
      if (hs_opt[[1]][2] > max(X, Y)) {
        hs_opt[[1]][2] <- max(X, Y)
      }
      
      bh_r_guess_mat[j, k] <- bh_opt[[1]][1]
      hs_r_guess_mat[j, k] <- hs_opt[[1]][1]
      ri_r_guess_mat[j, k] <- ri_opt[[1]][1]
      pt_r_guess_mat[j, k] <- pt_opt[[1]][1]
      
      bh_k_guess_mat[j, k] <- bh_opt[[1]][2]
      hs_k_guess_mat[j, k] <- hs_opt[[1]][2]
      ri_k_guess_mat[j, k] <- ri_opt[[1]][2]
      pt_k_guess_mat[j, k] <- pt_opt[[1]][2]
      
      bh_r2_guess_mat[j, k] <-
        r_squared(BH(bh_opt[[1]][1], bh_opt[[1]][2], X), Y)
      hs_r2_guess_mat[j, k] <-
        r_squared(HS(hs_opt[[1]][1], hs_opt[[1]][2], X), Y)
      ri_r2_guess_mat[j, k] <-
        r_squared(RI(ri_opt[[1]][1], ri_opt[[1]][2], X), Y)
      pt_r2_guess_mat[j, k] <-
        r_squared(PT(pt_opt[[1]][1], pt_opt[[1]][2], X), Y)
    }
  }
  
  bh_jk <-
    which(bh_r2_guess_mat == max(bh_r2_guess_mat), arr.ind <-
            TRUE)[1, ]
  hs_jk <-
    which(hs_r2_guess_mat == max(hs_r2_guess_mat), arr.ind <-
            TRUE)[1, ]
  ri_jk <-
    which(ri_r2_guess_mat == max(ri_r2_guess_mat), arr.ind <-
            TRUE)[1, ]
  pt_jk <-
    which(pt_r2_guess_mat == max(pt_r2_guess_mat), arr.ind <-
            TRUE)[1, ]
  
  bh_opt_r <- bh_r_guess_mat[bh_jk[1],][bh_jk[2]]
  hs_opt_r <- hs_r_guess_mat[hs_jk[1],][hs_jk[2]]
  ri_opt_r <- ri_r_guess_mat[ri_jk[1],][ri_jk[2]]
  pt_opt_r <- pt_r_guess_mat[pt_jk[1],][pt_jk[2]]
  
  bh_opt_k <- bh_k_guess_mat[bh_jk[1],][bh_jk[2]]
  hs_opt_k <- hs_k_guess_mat[hs_jk[1],][hs_jk[2]]
  ri_opt_k <- ri_k_guess_mat[ri_jk[1],][ri_jk[2]]
  pt_opt_k <- pt_k_guess_mat[pt_jk[1],][pt_jk[2]]
  
  
  #storing r-squared values
  r2_bh <-
    append(r2_bh, max(bh_r2_guess_mat))
  r2_hs <-
    append(r2_hs, max(hs_r2_guess_mat))
  r2_ri <-
    append(r2_ri, max(ri_r2_guess_mat))
  r2_pt <-
    append(r2_pt, max(pt_r2_guess_mat))
  
  
  #storing number of data points
  l <- append(l, length(X))
  
  #storing carrying capacity/max pops
  k_bh <- append(k_bh, bh_opt_k / max(Y))
  k_hs <- append(k_hs, hs_opt_k / max(Y))
  k_ri <- append(k_ri, ri_opt_k / max(Y))
  k_pt <- append(k_pt, pt_opt_k / max(Y))
  
  #storing growth rates
  r_bh <- append(r_bh, bh_opt_r)
  r_hs <- append(r_hs, hs_opt_r)
  r_ri <- append(r_ri, ri_opt_r)
  r_pt <- append(r_pt, pt_opt_r)
  
  #stores opitmal escapement as a proportion of the fitted K
  opt_esc_bh <- append(opt_esc_bh, bh_opt_esc(r_bh[i], 1))
  opt_esc_hs <- append(opt_esc_hs, hs_opt_esc(r_hs[i], 1))
  p <- r_ri[i]
  k <- 1
  result <- uniroot(derivative_ri, lower = 0, upper = 10 * k)
  opt_esc_ri <- append(opt_esc_ri, result$root)
  opt_esc_pt <- append(opt_esc_pt, (1 / (1 + phi)) ^ (1 / phi))
  
  #stores optimal escapement as a proportion of the max population
  opt_esc_bh_scaled <-
    append(opt_esc_bh_scaled, bh_opt_esc(r_bh[i], k_bh[i]))
  opt_esc_hs_scaled <-
    append(opt_esc_hs_scaled, hs_opt_esc(r_hs[i], k_hs[i]))
  p <- r_ri[i]
  k <- k_ri[i]
  result <- uniroot(derivative_ri, lower = 0, upper = 10 * k)
  opt_esc_ri_scaled <- append(opt_esc_ri_scaled, result$root)
  opt_esc_pt_scaled <-
    append(opt_esc_pt_scaled, k_pt[i] * (1 / (1 + phi)) ^
             (1 / phi))
  
  
  #counts the number of times each model is best-fit
  if (r2_bh[i] == max(c(r2_bh[i], r2_ri[i], r2_hs[i], r2_pt[i]))) {
    num_bh <- num_bh + 1
    opt_esc_bh_bf <-
      append(opt_esc_bh_bf, opt_esc_bh[i])
    opt_esc_bh_scaled_bf <-
      append(opt_esc_bh_scaled_bf, opt_esc_bh_scaled[i])
  }
  if (r2_ri[i] == max(c(r2_bh[i], r2_ri[i], r2_hs[i], r2_pt[i]))) {
    num_ri <- num_ri + 1
    opt_esc_ri_bf <-
      append(opt_esc_ri_bf, opt_esc_ri[i])
    opt_esc_ri_scaled_bf <-
      append(opt_esc_ri_scaled_bf, opt_esc_ri_scaled[i])
  }
  if (r2_hs[i] == max(c(r2_bh[i], r2_ri[i], r2_hs[i], r2_pt[i]))) {
    num_hs <- num_hs + 1
    opt_esc_hs_bf <-
      append(opt_esc_hs_bf, opt_esc_hs[i])
    opt_esc_hs_scaled_bf <-
      append(opt_esc_hs_scaled_bf, opt_esc_hs_scaled[i])
  }
  if (r2_pt[i] == max(c(r2_bh[i], r2_ri[i], r2_hs[i], r2_pt[i]))) {
    num_pt <- num_pt + 1
    opt_esc_pt_bf <-
      append(opt_esc_pt_bf, opt_esc_pt[i])
    opt_esc_pt_scaled_bf <-
      append(opt_esc_pt_scaled_bf, opt_esc_pt_scaled[i])
  }
  
  if (hs_opt_r < bh_opt_r) {
    num_hs_lower_bh_r = num_hs_lower_bh_r + 1
  }
  if (hs_opt_k < bh_opt_k) {
    num_hs_lower_bh_k = num_hs_lower_bh_k + 1
  }
  if (opt_esc_bh[i] < opt_esc_hs[i]) {
    num_hs_higher_bh_opt_esc_prop_k = num_hs_higher_bh_opt_esc_prop_k + 1
  }
  if (opt_esc_hs[i] > 0.6) {
    num_hs_higher_b60_opt_esc_prop_k = num_hs_higher_b60_opt_esc_prop_k + 1
  }
  if (opt_esc_hs_scaled[i] > opt_esc_bh_scaled[i]) {
    num_hs_higher_bh_opt_esc_prop_max_pop = num_hs_higher_bh_opt_esc_prop_max_pop + 1
  }
  if (k_bh[i] > 1) {
    num_bh_k_higher_1 = num_bh_k_higher_1 + 1
  }
  
  if (k_hs[i] > 0.99999 &
      opt_esc_hs_scaled[i] < opt_esc_bh_scaled[i]) {
    if (hs_opt_esc(r_hs[i], k_bh[i]) > opt_esc_bh_scaled[i]) {
      num_hs_higher_bh_opt_esc_replacing_1s = num_hs_higher_bh_opt_esc_replacing_1s + 1
    }
  }
}

#### Results ####
print("Results for non-deterministic data set:")
cat("\n")
cat("\n")
print(paste0("Number of data sets = ", num_cols))
cat("\n")
print("Best fit as a percentage.")
print(paste0(
  "BH: ",
  round(100 * num_bh / (num_cols)),
  "%, HS: ",
  round(100 * num_hs / (num_cols)),
  "%, RI: ",
  round(100 * num_ri / (num_cols)),
  "%, PT: ",
  round(100 * num_pt / (num_cols)),
  "%"
))
cat("\n")
print("Median r-squared values.")
print(paste0(
  "BH: ",
  round(median(r2_bh), digits <-
          2),
  ", HS: ",
  round(median(r2_hs), digits <-
          2),
  ", RI: ",
  round(median(r2_ri), digits <-
          2),
  ", PT: ",
  round(median(r2_pt), digits <- 2)
))
cat("\n")
print(paste0(
  "HS suggests a growth rate lower than BH ",
  round(100 * num_hs_lower_bh_r / (num_cols)),
  "% of the time."
))
cat("\n")
print(paste0(
  "HS suggests a carrying capacity lower than BH ",
  round(100 * num_hs_lower_bh_k / (num_cols)),
  "% of the time."
))
cat("\n")
print(
  paste0(
    "HS suggests higher optimal escapement than BH as a proportion of k ",
    round(100 * num_hs_higher_bh_opt_esc_prop_k / (num_cols)),
    "% of the time."
  )
)
cat("\n")
print(paste0(
  "HS suggests escapement higher than B60 ",
  round(100 * num_hs_higher_b60_opt_esc_prop_k / (num_cols)) ,
  "% of the time."
))
cat("\n")
print(
  paste0(
    "HS suggests higher optimal escapement than BH as a proportion of the max pop ",
    round(
      100 * num_hs_higher_bh_opt_esc_prop_max_pop / (num_cols)
    )  ,
    "% of the time."
  )
)
cat("\n")
print(
  paste0(
    "This increases to ",
    round(
      100 * (
        num_hs_higher_bh_opt_esc_prop_max_pop + num_hs_higher_bh_opt_esc_replacing_1s
      ) / (num_cols)
    ),
    "% if we take the BH k as true when HS k = 1."
  )
)
cat("\n")
print(
  paste0(
    "BH suggests a carrying capacity greater than the max pop ",
    round(100 * num_bh_k_higher_1 / (num_cols)),
    "% of the time."
  )
)
cat("\n")
cat("\n")
#### Fig 1 - Recruitment functions ####
pdf(file = './plots/recruitment_functions.pdf',
    width = 3.5,
    height = 3.5)
par(
  mfrow = c(1, 1),
  mar = c(1, 1, .7, .7),
  oma = c(2.51, 2.51, 0, 0)
)

r = 1.7

x = seq(0, 1.6, by = .005)
plot(
  x,
  HS(r = r, x = x),
  ylim = c(0, 1.2),
  xlab = "",
  ylab = "",
  lty = 2,
  type = 'l',
  xaxs = 'i',
  yaxs = 'i',
  cex.axis = ca,
  lwd = lw * 1.2,
  col = hs_col
)
lines(
  x,
  BH(r = r, x = x),
  ylim = c(0, 1),
  xlab = "",
  ylab = "",
  lty = 1,
  type = 'l',
  xaxs = 'i',
  yaxs = 'i',
  cex.axis = ca,
  lwd = lw,
  col = bh_col
)
lines(
  x,
  RI(r = r, x = x),
  ylim = c(0, 1),
  xlab = "",
  ylab = "",
  lty = 3,
  type = 'l',
  xaxs = 'i',
  yaxs = 'i',
  cex.axis = ca,
  lwd = lw * 1.7,
  col = ri_col
)
lines(
  x,
  PT(r = r, x = x),
  ylim = c(0, 1),
  xlab = "",
  ylab = "",
  lty = 4,
  type = 'l',
  xaxs = 'i',
  yaxs = 'i',
  cex.axis = ca,
  lwd = lw * 1.4,
  col = pt_col
)
legend(
  'bottomright',
  c('Beverton-Holt',
    'Hockey-stick',
    'Ricker',
    'Pella-Tomlinson'),
  col = c(bh_col, hs_col, ri_col, pt_col),
  lty = c(1, 2, 3, 4),
  bty = 'n',
  cex = .9,
  lwd = c(lw, lw * 1.2, lw * 1.7, lw * 1.5)
)
mtext(
  side = 1,
  "Population size, year t",
  outer = TRUE,
  padj = 1.9,
  cex = cl
)
mtext(
  side = 2,
  "Population size, year t+1",
  outer = TRUE,
  padj = -1.9,
  cex = cl
)
dev.off()



#### Fig 2 - Optimal escapement functions ####
pdf(file = './plots/opt_esc_curves.pdf',
    width = 3.5,
    height = 3.5)
par(
  mfrow = c(1, 1),
  mar = c(1, 1, .7, .7),
  oma = c(2.51, 2.51, 0, 0)
)

r = 1.7

x = seq(0, 1.6, by = .005)
R = seq(1.01, 4, by = .01)
e.B = (sqrt(R) - 1) / (R - 1)
e.H = 1 / R
B60 = rep(.6, length.out = length(R))

ri_opt_xs = c()
for (i in R) {
  p = i
  result =
    uniroot(derivative_ri, lower = 0, upper = 10)
  ri_opt_xs = append(ri_opt_xs, result$root)
}

plot(
  R,
  1 / R,
  ylim = c(0, 1.05),
  xlab = "",
  ylab = "",
  lty = 2,
  type = 'l',
  xaxs = 'i',
  yaxs = 'i',
  cex.axis = ca,
  lwd = lw * 1.2,
  col = hs_col
)
lines(
  R,
  e.B,
  ylim = c(0, 1),
  xlab = "",
  ylab = "",
  lty = 1,
  type = 'l',
  xaxs = 'i',
  yaxs = 'i',
  cex.axis = ca,
  lwd = lw,
  col = bh_col
)
lines(
  R,
  ri_opt_xs,
  ylim = c(0, 1),
  xlab = "",
  ylab = "",
  lty = 3,
  type = 'l',
  xaxs = 'i',
  yaxs = 'i',
  cex.axis = ca,
  lwd = lw * 1.7,
  col = ri_col
)
lines(
  R,
  rep(0.4, length(R)),
  ylim = c(0, 1),
  xlab = "",
  ylab = "",
  lty = 4,
  type = 'l',
  xaxs = 'i',
  yaxs = 'i',
  cex.axis = ca,
  lwd = lw * 1.4,
  col = pt_col
)
legend(
  'topright',
  c('Beverton-Holt',
    'Hockey-stick',
    'Ricker',
    'Pella-Tomlinson'),
  col = c(bh_col, hs_col, ri_col, pt_col),
  lty = c(1, 2, 3, 4),
  bty = 'n',
  cex = .9,
  lwd = c(lw, lw * 1.2, lw * 1.7, lw * 1.5)
)
mtext(
  side = 1,
  "Population growth multiplier, r",
  outer = TRUE,
  padj = 1.9,
  cex = cl
)
mtext(
  side = 2,
  "Optimal Escapement",
  outer = TRUE,
  padj = -1.9,
  cex = cl
)
dev.off()

#### Fig 3 - Catch when model is wrong ####

# hockey-stick true
#bev-holt decision given hockey-stick is true, biomass & catch
x2.B.H = HS(r = R, x = e.B)
C.B.H = x2.B.H - e.B
#B60 decision given hockey-stick is true, biomass & catch
x2.60.H = HS(r = R, x = B60)
C.60.H = x2.60.H - B60
#1/r hockey-stick decision given hockey-stick is true, biomass & catch
x2.H.H = HS(r = R, x = e.H)
C.H.H = x2.H.H - e.H

# Bev-Holt true
#bev-holt decision given Bev-Holt is true, biomass & catch
x2.B.B = BH(x = e.B, r = R)
C.B.B = x2.B.B - e.B
#B60 decision given Bev-Holt is true, biomass & catch
x2.60.B = BH(x = B60, r = R)
C.60.B = x2.60.B - B60
#1/r hockey-stick decision given Bev-Holt is true, biomass & catch
x2.H.B = BH(x = e.H, r = R)
C.H.B = x2.H.B - e.H

pdf(file = './plots/catch.pdf',
    width = 7,
    height = 3.5)
par(
  mfrow = c(1, 2),
  mar = c(1, 1, .5, .55),
  oma = c(2.51, 2.51, 0, 0)
)
ylimC = c(0, max(C.B.B, C.H.H))
ylimB = c(0, 1.01)
#BH true
plot(
  R,
  C.B.B,
  ylim = ylimC,
  xlab = "",
  ylab = "",
  type = 'l',
  xaxs = 'i',
  yaxs = 'i',
  cex.axis = ca,
  lwd = lw,
  col = bh_col
)
#lines(R, C.60.B, lty = 2, lwd = lw)
lines(R,
      C.H.B,
      lty = 2,
      lwd = lw * 1.5,
      col = hs_col)
text(1.92, .7, 'a) Beverton-Holt true')
#HS true
plot(
  R,
  C.H.H,
  ylim = ylimC,
  xlab = "",
  ylab = "",
  yaxt = 'n',
  type = 'l',
  lty = 2,
  col = hs_col,
  xaxs = 'i',
  yaxs = 'i',
  cex.axis = ca,
  lwd = lw * 1.5
)
#lines(R, C.60.H, lty = 2, lwd = lw)
lines(R,
      C.B.H,
      lty = 1,
      lwd = lw,
      col = bh_col)
text(1.92, .7, 'b) Hockey-Stick true')
mtext(
  side = 1,
  "Population growth multiplier, r",
  outer = TRUE,
  padj = 1.9,
  cex = cl
)
mtext(
  side = 2,
  "Catch",
  outer = TRUE,
  padj = -1.9,
  cex = cl
)
legend(
  'bottomright',
  c('Beverton-Holt harvest', 'Hockey-stick harvest'),
  col = c(bh_col, hs_col),
  lty = c(1, 2),
  bty = 'n',
  cex = .9,
  lwd = c(lw, lw * 1.2)
)
dev.off()


#### Fig 4 - Proportional catch when model is wrong ####
pdf(file = './plots/prop_catch.pdf',
    width = 3.5,
    height = 3.5)
par(
  mfrow = c(1, 1),
  mar = c(1, 1, .7, .7),
  oma = c(2.51, 2.51, 0, 0)
)

plot(
  R,
  C.H.B / C.B.B,
  ylim = c(0, 1),
  xlab = "",
  ylab = "",
  lty = 2,
  type = 'l',
  xaxs = 'i',
  yaxs = 'i',
  cex.axis = ca,
  lwd = lw * 1.2,
  col = hs_col
)
lines(
  R,
  C.B.H / C.H.H,
  ylim = c(0, 1),
  xlab = "",
  ylab = "",
  lty = 1,
  type = 'l',
  xaxs = 'i',
  yaxs = 'i',
  cex.axis = ca,
  lwd = lw,
  col = bh_col
)
legend(
  'bottomright',
  c('Beverton-Holt harvest', 'Hockey-stick harvest'),
  col = c(bh_col, hs_col),
  lty = c(1, 2),
  bty = 'n',
  cex = .9,
  lwd = c(lw, lw * 1.2)
)
mtext(
  side = 1,
  "Population growth multiplier, r",
  outer = TRUE,
  padj = 1.9,
  cex = cl
)
mtext(
  side = 2,
  "Proportion of optimal catch",
  outer = TRUE,
  padj = -1.9,
  cex = cl
)
dev.off()


#### Fig 5 - Growth rate and carrying capacity scatter plots####
R = seq(1, 5, length.out = 100)
K = seq(0, 5, length.out = 100)
pdf(file = './plots/scatter_combined.pdf',
    width = 7,
    height = 3.2)
par(
  mfrow = c(1, 2),
  mar = c(3.5, 3.5, 0.5, 2),
  oma = c(0, 0, 0, 0)
)
plot(
  r_bh,
  r_hs,
  cex = 0.3,
  pch = 19,
  xlim = c(1, 4),
  ylim = c(1, 4),
  xlab = "",
  ylab = "",
  xaxs = 'i',
  yaxs = 'i',
  cex.axis = ca,
  xaxt = 'n',
  yaxt = 'n',
)
lines(R, bh_vs_hs_r(R),
      lty = 1, lwd = lw * 1.2)
lines(R,
      R,
      lty = 2,
      lwd = lw * 1.5,
      col = 'red')
mtext(
  side = 1,
  "Beverton-Holt, r",
  outer = FALSE,
  padj = 3,
  cex = cl
)
mtext(
  side = 2,
  "Hockey-Stick, r",
  outer = FALSE,
  padj = -3,
  cex = cl
)
text(1.9, 3.8, 'a) Growth rate, r')
axis(
  side = 2,
  at = c(0, 1, 2, 3, 4),
  labels = c(0, 1, 2, 3, 4),
  cex.axis = ca
)
axis(
  side = 1,
  at = c(0, 1, 2, 3, 4),
  labels = c(0, 1, 2, 3, 4),
  cex.axis = ca
)
plot(
  k_bh,
  k_hs,
  cex = 0.3,
  pch = 19,
  xlim = c(0, 5.1),
  ylim = c(0, 1.2),
  xlab = "",
  ylab = "",
  xaxs = 'i',
  yaxs = 'i',
  cex.axis = ca,
  xaxt = 'n',
  yaxt = 'n'
)
lines(K,
      K,
      lty = 2,
      lwd = lw * 1.5,
      col = 'red')
mtext(
  side = 1,
  "Beverton-Holt, k",
  outer = FALSE,
  padj = 3,
  cex = cl
)
mtext(
  side = 2,
  "Hockey-Stick, k",
  outer = FALSE,
  padj = -3,
  cex = cl
)
text(2, 1.12, 'b) Carrying capacity, k')

axis(
  side = 2,
  at = c(0, 0.25, 0.5, 0.75, 1),
  labels = c(0, 0.25, 0.5, 0.75, 1),
  cex.axis = ca
)
axis(
  side = 1,
  at = c(0, 1, 2, 3, 4, 5),
  labels = c(0, 1, 2, 3, 4, 5),
  cex.axis = ca
)
dev.off()
#### Fig 6 - Optimal escapement histograms ####

pdf(file = './plots/hist_combined.pdf',
    width = 7,
    height = 3.5)
par(
  mfrow = c(4, 2),
  mar = c(0, 0, 0.7, 2),
  oma = c(3.5, 3.7, 0, 2)
)
hist(
  opt_esc_bh,
  yaxs = "i",
  yaxt = 'n',
  xlim = c(0, 1),
  ylim = c(0, 125),
  breaks = seq(from = -100, to = 100, by = 0.025),
  main = '',
  xlab = '',
  ylab = '',
  cex.lab = 1.5,
  xaxt = 'n',
  col = bh_col_tr,
  cex.axis = 2
)
hist(
  opt_esc_bh_bf,
  yaxs = "i",
  yaxt = 'n',
  xlim = c(0, 1),
  ylim = c(0, 125),
  breaks = seq(from = -100, to = 100, by = 0.025),
  main = '',
  xlab = '',
  ylab = '',
  cex.lab = 1.5,
  xaxt = 'n',
  col = bh_col,
  cex.axis = 2,
  add = TRUE
)
text(0.25, 100, 'a) Proportion of fitted k', cex = 1.2)
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
  opt_esc_bh_scaled,
  yaxs = "i",
  yaxt = 'n',
  xlim = c(0, 2.5),
  ylim = c(0, 100),
  breaks = seq(from = -100, to = 100, by = 0.0625),
  main = '',
  xlab = '',
  ylab = '',
  cex.lab = 1.5,
  xaxt = 'n',
  col = bh_col_tr,
  cex.axis = 2
)
hist(
  opt_esc_bh_scaled_bf,
  yaxs = "i",
  yaxt = 'n',
  xlim = c(0, 2.5),
  ylim = c(0, 100),
  breaks = seq(from = -100, to = 100, by = 0.0625),
  main = '',
  xlab = '',
  ylab = '',
  cex.lab = 1.5,
  xaxt = 'n',
  col = bh_col,
  cex.axis = 2,
  add = TRUE
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
text(0.7, 80, 'b) Proportion of max pop', cex = 1.2)
hist(
  opt_esc_ri,
  yaxs = "i",
  yaxt = 'n',
  xlim = c(0, 1),
  ylim = c(0, 250),
  breaks = seq(from = -100, to = 100, by = 0.025),
  main = '',
  xlab = '',
  ylab = '',
  cex.lab = 1.5,
  xaxt = 'n',
  col = ri_col_tr,
  cex.axis = 2
)
hist(
  opt_esc_ri_bf,
  yaxs = "i",
  yaxt = 'n',
  xlim = c(0, 1),
  ylim = c(0, 250),
  breaks = seq(from = -100, to = 100, by = 0.025),
  main = '',
  xlab = '',
  ylab = '',
  cex.lab = 1.5,
  xaxt = 'n',
  col = ri_col,
  cex.axis = 2,
  add = TRUE
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
  opt_esc_ri_scaled,
  yaxs = "i",
  yaxt = 'n',
  xlim = c(0, 2.5),
  ylim = c(0, 100),
  breaks = seq(from = -100, to = 100, by = 0.0625),
  main = '',
  xlab = '',
  ylab = '',
  cex.lab = 1.5,
  xaxt = 'n',
  col = ri_col_tr,
  cex.axis = 2
)
hist(
  opt_esc_ri_scaled_bf,
  yaxs = "i",
  yaxt = 'n',
  xlim = c(0, 2.5),
  ylim = c(0, 100),
  breaks = seq(from = -100, to = 100, by = 0.0625),
  main = '',
  xlab = '',
  ylab = '',
  cex.lab = 1.5,
  xaxt = 'n',
  col = ri_col,
  cex.axis = 2,
  add = TRUE
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
  opt_esc_hs,
  yaxs = "i",
  yaxt = 'n',
  xlim = c(0, 1),
  ylim = c(0, 45),
  breaks = seq(from = -100, to = 100, by = 0.025),
  main = '',
  xlab = '',
  ylab = '',
  xaxt = 'n',
  col = hs_col_tr,
  cex.axis = 2
)
hist(
  opt_esc_hs_bf,
  yaxs = "i",
  yaxt = 'n',
  xlim = c(0, 1),
  ylim = c(0, 450),
  breaks = seq(from = -100, to = 100, by = 0.025),
  main = '',
  xlab = '',
  ylab = '',
  xaxt = 'n',
  col = hs_col,
  cex.axis = 2,
  add = TRUE
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
  opt_esc_hs_scaled,
  yaxs = "i",
  yaxt = 'n',
  xlim = c(0, 2.5),
  ylim = c(0, 100),
  breaks = seq(from = -100, to = 100, by = 0.0625),
  main = '',
  xlab = '',
  ylab = '',
  xaxt = 'n',
  col = hs_col_tr,
  cex.axis = 2
)
hist(
  opt_esc_hs_scaled_bf,
  yaxs = "i",
  yaxt = 'n',
  xlim = c(0, 2.5),
  ylim = c(0, 100),
  breaks = seq(from = -100, to = 100, by = 0.0625),
  main = '',
  xlab = '',
  ylab = '',
  xaxt = 'n',
  col = hs_col,
  cex.axis = 2,
  add = TRUE
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
  opt_esc_pt,
  yaxs = "i",
  yaxt = 'n',
  xlim = c(0, 1),
  ylim = c(0, 700),
  breaks = seq(from = -100, to = 100, by = 0.025),
  main = '',
  xlab = '',
  ylab = '',
  xaxt = 'n',
  col = pt_col_tr,
  cex.axis = 2
)
hist(
  opt_esc_pt_bf,
  yaxs = "i",
  yaxt = 'n',
  xlim = c(0, 1),
  ylim = c(0, 700),
  breaks = seq(from = -100, to = 100, by = 0.025),
  main = '',
  xlab = '',
  ylab = '',
  xaxt = 'n',
  col = pt_col,
  cex.axis = 2,
  add = TRUE
)
axis(
  side = 2,
  at = c(0, 300),
  labels = c(0, 300),
  cex.axis = 1.25,
  mgp = c(2, 0.5, 0)
)
axis(1,
     pos = 0,
     cex.axis = 1.25,
     mgp = c(2, 0.5, 0))
hist(
  opt_esc_pt_scaled,
  yaxs = "i",
  yaxt = 'n',
  xlim = c(0, 2.5),
  ylim = c(0, 100),
  breaks = seq(from = -100, to = 100, by = 0.0625),
  main = '',
  xlab = '',
  ylab = '',
  xaxt = 'n',
  col = pt_col_tr,
  cex.axis = 2
)
hist(
  opt_esc_pt_scaled_bf,
  yaxs = "i",
  yaxt = 'n',
  xlim = c(0, 2.5),
  ylim = c(0, 100),
  breaks = seq(from = -100, to = 100, by = 0.0625),
  main = '',
  xlab = '',
  ylab = '',
  xaxt = 'n',
  col = pt_col,
  cex.axis = 2,
  add = TRUE
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
  "Pella-Tomlinson",
  outer = TRUE,
  padj = -3,
  adj = -1.3,
  at = 0.5,
  col = pt_col
)
mtext(
  side = 1,
  "Hockey-Stick",
  outer = TRUE,
  padj = -8.5,
  adj = -1.6,
  at = 0.5,
  col = hs_col
)
mtext(
  side = 1,
  "Ricker",
  outer = TRUE,
  padj = -15,
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
