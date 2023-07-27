#### Reading in files ####
load("./data_files/RAMCore[asmt][v4.495].rdata") #RAM DATABASE
#### All data - all_cols ####
all_cols <- c()
for (i in 1:length(Bio[1, ])) {
  are_there_only_nas <- c()
  for (j in 1:length(Bio[, i])) {
    na <- is.na(Bio[, i][j]) | is.na(MCatch[, i][j])
    are_there_only_nas <- append(are_there_only_nas , na)
  }
  if (FALSE %in% are_there_only_nas) {
    all_cols <- append(all_cols, i)
  }
}

#### Non-deterministic data - non_det_cols ####
#this was visually identified manually recorded in the below csv file
non_det_cols = read.csv("./data_files/non_det_cols.csv", header = FALSE)[, 1]
#### Filtered - filtered_cols ####
#less than 50% of net production values were negative in the two middle
#quadrants between 0 and carrying capacity, K
#OR
#the sum of net production values was positive in the two middle quadrants
#between 0 and K.
num_rows <- length(Bio[, 1]) #num rows to iterate through all data
num_cols <- length(all_cols) #num cols to iterate through all data
bio_mat <- matrix(NA, num_rows, num_cols) #will store biomass(MT)
catch_mat <- matrix(NA, num_rows, num_cols) #will store catch(MT)
filtered_cols <- c() #initialising

#Creates a matrix of Biomass data and a matrix of Catch data with only
#overlapping entries.
for (i in all_cols) {
  for (j in 1:num_rows) {
    if (FALSE %in% is.na(Bio[, i][j]) &
        FALSE %in% is.na(MCatch[, i][j])) {
      bio_mat[j, match(i, all_cols)] <- c(Bio[, i][j])
      catch_mat[j, match(i, all_cols)] <- c(MCatch[, i][j])
    }
  }
}

#Filters
for (i in 1:num_cols) {
  num_pos <- 0
  num_neg <- 0
  ps <- c()
  
  bios <- bio_mat[, i][!is.na(bio_mat[, i])]
  catchs <- catch_mat[, i][!is.na(catch_mat[, i])]
  
  X <- bios[1:(length(bios) - 1)] - catchs[1:(length(bios) - 1)]
  Y <- bios[2:length(bios)]
  P <-
    bios[2:length(bios)] - bios[1:(length(bios) - 1)] + catchs[2:(length(bios))]
  
  q1 = quantile(Y, prob = c(.25, .5, .75), type = 1)[1]
  q2 = quantile(Y, prob = c(.25, .5, .75), type = 1)[2]
  q3 = quantile(Y, prob = c(.25, .5, .75), type = 1)[3]
  
  for (j in 1:length(P)) {
    if (q1 < bios[j + 1] & bios[j + 1] < q3) {
      ps <- append(ps, P[j])
    }
    for (k in ps) {
      if (k > 0) {
        num_pos <- num_pos + 1
      }
      else if (k < 0) {
        num_neg <- num_neg + 1
      }
    }
  }
  if (num_pos > num_neg | sum(ps) > 0) {
    filtered_cols <- append(filtered_cols, all_cols[i])
  }
}

#### Constants ####
years = seq(1950, 2020, 1) #year of each row of data
num_rows = length(years) #num rows to iterate through all data
k_max = 5 #upper bound on carrying capacity (multiple of max observed population)
r_max = 100 #upper bound on growth rate
phi = 0.188
col_names = colnames(Bio)
bh_col = 'darkblue'
bh_col_tr = '#4281C388'
ri_col = 'red'
ri_col_tr = '#E4231388'
hs_col = 'darkgreen'
hs_col_tr = '#00843788'
pt_col = '#DB7093'
pt_col_tr = '#DB709388'
lw = 2
ca = 1.2
cl = 1.3

#### Functions ####

#Beverton-Holt
BH = function(r, k = 1, x) {
  return(r * x / (1 + (r - 1) * x / k))
}

#Hockey-Stick
HS = function(r, k = 1, x) {
  return(apply(cbind(r * x, k), 1, min))
}

#Ricker
RI = function(r, k = 1, x) {
  return(x * r ^ (1 - x / k))
}

#Pella-Tomlinson
PT = function(r, k = 1, x , phi = 0.188) {
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
pt_opt_esc = ((1 / (1 + phi)) ^ (1 / phi))


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

#### Main ####
for(which_data_set in list(all_cols, non_det_cols, filtered_cols)) {
  num_cols = length(which_data_set)
  # Initialising vectors and constants
  bio_mat <-
    matrix(NA, num_rows, num_cols) #biomass(MT)
  catch_mat <-
    matrix(NA, num_rows, num_cols) #catch(MT)
  
  r_bh <-
    rep(NA, num_cols) #estimated growth rates for Beverton-Holt
  r_hs <- rep(NA, num_cols) #estimated growth rates for Hockey-Stick
  r_ri <- rep(NA, num_cols) #estimated growth rates for Ricker
  r_pt <-
    rep(NA, num_cols) #estimated growth rates for Pella-Tomlinson
  
  k_bh <-
    rep(NA, num_cols) #estimated carrying capacity/max pop for Beverton-Holt
  k_hs <-
    rep(NA, num_cols) #estimated carrying capacity/max pop for Hockey-Stick
  k_ri <-
    rep(NA, num_cols) #estimated carrying capacity/max pop Ricker
  k_pt <-
    rep(NA, num_cols) #estimated carrying capacity for Pella-Tomlinson
  
  r2_bh <- rep(NA, num_cols) #r-sqaured for Beverton-Holt
  r2_hs <- rep(NA, num_cols) #r-squared for Hockey-Stick
  r2_ri <- rep(NA, num_cols) #r-squared for Ricker
  r2_pt <- rep(NA, num_cols) #r-squared for Pella-Tomlinson
  
  num_bh <-
    0 #number of data sets where Beverton-Holt is best fit
  num_hs <-
    0 #number of data sets where Hockey-Stick is best fit
  num_ri <-
    0 #number of data sets where Ricker is best fit
  num_pt <-
    0 #number of data sets where Pella-Tomlinson is best fit
  
  opt_esc_bh <-
    rep(NA, num_cols) #estimated optimal escapement for Beverton-Holt (k = 1)
  opt_esc_hs <-
    rep(NA, num_cols) #estimated optimal escapement for Hockey-Stick (k = 1)
  opt_esc_ri <-
    rep(NA, num_cols) #estimated optimal escapement for Ricker (k = 1)
  opt_esc_pt <-
    rep(NA, num_cols) #estimated optimal escapement for Pella-Tomlinson (k = 1)
  
  opt_esc_bh_scaled <-
    rep(NA, num_cols) #estimated optimal escapement for Beverton-Holt (proportion of max pop)
  opt_esc_hs_scaled <-
    rep(NA, num_cols) #estimated optimal escapement for Hockey-Stick (proportion of max pop)
  opt_esc_ri_scaled <-
    rep(NA, num_cols) #estimated optimal escapement for Ricker (proportion of max pop)
  opt_esc_pt_scaled <-
    rep(NA, num_cols) #estimated optimal escapement for Pella-Tomlinson (proportion of max pop)
  
  num_hs_lower_bh_r <- 0 #number of times Hockey-Stick suggests a lower growth rate than Beverton-Holt
  num_hs_lower_bh_k <- 0 #number of times Hockey-Stick suggests a lower carrying capacity than Beverton-Holt
  num_hs_higher_bh_opt_esc_prop_k <- 0 #number of times Hockey-Stick suggests a higher optimal escapement than Beverton-Holt (as a proportion of the fitted carrying capacity)
  num_hs_higher_bh_opt_esc_prop_max_pop <- 0 #number of times Hockey-Stick suggests a higher optimal escapement than Beverton-Holt (as a proportion of the max population)
  num_hs_higher_b60_opt_esc_prop_k <- 0 #number of times Hockey-Stick suggests an escapement higher than 60% of its fitted carrying capacity
  num_bh_k_higher_1 <- 0 #number of times Beverton-Holt suggests a carrying capacity higher than 1
  num_hs_higher_bh_opt_esc_replacing_1s <- 0 #number of times Hockey-Stick suggests a higher optimal escapement than Beverton-Holt if we take the Beverton-Holt estimate for carrying capacity to be true when Hockey-Stick suggests a carrying capacity equal to 1. 
  num_neg = 0 #number of data sets with negative values for biomass
  
  #Sorting data
  #Creates a matrix of Biomass data and a matrix of Catch data with only
  #overlapping entries.
  for (i in which_data_set) {
    for (j in 1:num_rows) {
      if (FALSE %in% is.na(Bio[, i][j]) &
          FALSE %in% is.na(MCatch[, i][j])) {
        bio_mat[j, match(i, which_data_set)] <- c(Bio[, i][j])
        catch_mat[j, match(i, which_data_set)] <-
          c(MCatch[, i][j])
      }
    }
  }
  
  for (i in  1:num_cols) {
    bios <- bio_mat[, i][!is.na(bio_mat[, i])]
    catchs <- catch_mat[, i][!is.na(catch_mat[, i])]
    X <-
      bios[1:(length(bios) - 1)] - catchs[1:(length(bios) - 1)]
    Y <- bios[2:length(bios)]
    
    #Handling data sets with negative values of biomass
    if (min(X, Y) < 0) {
      k_bh[i] <- 0
      k_hs[i] <- 0
      k_ri[i] <- 0
      k_pt[i] <- 0
      
      r_bh[i] <- 0
      r_hs[i] <- 0
      r_ri[i] <- 0
      r_pt[i] <- 0
      
      r2_bh[i] <- 0
      r2_hs[i] <- 0
      r2_ri[i] <- 0
      r2_pt[i] <- 0
      
      num_neg = num_neg + 1
      
      #stores opitmal escapement as a proportion of the fitted K
      opt_esc_bh[i] <- 0
      opt_esc_hs[i] <- 0
      opt_esc_ri[i] <- 0
      opt_esc_pt[i] <- 0
      
      #stores optimal escapement as a proportion of the max population
      opt_esc_bh_scaled[i] <- 0
      opt_esc_hs_scaled[i] <- 0
      opt_esc_ri_scaled[i] <- 0
      opt_esc_pt_scaled[i] <- 0
    } else{
      #intial fact_guess before optimisation
      init_guess <- c(mean(Y / X), max(X, Y)) #init_guess
      fact_guess <-
        c(0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75) #taking factors of init_guess for robustness
      num_guess <- length(fact_guess) #size of grid of init guesses
      
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
      
      #sorting through grod of initial guesses to find the parameter estimates with the best fit for each model
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
          
          
          if (hs_opt[[1]][1] > hs_opt[[1]][2] / min(X)) {
            hs_opt[[1]][1] <- hs_opt[[1]][2] / min(X)
          }
          if (hs_opt[[1]][2] > max(X, Y)) {
            hs_opt[[1]][2] <- max(X, Y)
          }
          
          #estimated parameters stored
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
      bh_opt_r <- bh_r_guess_mat[bh_jk[1],][bh_jk[2]]
      bh_opt_k <- bh_k_guess_mat[bh_jk[1],][bh_jk[2]]
      k_bh[i] <- bh_opt_k / max(Y)
      r_bh[i] <- bh_opt_r
      r2_bh[i] <- max(bh_r2_guess_mat)
      
      
      
      hs_jk <-
        which(hs_r2_guess_mat == max(hs_r2_guess_mat), arr.ind <-
                TRUE)[1, ]
      hs_opt_r <- hs_r_guess_mat[hs_jk[1],][hs_jk[2]]
      hs_opt_k <- hs_k_guess_mat[hs_jk[1],][hs_jk[2]]
      k_hs[i] <- hs_opt_k / max(Y)
      r_hs[i] <- hs_opt_r
      r2_hs[i] <- max(hs_r2_guess_mat)
      
      ri_jk <-
        which(ri_r2_guess_mat == max(ri_r2_guess_mat), arr.ind <-
                TRUE)[1, ]
      ri_opt_r <- ri_r_guess_mat[ri_jk[1],][ri_jk[2]]
      ri_opt_k <- ri_k_guess_mat[ri_jk[1],][ri_jk[2]]
      k_ri[i] <- ri_opt_k / max(Y)
      r_ri[i] <- ri_opt_r
      r2_ri[i] <- max(ri_r2_guess_mat)
      
      pt_jk <-
        which(pt_r2_guess_mat == max(pt_r2_guess_mat), arr.ind <-
                TRUE)[1, ]
      pt_opt_r <- pt_r_guess_mat[pt_jk[1],][pt_jk[2]]
      pt_opt_k <- pt_k_guess_mat[pt_jk[1],][pt_jk[2]]
      k_pt[i] <- pt_opt_k / max(Y)
      r_pt[i] <- pt_opt_r
      r2_pt[i] <- max(pt_r2_guess_mat)
      
      #stores opitmal escapement as a proportion of the fitted K
      opt_esc_bh[i] <- bh_opt_esc(r_bh[i], 1)
      opt_esc_hs[i] <- hs_opt_esc(r_hs[i], 1)
      p = r_ri[i]
      k = 1
      result = uniroot(derivative_ri, lower = 0, upper = 10 * k)
      opt_esc_ri[i] <- result$root
      opt_esc_pt[i] <- (1 / (1 + phi)) ^ (1 / phi)
      
      #stores optimal escapement as a proportion of the max population
      opt_esc_bh_scaled[i] <- bh_opt_esc(r_bh[i], k_bh[i])
      opt_esc_hs_scaled[i] <- hs_opt_esc(r_hs[i], k_hs[i])
      p = r_ri[i]
      k = k_ri[i]
      result = uniroot(derivative_ri, lower = 0, upper = 10 * k)
      opt_esc_ri_scaled[i] <- result$root
      opt_esc_pt_scaled[i] <- k_pt[i] * (1 / (1 + phi)) ^
        (1 / phi)
      
      
      #counts the number of times each model is best-fit
      if (r2_bh[i] == max(c(r2_bh[i], r2_ri[i], r2_hs[i], r2_pt[i]))) {
        num_bh = num_bh + 1
      }
      if (r2_ri[i] == max(c(r2_bh[i], r2_ri[i], r2_hs[i], r2_pt[i]))) {
        num_ri = num_ri + 1
      }
      if (r2_hs[i] == max(c(r2_bh[i], r2_ri[i], r2_hs[i], r2_pt[i]))) {
        num_hs = num_hs + 1
      }
      if (r2_pt[i] == max(c(r2_bh[i], r2_ri[i], r2_hs[i], r2_pt[i]))) {
        num_pt = num_pt + 1
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
  }
  
  #### Results ####
  if(num_cols == length(all_cols)) {
    print("Results for entire data set:")
    cat("\n")
  }
  if (num_cols == length(filtered_cols)) {
    print("Results for filtered data set:")
    cat("\n")
  }
  if (num_cols == length(non_det_cols)) {
    print("Results for non-deterministic data set:")
    cat("\n")
  }
  print(paste0("Number of data sets = ", num_cols - num_neg))
  cat("\n")
  print("Best fit as a percentage.")
  print(paste0(
    "BH: ",
    round(100 * num_bh / (num_cols - num_neg)),
    "%, HS: ",
    round(100 * num_hs / (num_cols - num_neg)),
    "%, RI: ",
    round(100 * num_ri / (num_cols - num_neg)),
    "%, PT: ",
    round(100 * num_pt / (num_cols - num_neg)),
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
    round(100 * num_hs_lower_bh_r / (num_cols - num_neg)),
    "% of the time."
  ))
  cat("\n")
  print(paste0(
    "HS suggests a carrying capacity lower than BH ",
    round(100 * num_hs_lower_bh_k / (num_cols - num_neg)),
    "% of the time."
  ))
  cat("\n")
  print(
    paste0(
      "HS suggests higher optimal escapement than BH as a proportion of k ",
      round(100 * num_hs_higher_bh_opt_esc_prop_k / (num_cols - num_neg)),
      "% of the time."
    )
  )
  cat("\n")
  print(paste0(
    "HS suggests escapement higher than B60 ",
    round(100 * num_hs_higher_b60_opt_esc_prop_k / (num_cols - num_neg)) ,
    "% of the time."
  ))
  cat("\n")
  print(
    paste0(
      "HS suggests higher optimal escapement than BH as a proportion of the max pop ",
      round(
        100 * num_hs_higher_bh_opt_esc_prop_max_pop / (num_cols - num_neg)
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
        ) / (num_cols - num_neg)
      ),
      "% if we take the BH k as true when HS k = 1."
    )
  )
  cat("\n")
  print(
    paste0(
      "BH suggests a carrying capacity greater than the max pop ",
      round(100 * num_bh_k_higher_1 / (num_cols - num_neg)),
      "% of the time."
    )
  )
  cat("\n")
  cat("\n")
}

#### Printed results ####
###############################################################################
# "Results for entire data set:"
# 
# "Number of data sets = 476"
# 
# "Best fit as a percentage."
# "BH: 12%, HS: 31%, RI: 15%, PT: 43%"
# 
# "Median r-squared values."
# "BH: 0.91, HS: 0.9, RI: 0.91, PT: 0.91"
# 
# "HS suggests a growth rate lower than BH 99% of the time,"
# 
# "HS suggests a carrying capacity lower than BH 71% of the time."
# 
# "HS suggests higher optimal escapement than BH as a proportion of k 97% of the time."
# 
# "HS suggests escapement higher than B60 87% of the time."
# 
# "HS suggests higher optimal escapement than BH as a proportion of the max pop 76% of the time"
# 
# "This increases to 82% if we take the BH k as true when HS k = 1."
# 
# "BH suggests a carrying capacity greater than the max pop 53% of the time."
# 
###############################################################################
# "Results for non-deterministic data set:"
# 
# "Number of data sets = 284"
# 
# "Best fit as a percentage."
# "BH: 15%, HS: 35%, RI: 14%, PT: 37%"
# 
# "Median r-squared values."
# "BH: 0.83, HS: 0.81, RI: 0.83, PT: 0.82"
# 
# "HS suggests a growth rate lower than BH 99% of the time."
# 
# "HS suggests a carrying capacity lower than BH 76% of the time."
# 
# "HS suggests higher optimal escapement than BH as a proportion of k 97% of the time."
# 
# "HS suggests escapement higher than B60 82% of the time."
# 
# "HS suggests higher optimal escapement than BH as a proportion of the max pop 74% of the time."
# 
# "This increases to 79% if we take the BH k as true when HS k = 1."
# 
# "BH suggests a carrying capacity greater than the max pop 52% of the time."
# 
###############################################################################
# "Results for filtered data set:"
# 
# "Number of data sets = 452"
# 
# "Best fit as a percentage."
# "BH: 12%, HS: 30%, RI: 14%, PT: 43%"
# 
# "Median r-squared values."
# "BH: 0.91, HS: 0.9, RI: 0.91, PT: 0.91"
# 
# "HS suggests a growth rate lower than BH 100% of the time."
# 
# "HS suggests a carrying capacity lower than BH 75% of the time."
# 
# "HS suggests higher optimal escapement than BH as a proportion of k 98% of the time."
#
# "HS suggests escapement higher than B60 86% of the time."
# 
# "HS suggests higher optimal escapement than BH as a proportion of the max pop 75% of the time."
# 
# "This increases to 81% if we take the BH k as true when HS k = 1."
# 
# "BH suggests a carrying capacity greater than the max pop 56% of the time."
