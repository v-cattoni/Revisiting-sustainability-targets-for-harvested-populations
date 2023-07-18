#### Reading in files #### 
load("./data_files/RAMCore[asmt][v4.495].rdata") #RAM DATABASE 

#### All data - all_cols ####
all_cols <- c()
for(i in 1:length(Bio[1,])){
  are_there_only_nas <- c()
  for(j in 1:length(Bio[,i])){
    na <- is.na(Bio[,i][j]) | is.na(MCatch[,i][j])
    are_there_only_nas <- append(are_there_only_nas ,na)
  }
  if(FALSE %in% are_there_only_nas){
    all_cols <- append(all_cols, i)
  }
}

#### Non-deterministic data - non_det_cols ####
#this was visually identified manually recorded in the below csv file
non_det_cols = read.csv("./data_files/non_det_cols.csv", header = FALSE)[,1]

#### Costello filters- costello_cols ####
#less than 50% of net production values were negative in the two middle 
#quadrants between 0 and carrying capacity, K
#OR
#the sum of net production values was positive in the two middle quadrants 
#between 0 and K.
num_rows <- length(Bio[,1]) #num rows to iterate through all data
num_cols <- length(all_cols) #num cols to iterate through all data
bio_mat <- matrix(NA, num_rows, num_cols) #will store biomass(MT)
catch_mat <- matrix(NA, num_rows, num_cols) #will store catch(MT)
costello_cols <- c() #initialising

#Creates a matrix of Biomass data and a matrix of Catch data with only
#overlapping entries.
for(i in all_cols){
  for(j in 1:num_rows){
    if(FALSE %in% is.na(Bio[,i][j]) & FALSE %in% is.na(MCatch[,i][j])){
      bio_mat[j, match(i, all_cols)] <- c(Bio[,i][j])
      catch_mat[j, match(i, all_cols)] <- c(MCatch[,i][j])
    }
  }
}

#costello filters
for(i in 1:num_cols){
  num_pos <- 0
  num_neg <- 0
  ps <- c()
  
  bios <- bio_mat[, i][!is.na(bio_mat[,i])]
  catchs <- catch_mat[, i][!is.na(catch_mat[,i])]
  
  X <- bios[1:(length(bios)-1)] - catchs[1:(length(bios)-1)]
  Y <- bios[2:length(bios)]
  P <- bios[2:length(bios)] - bios[1:(length(bios)-1)] + catchs[2:(length(bios))]
  
  q1 = quantile(Y, prob=c(.25,.5,.75), type=1)[1]
  q2 = quantile(Y, prob=c(.25,.5,.75), type=1)[2]
  q3 = quantile(Y, prob=c(.25,.5,.75), type=1)[3]
  
  for(j in 1:length(P)){
    if(q1 < bios[j+1] & bios[j+1] < q3){
      ps <- append(ps, P[j])
    }
    for(k in ps){
      if(k > 0){
        num_pos <- num_pos + 1
      }
      else if(k < 0){
        num_neg <- num_neg + 1
      }
    }
  }
  if(num_pos > num_neg | sum(ps) > 0){
    costello_cols <- append(costello_cols, all_cols[i])
  }
}

#### Selecting a data set ####
#which data set would you like to examine? 1, 2, or 3
which_data_set <- 2

if(which_data_set == 1){
  cols = all_cols
} else if(which_data_set == 2){
  cols = non_det_cols
} else if(which_data_set == 3){
  cols = costello_cols
} else{
  stop('Please select a valid data set on line 81.')
}

#### Constants ####
years = seq(1950, 2020, 1) #year of each row of data
num_rows = length(years) #num rows to iterate through all data
num_cols = length(cols) #num cols to iterate through all data
k_max = 5 #upper bound on carrying capacity (multiple of max observed population)
r_max = 100 #upper bound on growth rate
phi = 0.188
col_names = colnames(Bio)
bh_col = 'darkblue'
bh_col_tr = '#4281C388'
rck_col = 'red'
rck_col_tr = '#E4231388'
hs_col = 'darkgreen'
hs_col_tr = '#00843788'
pt_col = '#DB7093'
pt_col_tr = '#DB709388'
lw = 2; ca=1.2; cl = 1.3;

#### Initialising vectors and constants ####
bio_mat <- matrix(NA, num_rows, num_cols) #biomass(MT)
catch_mat <- matrix(NA, num_rows, num_cols) #catch(MT)

r_bh <- c() #estimated growth rates for Beverton-Holt
r_hs <- c() #estimated growth rates for Hockey-Stick
r_rck <- c() #estimated growth rates for Ricker
r_pt <- c() #estimated growth rates for Pella-Tomlinson

k_bh <- c() #estimated carrying capacity/max pop for Beverton-Holt
k_hs <- c() #estimated carrying capacity/max pop for Hockey-Stick
k_rck <- c() #estimated carrying capacity/max pop Ricker
k_pt <- c() #estimated carrying capacity for Pella-Tomlinson

r2_bh <- c() #r-sqaured for Beverton-Holt
r2_hs <- c() #r-squared for Hockey-Stick
r2_rck <- c() #r-squared for Ricker
r2_pt <- c() #r-squared for Pella-Tomlinson

l <- c() #number of data points in each data set

num_bh <- 0 #number of data sets where Beverton-Holt is best fit
num_hs <- 0 #number of data sets where Hockey-Stick is best fit
num_rck <- 0 #number of data sets where Ricker is best fit
num_pt <- 0 #number of data sets where Pella-Tomlinson is best fit

num_bh_rs <-c() #list containing number of times Beverton-Holt is best for with increasing upper bound on r (used for cumulative frequency plot)
num_hs_rs <- c() #list containing number of times Hockey-Stick is best for with increasing upper bound on r (used for cumulative frequency plot)
num_rck_rs <- c() #list containing number of times Ricker is best for with increasing upper bound on r (used for cumulative frequency plot)
num_pt_rs <- c() #list containing number of times Pella-Tomlinson is best for with increasing upper bound on r (used for cumulative frequency plot)

opt_esc_bh <- c() #estimated optimal escapement for Beverton-Holt (k = 1)
opt_esc_hs <- c() #estimated optimal escapement for Hockey-Stick (k = 1)
opt_esc_rck <- c() #estimated optimal escapement for Ricker (k = 1)
opt_esc_pt <- c() #estimated optimal escapement for Pella-Tomlinson (k = 1)

opt_esc_bh_bf <- c() #estimated optimal escapement for Beverton-Holt (k = 1) (best fit)
opt_esc_hs_bf <- c() #estimated optimal escapement for Hockey-Stick (k = 1) (best fit)
opt_esc_rck_bf <- c() #estimated optimal escapement for Ricker (k = 1) (best fit)
opt_esc_pt_bf <- c() #estimated optimal escapement for Pella-Tomlinson (k = 1) (best fit)

opt_esc_bh_scaled <- c() #estimated optimal escapement for Beverton-Holt (proportion of max pop)
opt_esc_hs_scaled <- c() #estimated optimal escapement for Hockey-Stick (proportion of max pop)
opt_esc_rck_scaled <- c() #estimated optimal escapement for Ricker (proportion of max pop)
opt_esc_pt_scaled <- c() #estimated optimal escapement for Pella-Tomlinson (proportion of max pop)

opt_esc_bh_scaled_bf <- c() #estimated optimal escapement for Beverton-Holt (proportion of max pop) (best fit)
opt_esc_hs_scaled_bf <- c() #estimated optimal escapement for Hockey-Stick (proportion of max pop) (best fit)
opt_esc_rck_scaled_bf <- c() #estimated optimal escapement for Ricker (proportion of max pop) (best fit)
opt_esc_pt_scaled_bf <- c() #estimated optimal escapement for Pella-Tomlinson (proportion of max pop) (best fit)
opt_esc_data <- c() #records optimal escapement of the data

#### Functions ####

#Beverton-Holt
BH = function(r, k=1, x){
  return(r*x/(1+(r-1)*x/k))
}

#Hockey-Stick
HS = function(r, k=1, x){
  return(apply(cbind(r*x,k),1,min))
}

#Ricker
RI = function(r, k=1, x){
  return(x*r^(1-x/k))
}

#Pella-Tomlinson
PT = function(r, k=1, x ,phi=0.188){
  return(x+x*(r-1)*(1-(x/k)^phi))
}


#Beverton-Holt optimal escapement
bh_opt_esc <- function(r, k){
  return(k*(r**0.5 - 1)/(r-1))
}

#Hockey-Stick optimal escpaement
hs_opt_esc <- function(r, k){
  return(k/r)
}

#Ricker optimal escapement
# Define the function to find the root of
derivative_rck <- function(X) {
  # Calculate the derivative of rck with respect to X
  derivative <- (p^(1-(X/k)))*(k-X*log(p))/k -1
  
  # Return the derivative
  return(derivative)
}

#Pella-Tomlinson optimal escapement 
pt_opt_esc = ((1/(1+phi))^(1/phi))


#Sum of Squares Error
sse <- function(X, Y){
  s <- 0
  for(i in 1:length(X)){
    s <- s + (X[i] - Y[i])**2
  }
  return(s)
}

#Sum of Squares Total
sst <- function(Y){
  s <- 0
  ym <- mean(Y)
  for(i in 1:length(Y)){
    s <- s + (Y[i] - ym)**2
  }
  return(s)
}

#R squared
r_squared <- function(X, Y){
  return(1 - sse(X, Y)/sst(Y))
}

#Sum of Squares Error (Beverton-Holt) 
#This is the function which will optimised in terms of r and k
sse_bh <- function(p, X, Y){
  r <- p[1]
  k <- p[2]
  return(sse(Y, BH(r, k, X)))
}

#Sum of Squares Error (Hockey-Stick)
#This is the function which will optimised in terms of r and k
sse_hs <- function(p, X, Y){
  r <- p[1]
  k <- p[2]
  return(sse(Y, HS(r, k, X)))
}

#Sum of Squares Error (Ricker)
#This is the function which will optimised in terms of r and k
sse_rck <- function(p, X, Y){
  r <- p[1]
  k <- p[2]
  return(sse(Y, RI(r, k, X)))
}

#Sum of Squares Error (Pella-Tomlinson)
#This is the function which will optimised in terms of r and k
sse_pt <- function(p, X, Y){
  r <- p[1]
  k <- p[2]
  return(sse(Y, PT(r, k, X)))
}

#### Sorting data ####
#Creates a matrix of Biomass data and a matrix of Catch data with only
#overlapping entries.
for(i in cols){
  for(j in 1:num_rows){
    if(FALSE %in% is.na(Bio[,i][j]) & FALSE %in% is.na(MCatch[,i][j])){
      bio_mat[j, match(i, cols)] <- c(Bio[,i][j])
      catch_mat[j, match(i, cols)] <- c(MCatch[,i][j])
    }
  }
}

#### Main ####
for(i in 1:length(cols)){
  bios <- bio_mat[, i][!is.na(bio_mat[,i])]
  catchs <- catch_mat[, i][!is.na(catch_mat[,i])]
  X <- bios[1:(length(bios)-1)] - catchs[1:(length(bios)-1)]
  Y <- bios[2:length(bios)]

  #intial guesses before optimisation
  if(min(X) > 0){ 
    bh_p_guess <- c(mean(Y/X), max(Y))
    hs_p_guess <- c(mean(Y/X), max(Y))
    rck_p_guess <- c(mean(Y/X), max(Y))
    pt_p_guess <- c(mean(Y/X), max(Y))
  }
  if(min(X < 0)){
    bh_p_guess <- abs(c(mean(Y/X), max(Y)))
    hs_p_guess <- abs(c(mean(Y/X), max(Y)))
    rck_p_guess <- abs(c(mean(Y/X), max(Y)))
    pt_p_guess <- abs(c(mean(Y/X), max(Y)))
  }
  
  #parameters estimated 
  bh_opt <- optim(bh_p_guess, sse_bh, X = X, Y = Y, method="L-BFGS-B", lower = c(1.00001, mean(Y)), upper = c(min(r_max, max(Y)/min(X)), max(Y)*k_max), control = c(factr = 1e-15, factr = 1e-15*max(Y)/mean(Y/X)))
  hs_opt <- optim(hs_p_guess, sse_hs, X = X, Y = Y, method="L-BFGS-B", lower = c(1.00001, mean(Y)), upper = c(min(r_max, max(Y)/min(X)), max(Y)), control = c(factr = 1e-15, factr = 1e-15*max(Y)/mean(Y/X)))
  rck_opt <- optim(rck_p_guess, sse_rck, X = X, Y = Y, method="L-BFGS-B", lower = c(1.00001, mean(Y)), upper = c(min(r_max, max(Y)/min(X)), max(Y)*k_max), control = c(factr = 1e-15, factr = 1e-15*max(Y)/mean(Y/X)))
  pt_opt <- optim(pt_p_guess, sse_pt, X = X, Y = Y, method="L-BFGS-B", lower = c(1.00001, mean(Y)), upper = c(min(r_max, max(Y)/min(X)), max(Y)*k_max), control = c(factr = 1e-15, factr = 1e-15*max(Y)/mean(Y/X)))
  
  #estimated parameters stored
  bh_opt_r <- bh_opt[[1]][1]
  bh_opt_k <- bh_opt[[1]][2]
  hs_opt_r <- hs_opt[[1]][1]
  hs_opt_k <- hs_opt[[1]][2]
  rck_opt_r <- rck_opt[[1]][1]
  rck_opt_k <- rck_opt[[1]][2]
  pt_opt_r <- pt_opt[[1]][1]
  pt_opt_k <- pt_opt[[1]][2]
  
  if(hs_opt_r > hs_opt_k/min(X)){
    hs_opt_r = hs_opt_k/min(X)
  }
  
  #storing number of data points
  l <- append(l, length(X))
  
  #storing carrying capacity/max pops
  k_bh <- append(k_bh, bh_opt_k/max(Y))
  k_hs <- append(k_hs, hs_opt_k/max(Y))
  k_rck <- append(k_rck, rck_opt_k/max(Y))
  k_pt <- append(k_pt, pt_opt_k/max(Y))
  
  #storing growth rates
  r_bh <- append(r_bh, bh_opt_r)
  r_hs <- append(r_hs, hs_opt_r)
  r_rck <- append(r_rck, rck_opt_r)
  r_pt <- append(r_pt, pt_opt_r)
  
  #storing r-squared values
  r2_bh <- append(r2_bh, r_squared(BH(bh_opt_r, bh_opt_k, X), Y))
  r2_hs <- append(r2_hs, r_squared(HS(hs_opt_r, hs_opt_k, X), Y))
  r2_rck <- append(r2_rck, r_squared(RI(rck_opt_r, rck_opt_k, X), Y))
  r2_pt <- append(r2_pt, r_squared(PT(pt_opt_r, pt_opt_k, X), Y))
  
  if(is.na(r2_bh[i])){
    r2_bh[i] = 0
  }
  if(is.na(r2_hs[i])){
    r2_hs[i] = 0
  }
  if(is.na(r2_rck[i])){
    r2_rck[i] = 0
  }
  if(is.na(r2_pt[i])){
    r2_pt[i] = 0
  }
  
  #stores opitmal escapement as a proportion of the fitted K
  opt_esc_bh = append(opt_esc_bh, bh_opt_esc(r_bh[i], 1))
  opt_esc_hs = append(opt_esc_hs, hs_opt_esc(r_hs[i], 1))
  p = r_rck[i]
  k = 1
  result = uniroot(derivative_rck, lower = 0, upper = 10*k)
  opt_esc_rck = append(opt_esc_rck, result$root)
  opt_esc_pt = append(opt_esc_pt, (1/(1 + phi))^(1/phi))
  
  #stores optimal escapement as a proportion of the max population
  opt_esc_bh_scaled = append(opt_esc_bh_scaled, bh_opt_esc(r_bh[i], k_bh[i]))
  opt_esc_hs_scaled = append(opt_esc_hs_scaled, hs_opt_esc(r_hs[i], k_hs[i]))
  p = r_rck[i]
  k = k_rck[i]
  result = uniroot(derivative_rck, lower = 0, upper = 10*k)
  opt_esc_rck_scaled = append(opt_esc_rck_scaled, result$root)
  opt_esc_pt_scaled = append(opt_esc_pt_scaled, k_pt[i]*(1/(1 + phi))^(1/phi))
  
  
  #counts the number of times each model is best-fit
  if(r2_bh[i] == max(c(r2_bh[i], r2_rck[i], r2_hs[i], r2_pt[i]))){
    num_bh = num_bh + 1
    opt_esc_bh_bf <- append(opt_esc_bh_bf, opt_esc_bh[i])
    opt_esc_bh_scaled_bf <- append(opt_esc_bh_scaled_bf, opt_esc_bh_scaled[i])
  }
  if(r2_rck[i] == max(c(r2_bh[i], r2_rck[i], r2_hs[i], r2_pt[i]))){
    num_rck = num_rck + 1
    opt_esc_rck_bf <- append(opt_esc_rck_bf, opt_esc_rck[i])
    opt_esc_rck_scaled_bf <- append(opt_esc_rck_scaled_bf, opt_esc_rck_scaled[i])
  }
  if(r2_hs[i] == max(c(r2_bh[i], r2_rck[i], r2_hs[i], r2_pt[i]))){
    num_hs = num_hs + 1
    opt_esc_hs_bf <- append(opt_esc_hs_bf, opt_esc_hs[i])
    opt_esc_hs_scaled_bf <- append(opt_esc_hs_scaled_bf, opt_esc_hs_scaled[i])
  }
  if(r2_pt[i] == max(c(r2_bh[i], r2_rck[i], r2_hs[i], r2_pt[i]))){
    num_pt = num_pt + 1
    opt_esc_pt_bf <- append(opt_esc_pt_bf, opt_esc_pt[i])
    opt_esc_pt_scaled_bf <- append(opt_esc_pt_scaled_bf, opt_esc_pt_scaled[i])
  }
}
#### Figures ####
if(which_data_set == 2){ #only prints figures for non-deterministic data
#### Fig 1 - Recruitment functions ####
pdf(file = './plots/recruitment_functions.pdf', width = 3.5, height = 3.5)
par(mfrow = c(1,1), mar = c(1,1,.7,.7), oma= c(2.51,2.51,0,0)); 
r=1.7;
x=seq(0,1.6, by=.005)
plot(x, HS(r = r, x = x), ylim = c(0,1.2),
     xlab="", ylab="", lty=2,
     type = 'l', xaxs='i', yaxs='i', cex.axis = ca,
     lwd = lw*1.2, col = hs_col)
lines(x, BH(r = r, x = x), ylim = c(0,1),
      xlab="", ylab="", lty=1,
      type = 'l', xaxs='i', yaxs='i', cex.axis = ca,
      lwd = lw, col = bh_col)
lines(x, RI(r = r, x = x), ylim = c(0,1),
      xlab="", ylab="", lty=3,
      type = 'l', xaxs='i', yaxs='i', cex.axis = ca,
      lwd = lw*1.7, col = rck_col)
lines(x, PT(r = r, x = x), ylim = c(0,1),
      xlab="", ylab="", lty=4,
      type = 'l', xaxs='i', yaxs='i', cex.axis = ca,
      lwd = lw*1.4, col = pt_col)
legend('bottomright',
       c('Beverton-Holt', 'Hockey-stick','Ricker','Pella-Tomlinson'),
       col = c(bh_col,hs_col,rck_col,pt_col),
       lty = c(1,2,3,4),
       bty='n', cex=.9, lwd = c(lw,lw*1.2,lw*1.7,lw*1.5))
mtext(side=1, "Population size, year t", outer = TRUE, padj=1.9, cex = cl)
mtext(side=2, "Population size, year t+1", outer = TRUE, padj=-1.9, cex = cl)
dev.off()



#### Fig 2 - Optimal escapement functions ####
pdf(file = './plots/opt_esc_curves.pdf', width = 3.5, height = 3.5)
par(mfrow = c(1,1), mar = c(1,1,.7,.7), oma= c(2.51,2.51,0,0)); 
r=1.7;
x=seq(0,1.6, by=.005)
R=seq(1.01,4, by=.01)
e.B = (sqrt(R)-1)/(R-1)
e.H = 1/R
B60 = rep(.6, length.out=length(R))

rck_opt_xs = c()
for(i in R){
  p = i
  result <- uniroot(derivative_rck, lower = 0, upper = 10)
  rck_opt_xs = append(rck_opt_xs, result$root)
}

plot(R, 1/R, ylim = c(0,1.05),
     xlab="", ylab="", lty=2,
     type = 'l', xaxs='i', yaxs='i', cex.axis = ca,
     lwd = lw*1.2, col = hs_col)
lines(R, e.B, ylim = c(0,1),
      xlab="", ylab="", lty=1,
      type = 'l', xaxs='i', yaxs='i', cex.axis = ca,
      lwd = lw, col = bh_col)
lines(R, rck_opt_xs, ylim = c(0,1),
      xlab="", ylab="", lty=3,
      type = 'l', xaxs='i', yaxs='i', cex.axis = ca,
      lwd = lw*1.7, col = rck_col)
lines(R, rep(0.4,length(R)), ylim = c(0,1),
      xlab="", ylab="", lty=4,
      type = 'l', xaxs='i', yaxs='i', cex.axis = ca,
      lwd = lw*1.4, col = pt_col)
legend('topright',
       c('Beverton-Holt', 'Hockey-stick','Ricker','Pella-Tomlinson'),
       col = c(bh_col,hs_col,rck_col,pt_col),
       lty = c(1,2,3,4),
       bty='n', cex=.9, lwd = c(lw,lw*1.2,lw*1.7,lw*1.5))
mtext(side=1, "Population growth multiplier, r", outer = TRUE, padj=1.9, cex = cl)
mtext(side=2, "Optimal Escapement", outer = TRUE, padj=-1.9, cex = cl)
dev.off()

#### Fig 3 - Catch when model is wrong ####

# hockey-stick true
  #bev-holt decision given hockey-stick is true, biomass & catch
  x2.B.H = HS(r=R, x=e.B)
  C.B.H = x2.B.H - e.B
  #B60 decision given hockey-stick is true, biomass & catch
  x2.60.H = HS(r=R, x=B60)
  C.60.H = x2.60.H - B60
  #1/r hockey-stick decision given hockey-stick is true, biomass & catch
  x2.H.H = HS( r=R, x=e.H)
  C.H.H = x2.H.H - e.H

# Bev-Holt true
  #bev-holt decision given Bev-Holt is true, biomass & catch
  x2.B.B = BH(x=e.B, r=R)
  C.B.B = x2.B.B - e.B
  #B60 decision given Bev-Holt is true, biomass & catch
  x2.60.B = BH(x=B60, r=R)
  C.60.B = x2.60.B - B60
  #1/r hockey-stick decision given Bev-Holt is true, biomass & catch
  x2.H.B = BH(x=e.H, r=R)
  C.H.B = x2.H.B - e.H

pdf(file = './plots/catch.pdf', width = 7, height = 3.5)
  par(mfrow = c(1,2), mar = c(1,1,.5,.55), oma= c(2.51,2.51,0,0)); 
  lw = 2; ca=1.2; cl = 1.3;
  ylimC = c(0, max(C.B.B,C.H.H));
  ylimB = c(0,1.01)
  #BH true
  plot(R, C.B.B, ylim = ylimC,
       xlab="", ylab="",
       type = 'l', xaxs='i', yaxs='i', cex.axis = ca,
       lwd = lw, col = bh_col)
  #lines(R, C.60.B, lty = 2, lwd = lw)
  lines(R, C.H.B, lty = 2, lwd = lw*1.5, col = hs_col)
  text(1.92,.7,'a) Beverton-Holt true')
  #HS true
  plot(R, C.H.H,  ylim = ylimC,
       xlab="", ylab="", yaxt = 'n',
       type = 'l', lty = 2, col = hs_col,
       xaxs='i', yaxs='i', cex.axis = ca, lwd = lw*1.5)
  #lines(R, C.60.H, lty = 2, lwd = lw)
  lines(R, C.B.H, lty = 1, lwd = lw, col = bh_col)
  text(1.92,.7,'b) Hockey-Stick true')
  mtext(side=1, "Population growth multiplier, r", outer = TRUE, padj=1.9, cex = cl)
  mtext(side=2, "Catch", outer = TRUE, padj=-1.9, cex = cl)
  legend('bottomright',
         c('Beverton-Holt harvest', 'Hockey-stick harvest'),
         col = c(bh_col,hs_col),
         lty = c(1,2),
         bty='n', cex=.9, lwd = c(lw,lw*1.2))
dev.off()


#### Fig 4 - Proportional catch when model is wrong ####
pdf(file = './plots/prop_catch.pdf', width = 3.5, height = 3.5)
par(mfrow = c(1,1), mar = c(1,1,.7,.7), oma= c(2.51,2.51,0,0)); 
plot(R, C.H.B/C.B.B, ylim = c(0,1),
     xlab="", ylab="", lty=2,
     type = 'l', xaxs='i', yaxs='i', cex.axis = ca,
     lwd = lw*1.2, col = hs_col)
lines(R, C.B.H/C.H.H, ylim = c(0,1),
     xlab="", ylab="", lty=1,
     type = 'l', xaxs='i', yaxs='i', cex.axis = ca,
     lwd = lw, col = bh_col)
legend('bottomright',
       c('Beverton-Holt harvest', 'Hockey-stick harvest'),
       col = c(bh_col,hs_col),
       lty = c(1,2),
       bty='n', cex=.9, lwd = c(lw,lw*1.2))
mtext(side=1, "Population growth multiplier, r", outer = TRUE, padj=1.9, cex = cl)
mtext(side=2, "Proportion of optimal catch", outer = TRUE, padj=-1.9, cex = cl)
dev.off()


#### Fig 5 - Growth rates scatter plot####
bh_vs_hs_r <- function(r){
  return((r-1)/((r**0.5)-1))
}

R = seq(1, 5, length.out = 100)
pdf(file = './plots/scatter_plot_bh_hs_r.pdf', width = 3.5, height = 3.5)
par(mfrow = c(1,1), mar = c(1,1,.5,.55), oma= c(2.51,2.51,0,0)); 
plot(r_bh, r_hs, cex = 0.3, pch = 19, xlim = c(1, 4), ylim = c(1, 4),
     xlab="", ylab="", xaxs='i', yaxs='i', cex.axis = ca, xaxt = 'n', yaxt = 'n',)
lines(R, bh_vs_hs_r(R), 
      lty=1, lwd = lw*1.2)
lines(R, R, 
      lty = 2, lwd = lw*1.5, col = 'red')
mtext(side=1, "Beverton-Holt, r", outer = TRUE, padj=1.9, cex = cl)
mtext(side=2, "Hockey-Stick, r", outer = TRUE, padj=-1.9, cex = cl)
axis(side = 2, at = c(0, 1, 2, 3, 4), labels = c(0, 1, 2, 3, 4), cex.axis = ca)
axis(side = 1, at = c(0, 1, 2, 3, 4), labels = c(0, 1, 2, 3, 4), cex.axis = ca)
dev.off()


#### Fig 6 - Carrying capacity scatter plot ####
K = seq(0, 5, length.out = 100)
pdf(file = './plots/scatter_plot_bh_hs_k.pdf', width = 3.5, height = 3.5)
par(mfrow = c(1,1), mar = c(1,1,.5,.55), oma= c(2.51,2.51,0,0)); 
plot(k_bh, k_hs, cex = 0.3, pch = 19,
     xlim = c(0, 4), ylim = c(0, 1.1),
     xlab="", ylab="", xaxs='i', yaxs='i', cex.axis = ca, xaxt = 'n', yaxt = 'n')
lines(K, K,
      lty=2, lwd = lw*1.5, col = 'red')
lines(K, 0.694*K,
      lty=1, lwd = lw*1.2)
mtext(side=1, "Beverton-Holt, k", outer = TRUE, padj=1.9, cex = cl)
mtext(side=2, "Hockey-Stick, k", outer = TRUE, padj=-1.9, cex = cl)
axis(side = 2, at = c(0, 0.25, 0.5, 0.75, 1), labels = c(0, 0.25, 0.5, 0.75, 1), cex.axis = ca)
axis(side = 1, at = c(0, 1, 2, 3, 4), labels = c(0, 1, 2, 3, 4), cex.axis = ca)
dev.off()
#### Fig 7 - Histogram Prop fitted K Optimal Escapement####
pdf(file = './plots/opt_esc_hist_prop_of_fitted_k.pdf', width = 3.5, height = 3.5)
# par(mfrow=c(4,1), oma = c(6, 6, 0, 1), mar = c(0, 1, 0, 0)) 
par(mfrow = c(4,1), mar = c(0,0,.7,.7), oma= c(3.5,3.7,0,0));
hist(opt_esc_bh,
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
     cex.axis = 2)
hist(opt_esc_bh_bf,
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
     add = TRUE)
text(0.5, 80, 'Beverton-Holt', cex = 1.5)
axis(side = 2, at = c(0, 75), labels = c(0, 75), cex.axis = 1.25, mgp = c(2, 0.5, 0))
axis(1,pos=0, cex.axis = 1.25, mgp = c(2, 0.5, 0))
#lines(c(-1, 1), c(0, 0))
hist(opt_esc_rck,
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
     col = rck_col_tr,
     cex.axis = 2)
hist(opt_esc_rck_bf,
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
     col = rck_col,
     cex.axis = 2,
     add = TRUE)
text(0.5, 155, 'Ricker', cex = 1.5)
axis(side = 2, at = c(0, 150), labels = c(0, 150), cex.axis = 1.25, mgp = c(2, 0.5, 0))
axis(1,pos=0, cex.axis = 1.25, mgp = c(2, 0.5, 0))
hist(opt_esc_hs,
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
     cex.axis = 2)
hist(opt_esc_hs_bf,
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
     add = TRUE)
text(0.5, 25, 'Hockey-Stick', cex = 1.5)
axis(side = 2, at = c(0, 25), labels = c(0, 25), cex.axis = 1.25, mgp = c(2, 0.5, 0))
axis(1,pos=0, cex.axis = 1.25, mgp = c(2, 0.5, 0))
hist(opt_esc_pt,
     yaxs = "i",
     yaxt = 'n',
     xlim = c(0, 1),
     ylim = c(0, 800),
     breaks = seq(from = -100, to = 100, by = 0.025),
     main = '',
     xlab = '',
     ylab = '',
     xaxt = 'n',
     col = pt_col_tr,
     cex.axis = 2)
hist(opt_esc_pt_bf,
     yaxs = "i",
     yaxt = 'n',
     xlim = c(0, 1),
     ylim = c(0, 800),
     breaks = seq(from = -100, to = 100, by = 0.025),
     main = '',
     xlab = '',
     ylab = '',
     xaxt = 'n',
     col = pt_col,
     cex.axis = 2,
     add = TRUE)
text(0.5, 400, 'Pella-Tomlinson', cex = 1.5)
axis(side = 2, at = c(0,300), labels = c(0,300), cex.axis = 1.25, mgp = c(2, 0.5, 0))
axis(1,pos=0, cex.axis = 1.25, mgp = c(2, 0.5, 0))
mtext(side=1, "Optimal escapement", outer = TRUE, padj=1.9, cex = cl)
mtext(side=2, "Frequency", outer = TRUE, padj=-1.9, cex = cl)
dev.off()


#### Fig 8 - Histogram Prop Max Pop Optimal Escapement ####
pdf(file = './plots/opt_esc_hist_prop_of_max_pop.pdf', width = 3.5, height = 3.5)
par(mfrow = c(4,1), mar = c(0,0,.7,.7), oma= c(3.5,3.7,0,0));
hist(opt_esc_bh_scaled,
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
     cex.axis = 2)
hist(opt_esc_bh_scaled_bf,
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
     add = TRUE)
text(1.25, 50, 'Beverton-Holt', cex = 1.5)
axis(side = 2, at = c(0, 50), labels = c(0, 50), cex.axis = 1.25, mgp = c(2, 0.5, 0))
axis(1,pos=0, cex.axis = 1.25, mgp = c(2, 0.5, 0))
lines(c(-1, 1), c(0, 0))
hist(opt_esc_rck_scaled,
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
     col = rck_col_tr,
     cex.axis = 2)
hist(opt_esc_rck_scaled_bf,
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
     col = rck_col,
     cex.axis = 2,
     add = TRUE)
text(1.25, 50, 'Ricker', cex = 1.5)
axis(side = 2, at = c(0, 50), labels = c(0, 50), cex.axis = 1.25, mgp = c(2, 0.5, 0))
axis(1,pos=0, cex.axis = 1.25, mgp = c(2, 0.5, 0))
hist(opt_esc_hs_scaled,
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
     cex.axis = 2)
hist(opt_esc_hs_scaled_bf,
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
     add = TRUE)
text(1.25, 50, 'Hockey-Stick', cex = 1.5)
axis(side = 2, at = c(0, 50), labels = c(0, 50), cex.axis = 1.25, mgp = c(2, 0.5, 0))
axis(1,pos=0, cex.axis = 1.25, mgp = c(2, 0.5, 0))
hist(opt_esc_pt_scaled,
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
     cex.axis = 2)
hist(opt_esc_pt_scaled_bf,
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
     add = TRUE)
text(1.25, 50, 'Pella-Tomlinson', cex = 1.5)
axis(side = 2, at = c(0,50), labels = c(0, 50), cex.axis = 1.25, mgp = c(2, 0.5, 0))
axis(1,pos=0, cex.axis = 1.25, mgp = c(2, 0.5, 0))
mtext(side=1, "Optimal escapement", outer = TRUE, padj=1.9, cex = cl)
mtext(side=2, "Frequency", outer = TRUE, padj=-1.9, cex = cl)
dev.off()

#### Fig 9 - Accumulation best model fit####
tol = seq(from = 1, to = 5, length.out = 100)
for(j in tol){
  num_b_h <- 0
  num_h_s <- 0
  num_rck <- 0
  num_p_t <- 0
  for(i in 1:num_cols){
    if(r2_bh[i] == max(c(r2_bh[i], r2_rck[i], r2_hs[i], r2_pt[i])) & r_bh[i] < j){
      num_b_h = num_b_h + 1
    }
    else if(r2_rck[i] == max(c(r2_bh[i], r2_rck[i], r2_hs[i], r2_pt[i])) & r_rck[i] < j){
      num_rck = num_rck + 1
    }
    else if(r2_hs[i] == max(c(r2_bh[i], r2_rck[i], r2_hs[i], r2_pt[i])) & r_hs[i] < j){
      num_h_s = num_h_s + 1
    }
    else if(r2_pt[i] == max(c(r2_bh[i], r2_rck[i], r2_hs[i], r2_pt[i])) & r_pt[i] < j){
      num_p_t = num_p_t + 1
    }
  }
  num_bh_rs = append(num_bh_rs, num_b_h)
  num_hs_rs = append(num_hs_rs, num_h_s)
  num_rck_rs = append(num_rck_rs, num_rck)
  num_pt_rs = append(num_pt_rs, num_p_t)
}

pdf(file = './plots/cumulative_frequency_plot.pdf', width = 3.5, height = 3.5)
par(mfrow = c(1,1), mar = c(1,1,.7,.7), oma= c(2.51,2.51,0,0));
plot(tol, num_bh_rs, ylim = c(0,102),
     xlab="", ylab="", lty=1,
     type = 'l', xaxs='i', yaxs='i', cex.axis = ca,
     lwd = lw*1.2, col = bh_col)
lines(tol, num_hs_rs, ylim = c(0,102),
      xlab="", ylab="", lty=2,
      type = 'l', xaxs='i', yaxs='i', cex.axis = ca,
      lwd = lw, col = hs_col)
lines(tol, num_rck_rs, ylim = c(0,102),
      xlab="", ylab="", lty=3,
      type = 'l', xaxs='i', yaxs='i', cex.axis = ca,
      lwd = lw*1.7, col = rck_col)
lines(tol, num_pt_rs, ylim = c(0,102),
      xlab="", ylab="", lty=4,
      type = 'l', xaxs='i', yaxs='i', cex.axis = ca,
      lwd = lw*1.4, col = pt_col)
legend('bottomright',
       c('Beverton-Holt', 'Hockey-stick','Ricker','Pella-Tomlinson'),
       col = c(bh_col,hs_col,rck_col,pt_col),
       lty = c(1,2,3,4),
       bty='n', cex=.8, lwd = c(lw,lw*1.2,lw*1.7,lw*1.5))
mtext(side=1, "Population growth multiplier, r", outer = TRUE, padj=1.9, cex = cl)
mtext(side=2, "Cumulative Frequency", outer = TRUE, padj=-1.9, cex = cl)
dev.off()

}