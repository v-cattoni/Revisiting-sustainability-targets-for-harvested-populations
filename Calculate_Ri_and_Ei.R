# Calculate Recruitment and Escapement

#### Reading in files ####
load("./data_files/RAMCore[asmt][v4.495].rdata") #RAM DATABASE
#this was visually identified manually recorded in the below csv file
cols <- read.csv("./data_files/non_det_cols.csv", header <- FALSE)[, 1]

#### Constants ####
years <- seq(1950, 2020, 1) #year of each row of data
num_rows <- length(years) #num rows to iterate through all data
num_cols <- length(cols) #num cols to iterate through all data

#### Initialising vectors and constants ####
bio_mat <- matrix(NA, num_rows, num_cols) # biomass (MT)
catch_mat <- matrix(NA, num_rows, num_cols) # catch (MT)
n_years <- rep(NA, num_cols) # number of datapoints
Ri <- rep(NA, length(cols)) # recruitment (MT)
Ei <- rep(NA, length(cols)) # escapement (MT)

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

#### Calculating Escapement and Recruitment ####
for (i in 1:length(cols)) {
  bios <- bio_mat[, i][!is.na(bio_mat[, i])]
  catchs <- catch_mat[, i][!is.na(catch_mat[, i])]
  X <-
    bios[1:(length(bios) - 1)] - catchs[1:(length(bios) - 1)]
  Y <- bios[2:length(bios)]
  
  n_years[i] <- length(X)
  
  Ri[i] <- list(c(Y))
  Ei[i] <- list(c(X))
}

#### Saving Escapement and Recruitment ####
saveRDS(Ri, file = "./data_files/Ri.rds")
saveRDS(Ei, file = "./data_files/Ei.rds")
