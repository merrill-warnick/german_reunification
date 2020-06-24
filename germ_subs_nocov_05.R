## Clear
rm(list = ls())

## Options
#setwd("/Users/nikolayd/Dropbox/Research/Synth/Germany/")

## Libraries
library(lars)
library(glmnet)
library(foreign)
library(reshape2)
library(R.matlab)

######################
####### DATA #########
######################

## Load the data- from website in Rdata format
load("repgermany.Rdata")
d <- x

##############
# Use dataprep function from synthetic control package to extract X0, X1, etc.
dataprep.out <-
  dataprep(
    foo = d,
    predictors    = c("gdp","trade","infrate"), 
    dependent     = "gdp",
    unit.variable = 1,
    time.variable = 3,
    special.predictors = list(
      list("industry" ,1981:1990, c("mean")),
      list("schooling",c(1980,1985), c("mean")),
      list("invest80" ,1980, c("mean"))
    ),
    treatment.identifier = 7,
    controls.identifier = unique(d$index)[-7],
    time.predictors.prior = 1981:1990,
    time.optimize.ssr = 1960:1989,
    unit.names.variable = 2,
    time.plot = 1960:2003
  )

# Predictors will be GDP, TRADE, INFRATE from 1981:1990 plus mean industry
# from 1981:190, schooling in 1980 and 1985 and invest80 in 1980
# treatment.identifier 7 corresponds to Germany
# controls identifier (all minus 7) corresponds to all other countries
# optimize.ssr: choose V s.t. over 1960 to 1989 MSPE is minimized???? -> look into this


##############
# Extract X0,X1, etc. from dataprep object

X0 <- dataprep.out$X0
X1 <- dataprep.out$X1

Z1 <- dataprep.out$Z1
Z0 <- dataprep.out$Z0

Y1 <- dataprep.out$Y1plot
Y0 <- dataprep.out$Y0plot

##############
# Combine to get X, Z and Y matrix

# Make sure that treated unit is at beginning! Essential for Std Er function!
# [X1,X0]
X <- cbind(X1, X0)

# [Y1,Y0]
Y <- cbind(Y1, Y0)

# [Z1,Z0]
Z <- cbind(Z1, Z0)

############################
####### PARAMETERS #########
############################

## Define Parameters
K <- dim(X)[1] # Number of predictors
N <- dim(Y)[2] # Number of units
T <- dim(Y)[1] # Number of time periods
T0 <- dim(Z)[1] # Time of intervention
T1 <- T - T0 # Number of time periods after intervention
T0_tr <- floor(T0 * 2 / 3) # what is this?? - here it is actually used
T0_te <- T0 - T0_tr

# Normalize predictors
div <- as.matrix(apply(X, 1, sd)) # Matrix of standard deviations for each predictor
X <- X / div[,rep(1, N)] # Standardizes each predictor to have std 1

################################
####### OPTIMAL SUBSET #########
################################

## Find the optimal subset #!!
# Iterate over i
n_max <- N - 1 # Number of units in subset ?
n_grid <- c(0:min(T0_tr - 1, n_max, N - 2)) # What is this grid?? 
nn <- length(n_grid) # Number of points in n grid
err <- matrix(0, nrow = N - 1, ncol = nn) # Storage for errors for each unit and n
c <- matrix(1, nrow = T0, ncol = 1) # What is c??

for (i in 2:N) { # over units
  cat('Unit', toString(i), '\n')
  
  # Fix matrices 
  Y1 <- as.matrix(Y[,i])
  Y0 <- as.matrix(Y[,-c(1,i)]) # Exclude West Germany as well
  Z1 <- as.matrix(Z[,i])
  Z0 <- as.matrix(Z[,-c(1,i)])
  X1 <- as.matrix(X[,i])
  X0 <- as.matrix(X[,-c(1,i)])
  Z1_tr <- as.matrix(Z1)
  Z0_tr <- as.matrix(Z0)
  Z1_te <- as.matrix(Y1[-(1:T0),])
  Z0_te <- as.matrix(Y0[-(1:T0),])
  
  # Fit the best-subset model
  V1 <- Z1_tr
  
  for (n in 0:(nn - 1)) { # iterate over n's
    cat('Subs', toString(n), '\n')
    subs_n <- combn(c(1:(N - 2)), n, simplify = FALSE) # generates all combinations of c() taking n at a time
    int <- matrix(0, nrow = 1, ncol = length(subs_n))
    w <- matrix(0, nrow = N - 2, ncol = length(subs_n)) # weight matrix for each of the N-2 units and column: each combination
    err_cur <- matrix(0, nrow = length(subs_n), ncol = 1) # Error for each subset
    
    for (j in 1:length(subs_n)) { # all possible combinations for this length of subset
      sub <- subs_n[[j]]
      V0 <- cbind(c, Z0_tr[,sub])
      w_cur <- solve(t(V0) %*% V0, t(V0) %*% V1)
      int[1,j] <- w_cur[1]
      w[sub,j] <- w_cur[-1]
      err_cur[j,1] <- mean((V1 - V0 %*% w_cur) ^ 2)
    }
    # Choose the optimal subset of size n and compute the error
    j_opt <- which.min(err_cur)
    e <- Z1_te - int[rep(1, T1),j_opt] - Z0_te %*% w[,j_opt]
    err[i - 1,n + 1] <- mean(e ^ 2)
  }
}

# Optimal n
err <- apply(t(scale(t(err))), 2, mean) # Get average error over all units
ind_opt <- which.min(err)
n_opt <- n_grid[ind_opt]

###########################
####### Fit Model #########
###########################

## Fit the model

# Initialize storage matrices
int_subs <- matrix(0, nrow = 1, ncol = 1)
w_subs <- matrix(0, nrow = N - 1, ncol = 1)
Y_subs <- matrix(0, nrow = T, ncol = 1)
Y_true <- matrix(0, nrow = T, ncol = 1)
c <- matrix(1, nrow = T0, ncol = 1)

# West Germany
i <- 1
Y1 <- as.matrix(Y[,i])
Y0 <- as.matrix(Y[,-i])
Z1 <- as.matrix(Z[,i])
Z0 <- as.matrix(Z[,-i])
X1 <- as.matrix(X[,i])
X0 <- as.matrix(X[,-i])

######
# Fit the model

V1 <- Z1
subs_n <- combn(c(1:(N - 1)), n_opt, simplify = FALSE) # Get all the subset combinations of size n*
int <- matrix(0, nrow = 1, ncol = length(subs_n))
w <- matrix(0, nrow = N - 1, ncol = length(subs_n))
err_cur <- matrix(0, nrow = length(subs_n), ncol = 1)

# Fit for each combination to get std error
for (j in 1:length(subs_n)) {
  sub <- subs_n[[j]]
  V0 <- cbind(c, Z0[,sub])
  w_cur <- solve(t(V0) %*% V0, t(V0) %*% V1)
  int[1,j] <- w_cur[1]
  w[sub,j] <- w_cur[-1]
  err_cur[j,1] <- mean((V1 - V0 %*% w_cur) ^ 2)
}

# Extract which subset is best
j_opt <- which.min(err_cur)

###########
# Fit best subset
int_subs <- as.matrix(int[1,j_opt])
w_subs <- w[,j_opt]
Y_subs <- int[rep(1, T),j_opt] + Y0 %*% w[,j_opt] # Estimated Y (no treatment)
Y_true <- Y1

#################################
####### Standard Errors #########
#################################

## Compute the model errors
err_subs <- matrix(0, nrow = N - 1, ncol = 1)

# Iterate over i (units)
c <- matrix(1, nrow = T0, ncol = 1) # What is this for??
for (i in 2:N) {
  cat('(Error) Unit', toString(i), '\n')
  
  # Define matrices
  Y1 <- as.matrix(Y[,i])
  Y0 <- as.matrix(Y[,-c(1,i)])
  Z1 <- as.matrix(Z[,i])
  Z0 <- as.matrix(Z[,-c(1,i)])
  X1 <- as.matrix(X[,i])
  X0 <- as.matrix(X[,-c(1,i)])
  Z1_tr <- as.matrix(Z1)
  Z0_tr <- as.matrix(Z0)
  Z1_te <- as.matrix(Y1[-(1:T0),])
  Z0_te <- as.matrix(Y0[-(1:T0),])
  
  # Fit the best-subset model
  V1 <- Z1_tr
  n <- n_opt # Fix optimal n
  subs_n <- combn(c(1:(N - 2)), n, simplify = FALSE)
  int <- matrix(0, nrow = 1, ncol = length(subs_n))
  w <- matrix(0, nrow = N - 2, ncol = length(subs_n))
  err_cur <- matrix(0, nrow = length(subs_n), ncol = 1)
  for (j in 1:length(subs_n)) {
    sub <- subs_n[[j]]
    V0 <- cbind(c, Z0_tr[,sub])
    w_cur <- solve(t(V0) %*% V0, t(V0) %*% V1)
    int[1,j] <- w_cur[1]
    w[sub,j] <- w_cur[-1]
    err_cur[j,1] <- mean((V1 - V0 %*% w_cur) ^ 2)
  }
  # Choose the optimal subset of size n and compute the error
  j_opt <- which.min(err_cur)
  e <- Z1_te - int[rep(1, T1),j_opt] - Z0_te %*% w[,j_opt]
  err_subs[i - 1] <- mean(e ^ 2)
}

## Compute the standard errors
# Over units
c <- matrix(1, nrow = T0, ncol = 1) # What is this for??
std_err_i <- matrix(0, T1, N - 1)

for (i in 2:N) {
  cat('(Std. Error) Over Unit i =', toString(i), '\n')
  
  # Define matrices
  Y1 <- as.matrix(Y[,i])
  Y0 <- as.matrix(Y[,-c(1,i)])
  Z1 <- as.matrix(Z[,i])
  Z0 <- as.matrix(Z[,-c(1,i)])
  X1 <- as.matrix(X[,i])
  X0 <- as.matrix(X[,-c(1,i)])
  
  # Fit subset
  V1 <- Z1
  subs_n <- combn(c(1:(N - 2)), n_opt, simplify = FALSE) # Note that c goes from 1 to N-2 because Z does not include Germany and pseudo unit 
  int <- matrix(0, nrow = 1, ncol = length(subs_n))
  w <- matrix(0, nrow = N - 2, ncol = length(subs_n))
  err_cur <- matrix(0, nrow = length(subs_n), ncol = 1)
  
  for (j in 1:length(subs_n)) {
    sub <- subs_n[[j]]
    V0 <- cbind(c, Z0[,sub])
    w_cur <- solve(t(V0) %*% V0, t(V0) %*% V1)
    int[1,j] <- w_cur[1]
    w[sub,j] <- w_cur[-1]
    err_cur[j,1] <- mean((V1 - V0 %*% w_cur) ^ 2)
  }
  j_opt <- which.min(err_cur)
  std_err_i[,i - 1] <- (Y1[-c(1:T0),] - int[rep(1, T1),j_opt] - Y0[-c(1:T0),] %*% w[,j_opt]) ^ 2 # SSR for each time point in T1
}
std_err_i <- as.matrix(sqrt(apply(std_err_i, 1, mean))) # Std Err is sqrt of mean over units at one point in time

# Over time
s <- floor(T0 / 2)
std_err_t <- matrix(0, s, 1)

# Fix X and Y matrices
Y1 <- as.matrix(Y[,1])
Y0 <- as.matrix(Y[,-1])
X1 <- as.matrix(X[,1])
X0 <- as.matrix(X[,-1])

for (t in 1:s) {
  cat('(Std. Error) Over Time t =', toString(t), '\n')
  
  # Varying Z1 and Z0: change intervention point
  c <- matrix(1, nrow = T0 - t, ncol = 1) # What is this for??
  Z1 <- as.matrix(Z[c(1:(T0 - t)),1])
  Z0 <- as.matrix(Z[c(1:(T0 - t)),-1])
  
  # Fit subs
  V1 <- Z1
  subs_n <- combn(c(1:(N - 1)), n_opt, simplify = FALSE) # N-1 again because only Germany drops out
  int <- matrix(0, nrow = 1, ncol = length(subs_n))
  w <- matrix(0, nrow = N - 1, ncol = length(subs_n))
  err_cur <- matrix(0, nrow = length(subs_n), ncol = 1)
  
  for (j in 1:length(subs_n)) {
    sub <- subs_n[[j]]
    V0 <- cbind(c, Z0[,sub])
    w_cur <- solve(t(V0) %*% V0, t(V0) %*% V1)
    int[1,j] <- w_cur[1]
    w[sub,j] <- w_cur[-1]
    err_cur[j,1] <- mean((V1 - V0 %*% w_cur) ^ 2)
  }
  j_opt <- which.min(err_cur)
  std_err_t[t,1] <- (Y1[T0 - t + 1,] - int[,j_opt] - Y0[T0 - t + 1,] %*% w[,j_opt]) ^ 2 # SSR for each t in s
}
std_err_t <- as.matrix(sqrt(apply(std_err_t, 2, mean))) # Sqrt of SSR over time

# Over units and time
std_err_it <- matrix(0, N - 1, 1)

for (i in 2:N) {
  
  # Over units
  Y1 <- as.matrix(Y[,i])
  Y0 <- as.matrix(Y[,-c(1,i)])
  X1 <- as.matrix(X[,i])
  X0 <- as.matrix(X[,-c(1,i)])
  
  std_err_temp <- matrix(0, s, 1)
  
  # Over time
  for (t in 1:s) {
    cat('(Std. Error) Over Unit and Time ( i , t ) = (',toString(i), ',', 
        toString(t), ')\n')
    
    # Time-varying Z matrix (pre-outcome matrices)
    c <- matrix(1, nrow = T0 - t, ncol = 1) # Why do we need this??
    Z1 <- as.matrix(Z[c(1:(T0 - t)),i])
    Z0 <- as.matrix(Z[c(1:(T0 - t)),-c(1,i)])
    
    # Fit subs
    V1 <- Z1
    subs_n <- combn(c(1:(N - 2)), n_opt, simplify = FALSE)
    int <- matrix(0, nrow = 1, ncol = length(subs_n))
    w <- matrix(0, nrow = N - 2, ncol = length(subs_n))
    err_cur <- matrix(0, nrow = length(subs_n), ncol = 1)
    for (j in 1:length(subs_n)) {
      sub <- subs_n[[j]]
      V0 <- cbind(c, Z0[,sub])
      w_cur <- solve(t(V0) %*% V0, t(V0) %*% V1)
      int[1,j] <- w_cur[1]
      w[sub,j] <- w_cur[-1]
      err_cur[j,1] <- mean((V1 - V0 %*% w_cur) ^ 2)
    }
    j_opt <- which.min(err_cur)
    std_err_temp[t,1] <- (Y1[T0 - t + 1,] - int[,j_opt] - Y0[T0 - t + 1,] %*% w[,j_opt]) ^ 2 # SSR for each time point t
  }
  std_err_temp <- as.matrix(apply(std_err_temp, 2, mean)) # Take mean over time points to get std err
  std_err_it[i - 1,1] <- std_err_temp
}
std_err_it <- as.matrix(sqrt(apply(std_err_it, 2, mean))) # Take mean over units to get finalized stad error

# Copy the standard errors
std_err_subs_i <- std_err_i
std_err_subs_t <- std_err_t
std_err_subs_it <- std_err_it

##############################
####### Save Results #########
##############################

## Save the results
save(list = c("w_subs", "int_subs", "Y_subs", "Y_true", "n_opt", 
              "std_err_subs_i", "std_err_subs_t", "std_err_subs_it",
              "err_subs"), 
     file = "germ_subs_nocov_05.RData")
writeMat("germ_subs_nocov_05.mat", 
         w_subs = w_subs, int_subs = int_subs, 
         Y_subs = Y_subs, Y_true = Y_true,
         n_opt = n_opt, 
         std_err_subs_i = std_err_i, 
         std_err_subs_t = std_err_t, 
         std_err_subs_it = std_err_it,
         err_subs = err_subs)