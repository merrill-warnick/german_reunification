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

## Load the data
#load("germany_data.RData")
d <- load("repgermany.Rdata")

#### figure it out

# Predictors
# "industry- mean 1971:1980
# schooling 1970 and 1975
# invest70 1980
# "gdp","trade","infrate"
# dependent "gdp"
# 1971:1980 time predictors prior
# time optimize ssr: 1981:1990
# treatment identifier: 7

# X1: matrix of treated predictor matrix (# predi x 1)
# X0: matrix of controls' predictor matrix(# predi x # of controls)

# X1
num70 <- which(d$year==1970)
num75 <- which(d$year==1975)
num80 <- which(d$year==1980)
indus_means <- rep(NA, length(unique(d$country)))
for(i in 1:unique(d$country)){
  indus_means[i] = mean(d$industry[(num70[i]+1):(num70[i]+10)])
}

special_predictors <- cbind(indus_means, d$schooling[num70], d$schooling[num75], d$invest70[num80])
X1 <- matrix(c(d$gdp[(num70[7]+1):(num70[7]+10)], d$trade[(num70[7]+1):(num70[7]+10)], d$infrate[(num70[7]+1):(num70[7]+10)], special_predictors[7,]), ncol = 1)

# X0
X0 <- matrix(NA, nrow = length(X1), ncol = length(unique(d$country))-1)
control_units <- unique(d$index)[-7]
j <- 1
for(i in control_units){
  X0[,j] = c(d$gdp[(num70[i]+1):(num70[i]+10)], d$trade[(num70[i]+1):(num70[i]+10)], d$infrate[(num70[i]+1):(num70[i]+10)], special_predictors[i,])
  j = j+1
}

# Z1: matrix of treated outcome data for the pre-treatment periods over which MPSE is to be minimized
Z1 <- d$gdp[(num70[7]-10):(num70[7]+20)]

# Z0
Z0 <- matrix(NA, nrow = length(Z1), ncol = length(unique(d$country))-1)
j <- 1
for(i in control_units){
  Z0[,j] = d$gdp[(num70[i]-10):(num70[i]+20)]
  j = j+1
}

# Y1: matrix of outcome data for treated unit
Y1 <- d$gdp[(num70[7]+21):(num70[7]+33)]

# Y0
Y0 <- matrix(NA, nrow = length(Y1), ncol = length(unique(d$country))-1)
j <- 1
for(i in control_units){
  Y0[,j] = d$gdp[(num70[i]+21):(num70[i]+33)]
  j = j+1
}


#################################
#################################

load("repgermany.Rdata")
d <- x

#########################
#########################
#########################
# TRY 2

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

#######################################

X0 <- dataprep.out$X0
X1 <- dataprep.out$X1

Z1 <- dataprep.out$Z1
Z0 <- dataprep.out$Z0

Y1 <- dataprep.out$Y1plot
Y0 <- dataprep.out$Y0plot


# [X1,X0]
X <- cbind(X1, X0)

# [Y1,Y0]
Y <- cbind(Y1, Y0)

# [Z1,Z0]
Z <- cbind(Z1, Z0)

## Parameters
K <- dim(X)[1]
N <- dim(Y)[2]
T <- dim(Y)[1]
T0 <- dim(Z)[1]
T1 <- T - T0
T0_tr <- floor(T0 * 2 / 3)
T0_te <- T0 - T0_tr

# Normalize
div <- as.matrix(apply(X, 1, sd))
X <- X / div[,rep(1, N)]

## Find the optimal subset
# Iterate over i
n_max <- N - 1
n_grid <- c(0:min(T0_tr - 1, n_max, N - 2))
nn <- length(n_grid) 
err <- matrix(0, nrow = N - 1, ncol = nn)
c <- matrix(1, nrow = T0, ncol = 1)
for (i in 2:N) {
  cat('Unit', toString(i), '\n')
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
  for (n in 0:(nn - 1)) {
    cat('Subs', toString(n), '\n')
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
    err[i - 1,n + 1] <- mean(e ^ 2)
  }
}
# Optimal n
err <- apply(t(scale(t(err))), 2, mean)
ind_opt <- which.min(err)
n_opt <- n_grid[ind_opt]

## Fit the model
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
# Fit the model
V1 <- Z1
subs_n <- combn(c(1:(N - 1)), n_opt, simplify = FALSE)
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
int_subs <- as.matrix(int[1,j_opt])
w_subs <- w[,j_opt]
Y_subs <- int[rep(1, T),j_opt] + Y0 %*% w[,j_opt]
Y_true <- Y1

## Compute the model errors
err_subs <- matrix(0, nrow = N - 1, ncol = 1)
# Iterate over i
c <- matrix(1, nrow = T0, ncol = 1)
for (i in 2:N) {
  cat('(Error) Unit', toString(i), '\n')
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
  n <- n_opt
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
c <- matrix(1, nrow = T0, ncol = 1)
std_err_i <- matrix(0, T1, N - 1)
for (i in 2:N) {
  cat('(Std. Error) Over Unit i =', toString(i), '\n')
  Y1 <- as.matrix(Y[,i])
  Y0 <- as.matrix(Y[,-c(1,i)])
  Z1 <- as.matrix(Z[,i])
  Z0 <- as.matrix(Z[,-c(1,i)])
  X1 <- as.matrix(X[,i])
  X0 <- as.matrix(X[,-c(1,i)])
  
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
  std_err_i[,i - 1] <- (Y1[-c(1:T0),] - int[rep(1, T1),j_opt] - Y0[-c(1:T0),] %*% w[,j_opt]) ^ 2
}
std_err_i <- as.matrix(sqrt(apply(std_err_i, 1, mean)))

# Over time
s <- floor(T0 / 2)
std_err_t <- matrix(0, s, 1)
Y1 <- as.matrix(Y[,1])
Y0 <- as.matrix(Y[,-1])
X1 <- as.matrix(X[,1])
X0 <- as.matrix(X[,-1])
for (t in 1:s) {
  cat('(Std. Error) Over Time t =', toString(t), '\n')
  c <- matrix(1, nrow = T0 - t, ncol = 1)
  Z1 <- as.matrix(Z[c(1:(T0 - t)),1])
  Z0 <- as.matrix(Z[c(1:(T0 - t)),-1])
  
  # Fit subs
  V1 <- Z1
  subs_n <- combn(c(1:(N - 1)), n_opt, simplify = FALSE)
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
  std_err_t[t,1] <- (Y1[T0 - t + 1,] - int[,j_opt] - Y0[T0 - t + 1,] %*% w[,j_opt]) ^ 2
}
std_err_t <- as.matrix(sqrt(apply(std_err_t, 2, mean)))

# Over units and time
std_err_it <- matrix(0, N - 1, 1)
for (i in 2:N) {
  Y1 <- as.matrix(Y[,i])
  Y0 <- as.matrix(Y[,-c(1,i)])
  X1 <- as.matrix(X[,i])
  X0 <- as.matrix(X[,-c(1,i)])
  
  std_err_temp <- matrix(0, s, 1)
  for (t in 1:s) {
    cat('(Std. Error) Over Unit and Time ( i , t ) = (',toString(i), ',', 
        toString(t), ')\n')
    c <- matrix(1, nrow = T0 - t, ncol = 1)
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
    std_err_temp[t,1] <- (Y1[T0 - t + 1,] - int[,j_opt] - Y0[T0 - t + 1,] %*% w[,j_opt]) ^ 2
  }
  std_err_temp <- as.matrix(apply(std_err_temp, 2, mean))
  std_err_it[i - 1,1] <- std_err_temp
}
std_err_it <- as.matrix(sqrt(apply(std_err_it, 2, mean)))

# Copy the standard errors
std_err_subs_i <- std_err_i
std_err_subs_t <- std_err_t
std_err_subs_it <- std_err_it

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

