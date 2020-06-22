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
a <- 1/3
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

## Fit the model
int_did <- matrix(0, nrow = 1, ncol = 1)
w_did <- matrix(1 / (N - 1), nrow = N - 1, ncol = 1)
Y_did <- matrix(0, nrow = T, ncol = 1)
Y_did <- matrix(0, nrow = T, ncol = 1)
# West Germany
i <- 1
Y1 <- as.matrix(Y[,i])
Y0 <- as.matrix(Y[,-i])
Z1 <- as.matrix(Z[,i])
Z0 <- as.matrix(Z[,-i])
X1 <- as.matrix(X[,i])
X0 <- as.matrix(X[,-i])
# Fit did
w <- w_did
int_did <- as.matrix(mean(Z1) - mean(Z0))
Y_did <- int_did[rep(1, T),] + Y0 %*% w
Y_true <- Y1

## Compute the standard errors
# Over units
std_err_i <- matrix(0, T1, N - 1)
for (i in 2:N) {
  Y1 <- as.matrix(Y[,i])
  Y0 <- as.matrix(Y[,-c(1,i)])
  Z1 <- as.matrix(Z[,i])
  Z0 <- as.matrix(Z[,-c(1,i)])
  X1 <- as.matrix(X[,i])
  X0 <- as.matrix(X[,-c(1,i)])
  
  # Fit did
  w <- matrix(1 / (N - 2), N - 2, 1)
  int <- as.matrix(mean(Z1) - mean(Z0))
  std_err_i[,i - 1] <- (Y1[-c(1:T0),] - int[rep(1, T1),] - Y0[-c(1:T0),] %*% w) ^ 2
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
  Z1 <- as.matrix(Z[c(1:(T0 - t)),1])
  Z0 <- as.matrix(Z[c(1:(T0 - t)),-1])
  
  # Fit did
  w <- matrix(1 / (N - 1), N - 1, 1)
  int <- as.matrix(mean(Z1) - mean(Z0))
  std_err_t[t,1] <- (Y1[T0 - t + 1,] - int - Y0[T0 - t + 1,] %*% w) ^ 2
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
    Z1 <- as.matrix(Z[c(1:(T0 - t)),i])
    Z0 <- as.matrix(Z[c(1:(T0 - t)),-c(1,i)])
    
    # Fit did
    w <- matrix(1 / (N - 2), N - 2, 1)
    int <- as.matrix(mean(Z1) - mean(Z0))
    std_err_temp[t,1] <- (Y1[T0 - t + 1,] - int - Y0[T0 - t + 1,] %*% w) ^ 2
  }
  std_err_temp <- as.matrix(apply(std_err_temp, 2, mean))
  std_err_it[i - 1,1] <- std_err_temp
}
std_err_it <- as.matrix(sqrt(apply(std_err_it, 2, mean)))

# Copy the standard errors
std_err_did_i <- std_err_i
std_err_did_t <- std_err_t
std_err_did_it <- std_err_it

## Save the results
save(list = c("w_did", "int_did", "Y_did", "Y_true",
              "std_err_did_i", "std_err_did_t", "std_err_did_it"), 
     file = "germ_did_nocov_05.RData")
writeMat("germ_did_nocov_05.mat", 
         w_did = w_did, int_did = int_did, 
         Y_did = Y_did, Y_true = Y_true,
         std_err_did_i = std_err_i, 
         std_err_did_t = std_err_t, 
         std_err_did_it = std_err_it)

