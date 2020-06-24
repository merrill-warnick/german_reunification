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

##############################
##############################
##############################
# DATA TRIAL 
##############################
##############################
##############################

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
    
    #this matches what we're doing in synth!
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

#okay so counterfactual: we're setting year 1980 for the counterfactual exercise.
## Counterfactual
T0_co <- 21
T1_co <- T0 - T0_co

# Normalize
#hmmmmm why are we normalizing? In the paper we say that we don't want to normalize Y^obs_c,pre and I thinkt hat that's what we're doing right here.
#I'll talk to Lea about this one also
div <- as.matrix(apply(X, 1, sd))
X <- X / div[,rep(1, N)]

## Find the optimal elasticity

#I'm a teeny bit confused about how this code works exactly but I get the idea, I just don't get all of the mechanics

# Iterate over i
lambda_grid <- c(seq(from = 1e-02, to = 1e-01, by = 1e-02),
                 seq(from = 2e-01, to = 100, by = 1e-01), 
                 seq(from = 200, to = 50000, by = 100))
a_grid <- seq(from = 0.1, to = 0.9, by = 0.1)
nlambda <- length(lambda_grid)
na <- length(a_grid)
err_alpha <- matrix(0, nrow = na, ncol = 1)
lambda_opt_alpha <- matrix(0, nrow = na, ncol = 1)
cat('*** Main ***\n')
for (j in 1:na) {
  a <- a_grid[j]
  cat('a =', toString(a), '\n')
  err <- matrix(0, nrow = N - 1, ncol = nlambda)
  for (i in 2:N) {
    Y1 <- as.matrix(Y[,i])
    Y0 <- as.matrix(Y[,-c(1,i)])
    Z1 <- as.matrix(Z[,i])
    Z0 <- as.matrix(Z[,-c(1,i)])
    X1 <- as.matrix(X[,i])
    X0 <- as.matrix(X[,-c(1,i)])
    Z1_tr <- Z1
    Z0_tr <- Z0
    Z1_te <- as.matrix(Y1[-(1:T0),])
    Z0_te <- as.matrix(Y0[-(1:T0),])
    
    # Fit elast
    V1 <- scale(Z1_tr, scale = FALSE)
    V0 <- scale(Z0_tr, scale = FALSE)
    fit <- glmnet(x = V0, y = V1,
                  alpha = a,
                  lambda = lambda_grid,
                  standardize = FALSE,
                  intercept = FALSE)
    w <- as.matrix(coef(fit, s = lambda_grid))
    w <- w[-1,]
    int <- t(as.matrix(apply(Z1_tr[,rep(1, nlambda)] - Z0_tr %*% w, 2, mean)))
    e <- Z1_te[,rep(1, nlambda)] - int[rep(1, T1),] - Z0_te %*% w
    err[i - 1,] <- colMeans(e ^ 2)
  }
  # Optimal lambda
  #err <- apply(t(scale(t(err))), 2, mean)
  err <- apply(err, 2, mean)
  #ind_opt <- max(which(err <= min(err) + 1 * sd(err)))
  ind_opt <- which.min(err)
  err_alpha[j] <- err[ind_opt]
  lambda_opt_alpha[j] <- lambda_grid[ind_opt]
}
# Optimal a
ind_opt <- which.min(err_alpha)
a_opt <- a_grid[ind_opt]
lambda_opt <- lambda_opt_alpha[ind_opt]

## Fit the model
int_elast <- matrix(0, nrow = 1, ncol = 1)
w_elast <- matrix(0, nrow = N - 1, ncol = 1)
Y_elast <- matrix(0, nrow = T, ncol = 1)
Y_true <- matrix(0, nrow = T, ncol = 1)
# West Germany
i <- 1
Y1 <- as.matrix(Y[,i])
Y0 <- as.matrix(Y[,-i])
Z1 <- as.matrix(Z[,i])
Z0 <- as.matrix(Z[,-i])
X1 <- as.matrix(X[,i])
X0 <- as.matrix(X[,-i])

#hmmm I guess it's fine since the xs get scaled but the zs never get scaled
V1 <- scale(Z1, scale = FALSE)
V0 <- scale(Z0, scale = FALSE)
# Fit elast
fit <- glmnet(x = V0, y = V1,
              alpha = a_opt,
              lambda = lambda_grid,
              standardize = FALSE,
              intercept = FALSE)
w <- as.matrix(coef(fit, s = lambda_opt))
w <- w[-1,]
int_elast <- as.matrix(apply(Z1 - Z0 %*% w, 2, mean))
w_elast <- w
Y_elast <- int_elast[rep(1, T),] + Y0 %*% w
Y_true <- Y1

#this standard errors procedure looks the same in all the files.

## Compute the standard errors
# Over units
std_err_i <- matrix(0, N - 1, T1)
for (i in 2:N) {
  cat('(Std. Error) Over Unit i =', toString(i), '\n')
  Y1 <- as.matrix(Y[,i])
  Y0 <- as.matrix(Y[,-c(1,i)])
  Z1 <- as.matrix(Z[,i])
  Z0 <- as.matrix(Z[,-c(1,i)])
  X1 <- as.matrix(X[,i])
  X0 <- as.matrix(X[,-c(1,i)])
  
  # Fit elast
  V1 <- scale(Z1, scale = FALSE)
  V0 <- scale(Z0, scale = FALSE)
  fit <- glmnet(x = V0, y = V1,
                alpha = a_opt,
                lambda = lambda_grid,
                standardize = FALSE,
                intercept = FALSE)
  w <- as.matrix(coef(fit, s = lambda_opt))
  w <- w[-1,]
  int <- as.matrix(apply(Z1 - Z0 %*% w, 2, mean))
  std_err_i[i - 1,] <- (Y1[-c(1:T0),] - int[rep(1,T1),] - Y0[-c(1:T0),] %*% w) ^ 2
}
std_err_i <- as.matrix(sqrt(apply(std_err_i, 2, mean)))

# Over time
s <- floor(T0 / 2)
std_err_t <- matrix(0, s, 1)
Y1 <- as.matrix(Y[,1])
Y0 <- as.matrix(Y[,-1])
X1 <- as.matrix(X[,1])
X0 <- as.matrix(X[,-1])
for (t in 1:s) {
  cat('(Std. Error) Over Time t =', toString(t), '\n')
  Z1 <- as.matrix(Z[c(1:(T0 - t)),1])
  Z0 <- as.matrix(Z[c(1:(T0 - t)),-1])
  
  # Fit elast
  V1 <- scale(Z1, scale = FALSE)
  V0 <- scale(Z0, scale = FALSE)
  fit <- glmnet(x = V0, y = V1,
                alpha = a_opt,
                lambda = lambda_grid,
                standardize = FALSE,
                intercept = FALSE)
  w <- as.matrix(coef(fit, s = lambda_opt))
  w <- w[-1,]
  int <- as.matrix(apply(Z1 - Z0 %*% w, 2, mean))
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
    cat('(Std. Error) Over Unit and Time ( i , t ) = (',toString(i), ',', 
        toString(t), ')\n')
    Z1 <- as.matrix(Z[c(1:(T0 - t)),i])
    Z0 <- as.matrix(Z[c(1:(T0 - t)),-c(1,i)])
    
    # Fit elast
    V1 <- scale(Z1, scale = FALSE)
    V0 <- scale(Z0, scale = FALSE)
    fit <- glmnet(x = V0, y = V1,
                  alpha = a_opt,
                  lambda = lambda_grid,
                  standardize = FALSE,
                  intercept = FALSE)
    w <- as.matrix(coef(fit, s = lambda_opt))
    w <- w[-1,]
    int <- as.matrix(apply(Z1 - Z0 %*% w, 2, mean))
    std_err_temp[t,1] <- (Y1[T0 - t + 1,] - int - Y0[T0 - t + 1,] %*% w) ^ 2
  }
  std_err_temp <- as.matrix(apply(std_err_temp, 2, mean))
  std_err_it[i - 1,1] <- std_err_temp
}
std_err_it <- as.matrix(sqrt(apply(std_err_it, 2, mean)))

# Copy the standard errors
std_err_elast_i <- std_err_i
std_err_elast_t <- std_err_t
std_err_elast_it <- std_err_it

#then we do the counterfactual exercize

## Find the optimal counterfactual elast
# Iterate over i
lambda_grid <- c(seq(from = 1e-02, to = 1e-01, by = 1e-02),
                 seq(from = 2e-01, to = 100, by = 1e-01), 
                 seq(from = 200, to = 50000, by = 100))
a_grid <- seq(from = 0.1, to = 0.9, by = 0.1)
nlambda <- length(lambda_grid)
na <- length(a_grid)
err_alpha <- matrix(0, nrow = na, ncol = 1)
lambda_opt_alpha <- matrix(0, nrow = na, ncol = 1)
cat('*** Counterfactual ***\n')
for (j in 1:na) {
  a <- a_grid[j]
  cat('a =', toString(a), '\n')
  err <- matrix(0, nrow = N - 1, ncol = nlambda)
  for (i in 2:N) {
    Y1 <- as.matrix(Y[1:T0,i])
    Y0 <- as.matrix(Y[1:T0,-c(1,i)])
    Z1 <- as.matrix(Z[1:T0_co,i])
    Z0 <- as.matrix(Z[1:T0_co,-c(1,i)])
    X1 <- as.matrix(X[,i])
    X0 <- as.matrix(X[,-c(1,i)])
    Z1_tr <- Z1
    Z0_tr <- Z0
    Z1_te <- as.matrix(Y1[-(1:T0_co),])
    Z0_te <- as.matrix(Y0[-(1:T0_co),])
    
    # Fit elast
    V1 <- scale(Z1_tr, scale = FALSE)
    V0 <- scale(Z0_tr, scale = FALSE)
    fit <- glmnet(x = V0, y = V1,
                  alpha = a,
                  lambda = lambda_grid,
                  standardize = FALSE,
                  intercept = FALSE)
    w <- as.matrix(coef(fit, s = lambda_grid))
    w <- w[-1,]
    int <- t(as.matrix(apply(Z1_tr[,rep(1, nlambda)] - Z0_tr %*% w, 2, mean)))
    e <- Z1_te[,rep(1, nlambda)] - int[rep(1, T1_co),] - Z0_te %*% w
    err[i - 1,] <- colMeans(e ^ 2)
  }
  # Optimal lambda
  #err <- apply(t(scale(t(err))), 2, mean)
  err <- apply(err, 2, mean)
  #ind_opt <- max(which(err <= min(err) + 1 * sd(err)))
  ind_opt <- which.min(err)
  err_alpha[j] <- err[ind_opt]
  lambda_opt_alpha[j] <- lambda_grid[ind_opt]
}
# Optimal a
ind_opt <- which.min(err_alpha)
a_opt_co <- a_grid[ind_opt]
lambda_opt_co <- lambda_opt_alpha[ind_opt]

## Compute the counterfactual treatment
int_elast_co <- matrix(0, nrow = 1, ncol = 1)
w_elast_co <- matrix(0, nrow = N - 1, ncol = 1)
Y_elast_co <- matrix(0, nrow = T0, ncol = 1)
Y_true_co <- matrix(0, nrow = T0, ncol = 1)
# West Germany
i <- 1
Y1 <- as.matrix(Y[1:T0,i])
Y0 <- as.matrix(Y[1:T0,-i])
Z1 <- as.matrix(Z[1:T0_co,i])
Z0 <- as.matrix(Z[1:T0_co,-i])
X1 <- as.matrix(X[,i])
X0 <- as.matrix(X[,-i])
V1 <- scale(Z1, scale = FALSE)
V0 <- scale(Z0, scale = FALSE)
# Fit elast
fit <- glmnet(x = V0, y = V1,
              alpha = a_opt_co,
              lambda = lambda_grid,
              standardize = FALSE,
              intercept = FALSE)
w <- as.matrix(coef(fit, s = lambda_opt_co))
w <- w[-1,]
int_elast_co <- as.matrix(apply(Z1 - Z0 %*% w, 2, mean))
w_elast_co <- w
Y_elast_co <- int_elast_co[rep(1, T0),] + Y0 %*% w
Y_true_co <- Y1

## Compute the standard errors
# Over units
std_err_i_co <- matrix(0, N - 1, T1_co)
for (i in 2:N) {
  cat('(Std. Error) Over Unit i =', toString(i), '\n')
  Y1 <- as.matrix(Y[1:T0,i])
  Y0 <- as.matrix(Y[1:T0,-c(1,i)])
  Z1 <- as.matrix(Z[1:T0_co,i])
  Z0 <- as.matrix(Z[1:T0_co,-c(1,i)])
  X1 <- as.matrix(X[,i])
  X0 <- as.matrix(X[,-c(1,i)])
  
  # Fit elast
  V1 <- scale(Z1, scale = FALSE)
  V0 <- scale(Z0, scale = FALSE)
  fit <- glmnet(x = V0, y = V1,
                alpha = a_opt_co,
                lambda = lambda_grid,
                standardize = FALSE,
                intercept = FALSE)
  w <- as.matrix(coef(fit, s = lambda_opt_co))
  w <- w[-1,]
  int <- as.matrix(apply(Z1 - Z0 %*% w, 2, mean))
  std_err_i_co[i - 1,] <- (Y1[-c(1:T0_co),] - int[rep(1,T1_co),] - Y0[-c(1:T0_co),] %*% w) ^ 2
}
std_err_i_co <- as.matrix(sqrt(apply(std_err_i_co, 2, mean)))

# Copy the counterfactual standard errors
std_err_elast_i_co <- std_err_i_co

## Save the results
save(list = c("w_elast", "int_elast", "Y_elast", "Y_true", "a_opt", "lambda_opt", 
              "std_err_elast_i", "std_err_elast_t", "std_err_elast_it",
              "Y_elast_co", "Y_true_co", "std_err_elast_i_co", "T0_co", "T1_co"), 
     file = "germ_elast_nocov_05.RData")
writeMat("germ_elast_nocov_05.mat", 
         w_elast = w_elast, int_elast = int_elast, 
         Y_elast = Y_elast, Y_true = Y_true,
         a_opt = a_opt,
         lambda_opt = lambda_opt,
         std_err_elast_i = std_err_i, 
         std_err_elast_t = std_err_t, 
         std_err_elast_it = std_err_it,
         Y_elast_co = Y_elast_co,
         Y_true_co = Y_true_co,
         std_err_elast_i_co = std_err_i_co, 
         T0_co = T0_co,
         T1_co = T1_co)

