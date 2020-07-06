## Clear
rm(list = ls())

## Libraries
library(lars)
library(glmnet)
library(foreign)
library(reshape2)
library(R.matlab)
library(Synth)

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
T1 <- T - T0 # Time periods after intervention

# Define time period of intervention and time after intervention for counterfactual - why?
## Counterfactual
T0_co <- 21
T1_co <- T0 - T0_co

#####################################
####### Optimal elastic net #########
#####################################

## Find the optimal elastic net:
# optimal alpha and lambda selected from minimal standard error

# Iterate over i (units)

# Definitions
lambda_grid <- c(seq(from = 1e-02, to = 1e-01, by = 1e-02),
                 seq(from = 2e-01, to = 100, by = 1e-01), 
                 seq(from = 200, to = 50000, by = 100)) # non-linear lambda grid- why not just log-linear??/ how did you determine it needs 5000
a_grid <- seq(from = 0.1, to = 0.9, by = 0.1) # alpha for convex combination between two penalty terms: exlude 0 and 1
nlambda <- length(lambda_grid) # length of lambda grid
na <- length(a_grid) # length of a grid

err_alpha <- matrix(0, nrow = na, ncol = 1) # Matrix to store error terms associated with each value of alpha
lambda_opt_alpha <- matrix(0, nrow = na, ncol = 1) # Matrix to store optimal lambda associated with each alpha

# Iteration
cat('*** Main ***\n')
for (j in 1:na) { # Iterate over alpha points
  a <- a_grid[j]
  cat('a =', toString(a), '\n')
  err <- matrix(0, nrow = N - 1, ncol = nlambda) # Matrix for storage of error terms for each control unit and all lambda values
  for (i in 2:N) { # iterate over units
    
    # Determine matrices appropriately
    Y1 <- as.matrix(Y[,i])
    Y0 <- as.matrix(Y[,-c(1,i)])
    Z1 <- as.matrix(Z[,i])
    Z0 <- as.matrix(Z[,-c(1,i)])
    X1 <- as.matrix(X[,i])
    X0 <- as.matrix(X[,-c(1,i)])
    Z1_tr <- Z1 # what does this stand for: tr??
    Z0_tr <- Z0 # pre-treatment outcomes?
    Z1_te <- as.matrix(Y1[-(1:T0),]) # what does this stand for: te??
    Z0_te <- as.matrix(Y0[-(1:T0),]) # post treatment outcomes?
    
    # Fit elastic net -> on pre-treatment outcomes
    V1 <- scale(Z1_tr, scale = FALSE)
    V0 <- scale(Z0_tr, scale = FALSE)
    fit <- glmnet(x = V0, y = V1,
                  alpha = a,
                  lambda = lambda_grid,
                  standardize = FALSE,
                  intercept = FALSE) # elastic net function: alpha = a and for all lambda points in lambda grid
    w <- as.matrix(coef(fit, s = lambda_grid)) # Save coefficients of fit for weight
    w <- w[-1,] # Delete intercept
    int <- t(as.matrix(apply(Z1_tr[,rep(1, nlambda)] - Z0_tr %*% w, 2, mean))) # For each lambda point, on pre-treatment outcome
    e <- Z1_te[,rep(1, nlambda)] - int[rep(1, T1),] - Z0_te %*% w # Dimensions?: Columns: lambdas, Rows: T1 (periods after intervention)
    err[i - 1,] <- colMeans(e ^ 2) # SSR for this error term for each lambda point
  }
  
  # Optimal lambda
  #err <- apply(t(scale(t(err))), 2, mean)
  err <- apply(err, 2, mean) # mean error over all control units
  #ind_opt <- max(which(err <= min(err) + 1 * sd(err)))
  ind_opt <- which.min(err) # Find lambda that minimizes error
  err_alpha[j] <- err[ind_opt] # Save that error for this alpha
  lambda_opt_alpha[j] <- lambda_grid[ind_opt] # Save the corresponding lambda value
}

# After loop:
# Optimal a
ind_opt <- which.min(err_alpha) # Find alpha that minimizes error
a_opt <- a_grid[ind_opt]
lambda_opt <- lambda_opt_alpha[ind_opt] # Find associated lambda value

###################################################################
####### Standard Errors without changing alpha and lambda #########
###################################################################

# Over units
std_err_i <- matrix(0, N - 1, T1) # Storage matrix
coef_units <- matrix(0, N-1, N)

for (i in 2:N) {
  cat('(Std. Error) Over Unit i =', toString(i), '\n')
  
  # Define matrices appropriately
  Y1 <- as.matrix(Y[,i])
  Y0 <- as.matrix(Y[,-c(1,i)]) # also delete West Germany here
  Z1 <- as.matrix(Z[,i])
  Z0 <- as.matrix(Z[,-c(1,i)])
  X1 <- as.matrix(X[,i])
  X0 <- as.matrix(X[,-c(1,i)])
  
  # Fit elastic net
  # Use same optimal alpha and lambda here
  V1 <- scale(Z1, scale = FALSE)
  V0 <- scale(Z0, scale = FALSE)
  fit <- glmnet(x = V0, y = V1,
                alpha = a_opt,
                lambda = lambda_grid,
                standardize = FALSE,
                intercept = FALSE) # Fit for optimal alpha and lambda grid
  w <- as.matrix(coef(fit, s = lambda_opt)) # Save coefficients for optimal lambda only
  w <- w[-1,]
  int <- as.matrix(apply(Z1 - Z0 %*% w, 2, mean))
  if(i ==2){
    coef_units[i-1,] <- c(int, 0,w)
  }else{
    if(i==N){
      coef_units[i-1,] <- c(int, w,0)
    }else{
      coef_units[i-1,] <- c(int, w[1:(i-2)],0,w[(i-1):length(w)])
    }
  }
  std_err_i[i - 1,] <- (Y1[-c(1:T0),] - int[rep(1,T1),] - Y0[-c(1:T0),] %*% w) ^ 2 # Save SSR for each unit
}
std_err_i <- as.matrix(sqrt(apply(std_err_i, 2, mean))) # Sqrt of mean of all SSR is Std Err over units

# Y fit
Y_fit <- matrix(0,N-1,14)
for(i in 2:N){
  Y0 <- as.matrix(Y[,-1])
  Y_fit[i-1,] = t(coef_units[i-1,1]+Y0[-c(1:T0),] %*% coef_units[i-1,-1])
}

# Over time
s <- floor(T0 / 2) # Number of periods for Std Errs
std_err_t <- matrix(0, s, 1) # Storage matrix 
coef_time <- matrix(0, s, N)


# Fix matrices for Y and X for time-varying standard errors
Y1 <- as.matrix(Y[,1])
Y0 <- as.matrix(Y[,-1])
X1 <- as.matrix(X[,1])
X0 <- as.matrix(X[,-1])

for (t in 1:s) {
  cat('(Std. Error) Over Time t =', toString(t), '\n')
  
  # Have varying Z matrices: Z = pre-treatment outcome over which MSPE to be minimized
  Z1 <- as.matrix(Z[c(1:(T0 - t)),1]) # Pretend intervention happens earlier
  Z0 <- as.matrix(Z[c(1:(T0 - t)),-1])
  
  # Fit elast: same optimal alpha and lambda
  V1 <- scale(Z1, scale = FALSE)
  V0 <- scale(Z0, scale = FALSE)
  fit <- glmnet(x = V0, y = V1,
                alpha = a_opt,
                lambda = lambda_grid,
                standardize = FALSE,
                intercept = FALSE) # Fit with optimal alpha and over lambda grid
  w <- as.matrix(coef(fit, s = lambda_opt)) # Save for optimal lambda only
  w <- w[-1,] # Delete intercept
  int <- as.matrix(apply(Z1 - Z0 %*% w, 2, mean))
  coef_time[t,] = c(int, w)
  
  std_err_t[t,1] <- (Y1[T0 - t + 1,] - int - Y0[T0 - t + 1,] %*% w) ^ 2 # Save SSR for each time point
}
std_err_t <- as.matrix(sqrt(apply(std_err_t, 2, mean)))

# Y fit
Y_fit_time <- rep(0,s)
for(t in 1:s){
  Y_fit_time[t] = coef_time[t,1] + Y0[T0 - t + 1,] %*% coef_time[t,-1]
}

# Over units and time
std_err_it <- matrix(0, N - 1, 1)

for (i in 2:N) { # units
  
  # Fix Y and X matrices
  Y1 <- as.matrix(Y[,i])
  Y0 <- as.matrix(Y[,-c(1,i)])
  X1 <- as.matrix(X[,i])
  X0 <- as.matrix(X[,-c(1,i)])
  
  std_err_temp <- matrix(0, s, 1) # Std Err for each number of varying periods
  
  for (t in 1:s) {
    cat('(Std. Error) Over Unit and Time ( i , t ) = (',toString(i), ',', 
        toString(t), ')\n')
    
    # Same time-varying Z, but now different unit treatment (i)
    Z1 <- as.matrix(Z[c(1:(T0 - t)),i])
    Z0 <- as.matrix(Z[c(1:(T0 - t)),-c(1,i)])
    
    # Fit elastic net
    V1 <- scale(Z1, scale = FALSE)
    V0 <- scale(Z0, scale = FALSE)
    fit <- glmnet(x = V0, y = V1,
                  alpha = a_opt,
                  lambda = lambda_grid,
                  standardize = FALSE,
                  intercept = FALSE) # Fit function same optimal alpha and lambda grid
    w <- as.matrix(coef(fit, s = lambda_opt)) # Save coefficient only for optimal lambda
    w <- w[-1,] # Delete intercept
    int <- as.matrix(apply(Z1 - Z0 %*% w, 2, mean))
    std_err_temp[t,1] <- (Y1[T0 - t + 1,] - int - Y0[T0 - t + 1,] %*% w) ^ 2 # Save SSR for each time point
  }
  std_err_temp <- as.matrix(apply(std_err_temp, 2, mean)) # Get average std err over time
  std_err_it[i - 1,1] <- std_err_temp # Save that as std err for that specific unit
}
std_err_it <- as.matrix(sqrt(apply(std_err_it, 2, mean))) # Sqrt of mean of these unit time std err is std err

# Copy the standard errors
std_err_elast_i <- std_err_i
std_err_elast_t <- std_err_t
std_err_elast_it <- std_err_it

########################################################################
####### Standard Errors with changing optimal alpha and lambda #########
########################################################################

# Over units
std_err_i_c <- matrix(0, N - 1, T1) # Storage matrix
coef_units <- matrix(0, N-1, N)
optimal_parameters <- matrix(0, N-1, 2)
for (i in 2:N) {
  cat('(Std. Error) Over Unit i =', toString(i), '\n')
  
  # Define matrices appropriately
  Y1 <- as.matrix(Y[,i])
  Y0 <- as.matrix(Y[,-c(1,i)]) # also delete West Germany here
  Z1 <- as.matrix(Z[,i])
  Z0 <- as.matrix(Z[,-c(1,i)])
  X1 <- as.matrix(X[,i])
  X0 <- as.matrix(X[,-c(1,i)])
  
  # Find optimal alpha and lambda in this iteration
  err_alpha <- matrix(0, nrow = na, ncol = 1) # Matrix to store error terms associated with each value of alpha
  lambda_opt_alpha <- matrix(0, nrow = na, ncol = 1) # Matrix to store optimal lambda associated with each alpha
  
  # new Y and Z matrix to make things easier
  Y_temp <- cbind(Y1,Y0)
  Z_temp <- cbind(Z1,Z0)
  X_temp <- cbind(X1,X0)
  
  # Iteration
  for (j in 1:na) { # Iterate over alpha points
    a <- a_grid[j]
    err <- matrix(0, nrow = N - 2, ncol = nlambda) # Matrix for storage of error terms for each control unit and all lambda values
    for (l in 2:(N-1)) { # iterate over units
      
      # Determine matrices appropriately
      Y1 <- as.matrix(Y_temp[,l])
      Y0 <- as.matrix(Y_temp[,-c(1,l)])
      Z1 <- as.matrix(Z_temp[,l])
      Z0 <- as.matrix(Z_temp[,-c(1,l)])
      X1 <- as.matrix(X_temp[,l])
      X0 <- as.matrix(X_temp[,-c(1,l)])
      Z1_tr <- Z1 # what does this stand for: tr??
      Z0_tr <- Z0 # pre-treatment outcomes?
      Z1_te <- as.matrix(Y1[-(1:T0),]) # what does this stand for: te??
      Z0_te <- as.matrix(Y0[-(1:T0),]) # post treatment outcomes?
      
      # Fit elastic net -> on pre-treatment outcomes
      V1 <- scale(Z1_tr, scale = FALSE)
      V0 <- scale(Z0_tr, scale = FALSE)
      fit <- glmnet(x = V0, y = V1,
                    alpha = a,
                    lambda = lambda_grid,
                    standardize = FALSE,
                    intercept = FALSE) # elastic net function: alpha = a and for all lambda points in lambda grid
      w <- as.matrix(coef(fit, s = lambda_grid)) # Save coefficients of fit for weight
      w <- w[-1,] # Delete intercept
      int <- t(as.matrix(apply(Z1_tr[,rep(1, nlambda)] - Z0_tr %*% w, 2, mean))) # For each lambda point, on pre-treatment outcome
      e <- Z1_te[,rep(1, nlambda)] - int[rep(1, T1),] - Z0_te %*% w # Dimensions?: Columns: lambdas, Rows: T1 (periods after intervention)
      err[l - 1,] <- colMeans(e ^ 2) # SSR for this error term for each lambda point
    }
    
    # Optimal lambda
    #err <- apply(t(scale(t(err))), 2, mean)
    err <- apply(err, 2, mean) # mean error over all control units
    #ind_opt <- max(which(err <= min(err) + 1 * sd(err)))
    ind_opt <- which.min(err) # Find lambda that minimizes error
    err_alpha[j] <- err[ind_opt] # Save that error for this alpha
    lambda_opt_alpha[j] <- lambda_grid[ind_opt] # Save the corresponding lambda value
  }
  
  # After loop:
  # Optimal a
  ind_opt <- which.min(err_alpha) # Find alpha that minimizes error
  a_opt_temp <- a_grid[ind_opt]
  lambda_opt_temp <- lambda_opt_alpha[ind_opt] # Find associated lambda value
  
  optimal_parameters[i-1,] <- c(a_opt_temp, lambda_opt_temp) # Save optimal parameters for each unit
  
  # Define matrices appropriately
  Y1 <- as.matrix(Y[,i])
  Y0 <- as.matrix(Y[,-c(1,i)]) # also delete West Germany here
  Z1 <- as.matrix(Z[,i])
  Z0 <- as.matrix(Z[,-c(1,i)])
  X1 <- as.matrix(X[,i])
  X0 <- as.matrix(X[,-c(1,i)])
  
  # Fit elastic net
  # Use new optimal alpha and lambda here
  V1 <- scale(Z1, scale = FALSE)
  V0 <- scale(Z0, scale = FALSE)
  fit <- glmnet(x = V0, y = V1,
                alpha = a_opt_temp,
                lambda = lambda_grid,
                standardize = FALSE,
                intercept = FALSE) # Fit for optimal alpha and lambda grid
  w <- as.matrix(coef(fit, s = lambda_opt_temp)) # Save coefficients for optimal lambda only
  w <- w[-1,]
  int <- as.matrix(apply(Z1 - Z0 %*% w, 2, mean))
  if(i ==2){
    coef_units[i-1,] <- c(int, 0,w)
  }else{
    if(i==N){
      coef_units[i-1,] <- c(int, w,0)
    }else{
      coef_units[i-1,] <- c(int, w[1:(i-2)],0,w[(i-1):length(w)])
    }
  }
  std_err_i_c[i - 1,] <- (Y1[-c(1:T0),] - int[rep(1,T1),] - Y0[-c(1:T0),] %*% w) ^ 2 # Save SSR for each unit
}
std_err_i_c <- as.matrix(sqrt(apply(std_err_i_c, 2, mean))) # Sqrt of mean of all SSR is Std Err over units

# Y fit
Y_fit <- matrix(0,N-1,14)
for(i in 2:N){
  Y0 <- as.matrix(Y[,-1])
  Y_fit[i-1,] = t(coef_units[i-1,1]+Y0[-c(1:T0),] %*% coef_units[i-1,-1])
}

# Over time
s <- floor(T0 / 2) # Number of periods for Std Errs
std_err_t_c <- matrix(0, s, 1) # Storage matrix 
optimal_parameters_t <- matrix(0, s, 2)
coef_time <- matrix(0, s, N)

# Fix matrices for Y and X for time-varying standard errors
Y1 <- as.matrix(Y[,1])
Y0 <- as.matrix(Y[,-1])
X1 <- as.matrix(X[,1])
X0 <- as.matrix(X[,-1])

for (t in 1:s) {
  cat('(Std. Error) Over Time t =', toString(t), '\n')
  
  # Have varying Z matrices: Z = pre-treatment outcome over which MSPE to be minimized
  Z1 <- as.matrix(Z[c(1:(T0 - t)),1]) # Pretend intervention happens earlier
  Z0 <- as.matrix(Z[c(1:(T0 - t)),-1])
  
  # Find optimal alpha and lambda in this iteration
  err_alpha <- matrix(0, nrow = na, ncol = 1) # Matrix to store error terms associated with each value of alpha
  lambda_opt_alpha <- matrix(0, nrow = na, ncol = 1) # Matrix to store optimal lambda associated with each alpha
  
  # new Y and Z matrix to make things easier
  Z_temp <- cbind(Z1,Z0)

  # Iteration
  for (j in 1:na) { # Iterate over alpha points
    a <- a_grid[j]
    err <- matrix(0, nrow = N - 1, ncol = nlambda) # Matrix for storage of error terms for each control unit and all lambda values
    for (l in 2:N) { # iterate over units
      
      # Determine matrices appropriately
      Y1 <- as.matrix(Y[,l])
      Y0 <- as.matrix(Y[,-c(1,l)])
      Z1 <- as.matrix(Z_temp[,l])
      Z0 <- as.matrix(Z_temp[,-c(1,l)])
      X1 <- as.matrix(X[,l])
      X0 <- as.matrix(X[,-c(1,l)])
      Z1_tr <- Z1 # what does this stand for: tr??
      Z0_tr <- Z0 # pre-treatment outcomes?
      Z1_te <- as.matrix(Y1[-(1:(T0 - t)),]) # what does this stand for: te??
      Z0_te <- as.matrix(Y0[-(1:(T0 - t)),]) # post treatment outcomes?
      
      # Fit elastic net -> on pre-treatment outcomes
      V1 <- scale(Z1_tr, scale = FALSE)
      V0 <- scale(Z0_tr, scale = FALSE)
      fit <- glmnet(x = V0, y = V1,
                    alpha = a,
                    lambda = lambda_grid,
                    standardize = FALSE,
                    intercept = FALSE) # elastic net function: alpha = a and for all lambda points in lambda grid
      w <- as.matrix(coef(fit, s = lambda_grid)) # Save coefficients of fit for weight
      w <- w[-1,] # Delete intercept
      int <- t(as.matrix(apply(Z1_tr[,rep(1, nlambda)] - Z0_tr %*% w, 2, mean))) # For each lambda point, on pre-treatment outcome
      e <- Z1_te[,rep(1, nlambda)] - int[rep(1, T1+t),] - Z0_te %*% w # Dimensions?: Columns: lambdas, Rows: T1 (periods after intervention)
      err[l - 1,] <- colMeans(e ^ 2) # SSR for this error term for each lambda point
    }
    
    # Optimal lambda
    #err <- apply(t(scale(t(err))), 2, mean)
    err <- apply(err, 2, mean) # mean error over all control units
    #ind_opt <- max(which(err <= min(err) + 1 * sd(err)))
    ind_opt <- which.min(err) # Find lambda that minimizes error
    err_alpha[j] <- err[ind_opt] # Save that error for this alpha
    lambda_opt_alpha[j] <- lambda_grid[ind_opt] # Save the corresponding lambda value
  }
  
  # After loop:
  # Optimal a
  ind_opt <- which.min(err_alpha) # Find alpha that minimizes error
  a_opt_temp <- a_grid[ind_opt]
  lambda_opt_temp <- lambda_opt_alpha[ind_opt] # Find associated lambda value
  
  optimal_parameters_t[t,] <- c(a_opt_temp, lambda_opt_temp) # Save optimal parameters for each unit
  
  # Fix matrices for Y and X for time-varying standard errors
  Y1 <- as.matrix(Y[,1])
  Y0 <- as.matrix(Y[,-1])
  X1 <- as.matrix(X[,1])
  X0 <- as.matrix(X[,-1])
  
  # Have varying Z matrices: Z = pre-treatment outcome over which MSPE to be minimized
  Z1 <- as.matrix(Z[c(1:(T0 - t)),1]) # Pretend intervention happens earlier
  Z0 <- as.matrix(Z[c(1:(T0 - t)),-1])
  
  
  # Fit elast: same optimal alpha and lambda
  V1 <- scale(Z1, scale = FALSE)
  V0 <- scale(Z0, scale = FALSE)
  fit <- glmnet(x = V0, y = V1,
                alpha = a_opt_temp,
                lambda = lambda_grid,
                standardize = FALSE,
                intercept = FALSE) # Fit with optimal alpha and over lambda grid
  w <- as.matrix(coef(fit, s = lambda_opt_temp)) # Save for optimal lambda only
  w <- w[-1,] # Delete intercept
  int <- as.matrix(apply(Z1 - Z0 %*% w, 2, mean))
  coef_time[t,] = c(int, w)
  std_err_t_c[t,1] <- (Y1[T0 - t + 1,] - int - Y0[T0 - t + 1,] %*% w) ^ 2 # Save SSR for each time point
}
std_err_t_c <- as.matrix(sqrt(apply(std_err_t_c, 2, mean)))

# Y fit
Y_fit_time <- rep(0,s)
for(t in 1:s){
  Y_fit_time[t] = coef_time[t,1] + Y0[T0 - t + 1,] %*% coef_time[t,-1]
}

# Over units and time
std_err_it_c <- matrix(0, N - 1, 1)
optimal_parameters_it <- array(0, dim= c(N-1, s, 2))

for (i in 2:N) { # units
  
  # Fix Y and X matrices
  Y1 <- as.matrix(Y[,i])
  Y0 <- as.matrix(Y[,-c(1,i)])
  X1 <- as.matrix(X[,i])
  X0 <- as.matrix(X[,-c(1,i)])
  
  # new Y and Z matrix to make things easier
  Y_temp <- cbind(Y1,Y0)
  X_temp <- cbind(X1,X0)
  
  std_err_temp <- matrix(0, s, 1) # Std Err for each number of varying periods
  
  for (t in 1:s) {
    cat('(Std. Error) Over Unit and Time ( i , t ) = (',toString(i), ',', 
        toString(t), ')\n')
    
    # Same time-varying Z, but now different unit treatment (i)
    Z1 <- as.matrix(Z[c(1:(T0 - t)),i])
    Z0 <- as.matrix(Z[c(1:(T0 - t)),-c(1,i)])
    
    # Find optimal alpha and lambda in this iteration
    err_alpha <- matrix(0, nrow = na, ncol = 1) # Matrix to store error terms associated with each value of alpha
    lambda_opt_alpha <- matrix(0, nrow = na, ncol = 1) # Matrix to store optimal lambda associated with each alpha
    
    # new Y and Z matrix to make things easier
    Z_temp <- cbind(Z1,Z0)
    
    # Iteration
    for (j in 1:na) { # Iterate over alpha points
      a <- a_grid[j]
      err <- matrix(0, nrow = N - 2, ncol = nlambda) # Matrix for storage of error terms for each control unit and all lambda values
      for (l in 2:(N-1)) { # iterate over units
        
        # Determine matrices appropriately
        Y1 <- as.matrix(Y_temp[,l])
        Y0 <- as.matrix(Y_temp[,-c(1,l)])
        Z1 <- as.matrix(Z_temp[,l])
        Z0 <- as.matrix(Z_temp[,-c(1,l)])
        X1 <- as.matrix(X_temp[,l])
        X0 <- as.matrix(X_temp[,-c(1,l)])
        Z1_tr <- Z1 # what does this stand for: tr??
        Z0_tr <- Z0 # pre-treatment outcomes?
        Z1_te <- as.matrix(Y1[-(1:(T0 - t)),]) # what does this stand for: te??
        Z0_te <- as.matrix(Y0[-(1:(T0 - t)),]) # post treatment outcomes?
        
        # Fit elastic net -> on pre-treatment outcomes
        V1 <- scale(Z1_tr, scale = FALSE)
        V0 <- scale(Z0_tr, scale = FALSE)
        fit <- glmnet(x = V0, y = V1,
                      alpha = a,
                      lambda = lambda_grid,
                      standardize = FALSE,
                      intercept = FALSE) # elastic net function: alpha = a and for all lambda points in lambda grid
        w <- as.matrix(coef(fit, s = lambda_grid)) # Save coefficients of fit for weight
        w <- w[-1,] # Delete intercept
        int <- t(as.matrix(apply(Z1_tr[,rep(1, nlambda)] - Z0_tr %*% w, 2, mean))) # For each lambda point, on pre-treatment outcome
        e <- Z1_te[,rep(1, nlambda)] - int[rep(1, T1+t),] - Z0_te %*% w # Dimensions?: Columns: lambdas, Rows: T1 (periods after intervention)
        err[l - 1,] <- colMeans(e ^ 2) # SSR for this error term for each lambda point
      }
      
      # Optimal lambda
      #err <- apply(t(scale(t(err))), 2, mean)
      err <- apply(err, 2, mean) # mean error over all control units
      #ind_opt <- max(which(err <= min(err) + 1 * sd(err)))
      ind_opt <- which.min(err) # Find lambda that minimizes error
      err_alpha[j] <- err[ind_opt] # Save that error for this alpha
      lambda_opt_alpha[j] <- lambda_grid[ind_opt] # Save the corresponding lambda value
    }
    
    # After loop:
    # Optimal a
    ind_opt <- which.min(err_alpha) # Find alpha that minimizes error
    a_opt_temp <- a_grid[ind_opt]
    lambda_opt_temp <- lambda_opt_alpha[ind_opt] # Find associated lambda value
    
    optimal_parameters_it[i-1,t,] <- c(a_opt_temp, lambda_opt_temp) # Save optimal parameters for each unit
    
    # Same time-varying Z, but now different unit treatment (i)
    Z1 <- as.matrix(Z[c(1:(T0 - t)),i])
    Z0 <- as.matrix(Z[c(1:(T0 - t)),-c(1,i)])
    
    # Fit elastic net
    V1 <- scale(Z1, scale = FALSE)
    V0 <- scale(Z0, scale = FALSE)
    fit <- glmnet(x = V0, y = V1,
                  alpha = a_opt_temp,
                  lambda = lambda_grid,
                  standardize = FALSE,
                  intercept = FALSE) # Fit function same optimal alpha and lambda grid
    w <- as.matrix(coef(fit, s = lambda_opt_temp)) # Save coefficient only for optimal lambda
    w <- w[-1,] # Delete intercept
    int <- as.matrix(apply(Z1 - Z0 %*% w, 2, mean))
    std_err_temp[t,1] <- (Y1[T0 - t + 1,] - int - Y0[T0 - t + 1,] %*% w) ^ 2 # Save SSR for each time point
  }
  std_err_temp <- as.matrix(apply(std_err_temp, 2, mean)) # Get average std err over time
  std_err_it_c[i - 1,1] <- std_err_temp # Save that as std err for that specific unit
}
std_err_it_c <- as.matrix(sqrt(apply(std_err_it_c, 2, mean))) # Sqrt of mean of these unit time std err is std err

# Copy the standard errors
std_err_elast_i_c <- std_err_i_c
std_err_elast_t_c <- std_err_t_c
std_err_elast_it_c <- std_err_it_c
