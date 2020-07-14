#######################################
########### FUNCTION SCRIPT ###########
#######################################

## Libraries
library(lars)
library(glmnet)
library(foreign)
library(reshape2)
library(R.matlab)
library(Synth)
library(modopt.matlab)


# Function to get weights, fitted values and standard errors

general_estimate <- function(data, method = NULL, lambda_grid, alpha_grid, ind_treatment =1){
  
  ## INPUT:
  #
  # data:
  # method:
  # lambda_grid: pre-specified lambda grid over which we are optimizing
  # alpha_grid: pre-specified alpha grid over which we are optimizing
  # ind_treatment: indicator which unit is the treatment unit; default is column 1
  
  ## OUTPUT:
  #
  # int: estimated intercept
  # w: estimated weights
  # Y_est: fitted Y values
  # Y_true: true Y-values
  # alpha_opt: optimal alpha (if method is elastic net)
  # lambda_opt: optimal lambda (if method is elastic net)
  # std_err_i: standard error over units
  # std_err_t: standard error over time
  # std_err_it: standard error over unit and time
  # T_0: time point of intervention
  # T_1: time points after intervention
  
  # Make sure that one of the methods that is feasible is specified
  if(is.null(method) | method != "diff_in_diff" | method != "elastic_net" | method != "constr_reg" | method != "synth" | method != "best_subset"){
    stop('Please specify one of the following methods: "diff_in_diff", "elastic_net", "constr_reg", "synth" or "best_subset"!')
  }else{
    ####
    # Dataprep
    # Here or before?
    
    ####
    # Get weights
    
    # Elastic Net: Find tuning parameter and weights
    if(method == "elastic_net"){
      params <- tuning_parameters_elastic_net(Y,Z,X,lambda_grid, alpha_grid, ind_treatment)
      w <- find_weights_elastic_net(Y, Z, X, params$alpha, params$lambda, lambda_grid, ind_treatment) 
    }
    
    # Best subset: Find weights
    if(method == "best_subset"){
      w <- find_weights_subset(Y,Z,X, ind_treatment)
    }
    
    # Constrained regression: Find weights
    if(method == "constr_reg"){
      w <- find_weights_constr_reg(Y,Z,X, ind_treatment)
    }
    
    ###
    # Get estimate
    Y_est = w$intercept + Y0%*%w$weights
    
    ###
    # Get standard error
    
    ###
    # Output
    if(method == "elastic_net"){
      out <- list("int" = w$intercept, "w" = w$weights, "Y_est" = Y_est, "Y_true" = Y_true, "alpha_opt" = params$alpha, "lambda_opt" = params$lambda,"std_err_i" = std_err_i, "std_err_t" = std_err_t, "std_err_it" = std_err_it, "T_0"= T0,"T_1"=T1)
    }else{
      out <- list("int" = w$intercept, "w" = w$weights, "Y_est" = Y_est, "Y_true" = Y_true,"std_err_i" = std_err_i, "std_err_t" = std_err_t, "std_err_it" = std_err_it, "T_0"= T0,"T_1"=T1)
    }
  }
}


###############################
##### Auxiliary functions #####
###############################

tuning_parameters_elastic_net <- function(Y,Z,X,lambda_grid, alpha_grid, ind_treatment){
  
  ## INPUT:
  #
  # Y: matrix of outcome variables for treated unit and controls 
  # Z: matrix of pre-treatment outcomes for treated unit and controls 
  # X: matrix of covariates for treated unit and controls 
  # lambda_grid: pre-specified lambda grid over which we are optimizing
  # alpha_grid: pre-specified alpha grid over which we are optimizing
  # ind_treatment: indicator which unit is the treatment unit
  
  ## OUTPUT:
  #
  # alpha: optimal alpha
  # lambda: optimal lambda
  
  # Parameters
  N <- dim(Y)[2] # Number of units
  T <- dim(Y)[1] # Number of time periods
  T0 <- dim(Z)[1] # Time of intervention
  T1 <- T - T0 # Number of time periods after intervention
  nlambda <- length(lambda_grid) # length of lambda grid
  na <- length(alpha_grid) # length of alpha grid
  err_alpha <- matrix(0, nrow = na, ncol = 1) # Matrix to store error terms associated with each value of alpha
  lambda_opt_alpha <- matrix(0, nrow = na, ncol = 1) # Matrix to store optimal lambda associated with each alpha
  
  # Rearrange matrices to have treated unit first
  Y <- as.matrix(cbind(Y[,ind_treatment], Y[,-ind_treatment]))
  Z <- as.matrix(cbind(Z[,ind_treatment], Z[,-ind_treatment]))
  X <- as.matrix(cbind(X[,ind_treatment], X[,-ind_treatment]))
  
  # LOOP:
  for (j in 1:na) { # Iterate over alpha points
    a <- alpha_grid[j]
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
    err <- apply(err, 2, mean) # mean error over all control units
    ind_opt <- which.min(err) # Find lambda that minimizes error
    err_alpha[j] <- err[ind_opt] # Save that error for this alpha
    lambda_opt_alpha[j] <- lambda_grid[ind_opt] # Save the corresponding lambda value
  }
  
  
  # Determine optimal alpha and lambda
  ind_opt <- which.min(err_alpha) # Find alpha that minimizes error
  alpha_opt <- alpha_grid[ind_opt]
  lambda_opt <- lambda_opt_alpha[ind_opt] # Find associated lambda value
  
  # OUTPUT
  out <- list("alpha"= alpha_opt, "lambda"= lambda_opt)
}

find_weights_elastic_net <- function(Y, Z, X, alpha, lambda, lambda_grid, ind_treatment){
  
  ## INPUT:
  #
  # Y: matrix of outcome variables for treated unit and controls 
  # Z: matrix of pre-treatment outcomes for treated unit and controls 
  # X: matrix of covariates for treated unit and controls 
  # alpha: optimal value for elastic fit alpha
  # lambda: optimal value for elastic fit lambda
  # lambda_grid: lambda_grid used to find optimal lambda -> faster fit if whole sequence of lambdas is fitted instead of single one
  # ind_treatment: indicator which unit is the treatment unit
  
  ## OUTPUT:
  #
  # intercept: intercept resulting from elastic net fit
  # weights: weights resulting from elastic net fit
  
  # Parameters
  N <- dim(Y)[2] # Number of units
  
  # Normalize predictors
  div <- as.matrix(apply(X, 1, sd)) # Matrix of standard deviations for each predictor
  X <- X / div[,rep(1, N)] # Standardizes each predictor to have std 1
  ## SKIP THIS IF NOT NEEDED!
  
  # Matrices for storage
  int <- matrix(0, nrow = 1, ncol = 1)
  w <- matrix(0, nrow = N - 1, ncol = 1)
  
  # Fix treatment unit
  i <- ind_treatment
  Y1 <- as.matrix(Y[,i])
  Y0 <- as.matrix(Y[,-i])
  Z1 <- as.matrix(Z[,i])
  Z0 <- as.matrix(Z[,-i])
  X1 <- as.matrix(X[,i])
  X0 <- as.matrix(X[,-i])
  V1 <- scale(Z1, scale = FALSE)
  V0 <- scale(Z0, scale = FALSE)
  
  # Fit elastic net
  fit <- glmnet(x = V0, y = V1,
                alpha = alpha,
                lambda = lambda_grid,
                standardize = FALSE, # Don't standardize
                intercept = FALSE) # Function to fit elastic net with optimal alpha over lambda grid
  w <- as.matrix(coef(fit, s = lambda)) # only save coefficients for optimal lambda
  w <- w[-1,] # Delete intercept value (0.0000); elastic net weights for all 16 control units
  int <- as.matrix(apply(Z1 - Z0 %*% w, 2, mean)) # Intercept for elastic net fit
  
  # Output
  out <- list("intercept" = int, "weights" =w)
  
}

find_weights_subset <- function(Y,Z,X,ind_treatment){
  
  ## INPUT:
  #
  # Y: matrix of outcome variables for treated unit and controls 
  # Z: matrix of pre-treatment outcomes for treated unit and controls 
  # X: matrix of covariates for treated unit and controls 
  # ind_treatment: indicator which unit is the treatment unit
  
  ## OUTPUT:
  #
  # intercept: intercept resulting from best subset selection fit
  # weights: weights resulting from best subset selection fit
  
  # Parameters
  N <- dim(Y)[2] # Number of units
  T <- dim(Y)[1] # Number of time periods
  T0 <- dim(Z)[1] # Time of intervention
  T1 <- T - T0 # Number of time periods after intervention
  T0_tr <- floor(T0 * 2 / 3) # what is this?? - here it is actually used
  
  
  # Normalize predictors
  div <- as.matrix(apply(X, 1, sd)) # Matrix of standard deviations for each predictor
  X <- X / div[,rep(1, N)] # Standardizes each predictor to have std 1
  
  # Matrices for storage
  int <- matrix(0, nrow = 1, ncol = 1)
  w <- matrix(0, nrow = N - 1, ncol = 1)
  
  #########################################
  ####### OPTIMAL UNITS IN SUBSET #########
  #########################################
  
  # Rearrange matrices to have treated unit first
  Y <- as.matrix(cbind(Y[,ind_treatment], Y[,-ind_treatment]))
  Z <- as.matrix(cbind(Z[,ind_treatment], Z[,-ind_treatment]))
  X <- as.matrix(cbind(X[,ind_treatment], X[,-ind_treatment]))
  
  ## Find the optimal subset #!!
  # Iterate over i
  
  # Still figure out what is going on here!!
  n_max <- N - 1 # Number of units in subset ?
  n_grid <- c(0:min(T0_tr - 1, n_max, N - 2)) # What is this grid?? Why N-2 and N-1?
  nn <- length(n_grid) # Number of points in n grid
  err <- matrix(0, nrow = N - 1, ncol = nn) # Storage for errors for each unit and n
  c <- matrix(1, nrow = T0, ncol = 1) # Just a constant! Needed for fit here later
  
  for (i in 2:N) { # over units
    
    # Fix matrices 
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
    
    for (n in 0:(nn - 1)) { # iterate over n's
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
      # Choose the optimal subset of size n and compute the error over time
      j_opt <- which.min(err_cur)
      e <- Z1_te - int[rep(1, T1),j_opt] - Z0_te %*% w[,j_opt]
      err[i - 1,n + 1] <- mean(e ^ 2)
    }
  }
  
  # Determine optimal n
  err <- apply(t(scale(t(err))), 2, mean) # Get average error over all units
  ind_opt <- which.min(err)
  n_opt <- n_grid[ind_opt]
  
  ################################
  ####### OPTIMAL SUBSET #########
  ################################
  
  # Matrices for storage
  int_subs <- matrix(0, nrow = 1, ncol = 1)
  w_subs <- matrix(0, nrow = N - 1, ncol = 1)
  Y_subs <- matrix(0, nrow = T, ncol = 1)
  Y_true <- matrix(0, nrow = T, ncol = 1)
  c <- matrix(1, nrow = T0, ncol = 1)
  
  # Fix treatment unit
  i <- 1
  Y1 <- as.matrix(Y[,i])
  Y0 <- as.matrix(Y[,-i])
  Z1 <- as.matrix(Z[,i])
  Z0 <- as.matrix(Z[,-i])
  X1 <- as.matrix(X[,i])
  X0 <- as.matrix(X[,-i])
  
  ## Loop over each possible combination
  V1 <- Z1
  subs_n <- combn(c(1:(N - 1)), n_opt, simplify = FALSE) # Get all the subset combinations of size n*
  int_s <- matrix(0, nrow = 1, ncol = length(subs_n))
  w_s <- matrix(0, nrow = N - 1, ncol = length(subs_n))
  err_cur <- matrix(0, nrow = length(subs_n), ncol = 1)
  
  for (j in 1:length(subs_n)) {
    sub <- subs_n[[j]]
    V0 <- cbind(c, Z0[,sub])
    w_cur <- solve(t(V0) %*% V0, t(V0) %*% V1)
    int_s[1,j] <- w_cur[1]
    w_s[sub,j] <- w_cur[-1]
    err_cur[j,1] <- mean((V1 - V0 %*% w_cur) ^ 2)
  }
  
  # Determine indicator for which subset is best
  j_opt <- which.min(err_cur)
  
  ###########################
  ####### Fit Model #########
  ###########################
  
  int <- as.matrix(int_s[1,j_opt])
  w <- w_s[,j_opt]
  
  # Output
  out <- list("intercept" = int, "weights" =w)
}

find_weights_constr_reg <- function(Y,Z,X,ind_treatment){
  ## INPUT:
  #
  # Y: matrix of outcome variables for treated unit and controls 
  # Z: matrix of pre-treatment outcomes for treated unit and controls
  # X: matrix of covariates for treated unit and controls 
  # ind_treatment: indicator which unit is the treatment unit
  
  ## OUTPUT:
  #
  # intercept: intercept resulting from constrained regression fit
  # weights: weights resulting from constrained regression fit
  
  # Parameters
  N <- dim(Y)[2] # Number of units
  T <- dim(Y)[1] # Number of time periods
  T0 <- dim(Z)[1] # Time of intervention
  T1 <- T - T0 # Number of time periods after intervention
  
  # Matrices for storage
  int <- matrix(0, nrow = 1, ncol = 1)
  w <- matrix(0, nrow = N - 1, ncol = 1)
  
  # Fix treatment unit
  i <- ind_treatment
  Y1 <- as.matrix(Y[,i])
  Y0 <- as.matrix(Y[,-i])
  Z1 <- as.matrix(Z[,i])
  Z0 <- as.matrix(Z[,-i])
  X1 <- as.matrix(X[,i])
  X0 <- as.matrix(X[,-i])
  V1 <- Z1round
  V0 <- Z0
  
  # Fit constrained regression
  H = t(V0)%*%V0
  f = as.vector(-t(V0)%*%V1)
  
  Aeq = rep(1,N-1)
  beq = 1
  lb = rep(0, N-1)
  ub = rep(1,N-1)
  
  w = quadprog(H, f, NULL, NULL, Aeq, beq, lb, ub)
  
  out <- list("intercept" = 0,"weights"= w$x)
}