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
library(LowRankQP)


# Function to get weights, fitted values and standard errors

#note: I think that ind_treatment shouldn't be defaulted as 1 here, since the only place we use ind_treatment after prep_data is 
#in synthetic control stuff and ind_treatment is going to be not necessarily one.
# Yes, agreed! The way we have it now, it should not be defualted to 1

general_estimate <- function(Y, Z, X, W, method = NULL, tune_params = NULL){
  
  ## INPUT:
  #
  # Y:              matrix of outcome variables for treated unit and controls.
  #                   (Note that currently our code only will work for a single treated unit and many control units)
  # Z:              matrix of pre-treatment outcomes for treated unit and controls 
  # X:              matrix of covariates for treated unit and controls 
  # W:              Vector that specifies which unit is treated
  # method:         "diff_in_diff", "elastic_net", "constr_reg", "synth","best_subset"; 
  # tune_params:    vector of parameters needed for methods that need cross-validation for tuning parameters.
  #                   for elastic_net, it should have a lambda grid and an alpha grid 
  #                   for synth, it should be Y_tune, Z_tune, X_tune matrices that have been previously prepared like X, Y, Z but use the years and predictors 
  #                     for synthetic control cross-validations
  # ind_treatment:  the column that corresponds to the treated unit
  
  ## OUTPUT:
  #
  # int:            estimated intercept
  # w:              estimated weights
  # Y_est:          fitted Y values
  # Y_true:         true Y-values
  # alpha_opt:      optimal alpha (if method is elastic net)
  # lambda_opt:     optimal lambda (if method is elastic net)
  # n_opt:          optimal number of control units (if method is best subset)
  # std_err_i:      standard error over units
  # std_err_t:      standard error over time
  # std_err_it:     standard error over unit and time
  # T_0:            time point of intervention
  # T_1:            time points after intervention
  
  # Make sure that one of the methods that is feasible is specified
  if(is.null(method) || (method != "diff_in_diff" && method != "elastic_net" && method != "constr_reg" && method != "synth" && method != "best_subset")){
    stop('Please specify one of the following methods: "diff_in_diff", "elastic_net", "constr_reg", "synth" or "best_subset"!')
  }else{
    
    
    ################ Find Treated Unit ##############

    if (sum(W!=0>1)){
      stop('More than one unit is specified as treated in W. Please check the assigment matrix.')
    }
    
    ind_treatment <- which.max(colSums(W!=0))
    
    ################ Find weights ################
    
    # Elastic Net
    if(method == "elastic_net"){
      
      #find tuning parameters using cross-validation
      params <- tuning_parameters_elastic_net(Y, Z, X, tune_params[[1]], tune_params[[2]])
      
      #find imputed Y weights using the tuning parameters we found
      
      w <- find_weights_elastic_net(Y, Z, X, params$alpha, params$lambda, tune_params[[1]]) 
      
      #reassign tune_params to plug into the standard error functions
      tune_params <- list(params$alpha, params$lambda, tune_params[[1]])
    }
    
    #synthetic control
    if(method == "synth"){
   
      #find tuning parameters (vweights)
      v <- tuning_parameters_synth(tune_params$Y, tune_params$Z, tune_params$X)
      
      #find synthetic control weights
      w <- find_weights_synth(Y, Z, X, vweight = v)
      
      
    }
    
    # Best subset
    if(method == "best_subset"){
      
      #find the optimal value for n
      n_opt <- tuning_parameters_best_subset(Y,Z,X)
      
      #find the best subset using the n_opt from the previous part
      w <- find_weights_subset(Y,Z,X,n_opt)
      
      #reassign tune_params to plug into standard errors
      tune_params <- n_opt
    }
    
    # Diff-in-diff: 
    if (method == "diff_in_diff"){
      
      #find diff-in-diff weights (simply by using 1/N-1). We use the extra function step for consistency and readability
      w <- find_weights_did(Y, Z, X)
    }
    
    # Constrained regression:
    if(method == "constr_reg"){
      
      #find weights using constrained regression
      w <- find_weights_constr_reg(Y,Z,X)
    }
    
    
    ################ Get estimate ################
    
    Y_est = rep(w$intercept,nrow(Y)) + Y[,-1]%*%w$weights
    Y_true = Y[,1]
    
    
    ################ Find Standard Error ################
    
    #don't need to specify ind_treat since they're always in the 1 slot
    std_err_i = general_se(Y, Z, X, method = method, se_method = "unit", tune_params = tune_params)
    std_err_t = general_se(Y, Z, X, method = method, se_method = "time", tune_params = tune_params)
    std_err_it = general_se(Y, Z, X, method = method, se_method = "unit_time", tune_params = tune_params)
    
    
    ################ Output ################
    
    if(method == "elastic_net"){
      out <- list("int" = w$intercept, "w" = w$weights, "Y_est" = Y_est, "Y_true" = Y_true, "alpha_opt" = params$alpha, "lambda_opt" = params$lambda,"std_err_i" = std_err_i, "std_err_t" = std_err_t, "std_err_it" = std_err_it, "T_0" = T_0, "T_1" = T_1)
    }else{
      if(method == "best_subset"){
        out <- list("int" = w$intercept, "w" = w$weights, "Y_est" = Y_est, "Y_true" = Y_true,"n_opt" = n_opt, "std_err_i" = std_err_i, "std_err_t" = std_err_t, "std_err_it" = std_err_it, "T_0" = T_0, "T_1" = T_1)
      }else{
        out <- list("int" = w$intercept, "w" = w$weights, "Y_est" = Y_est, "Y_true" = Y_true,"std_err_i" = std_err_i, "std_err_t" = std_err_t, "std_err_it" = std_err_it, "T_0" = T_0, "T_1" = T_1)
        
      }
    }
  }
}



###############################
##### Auxiliary functions #####
###############################


#Prepare the data for finding tuning parameters and weights
prep_data <- function(d, ind_treat, pred, dep, u, t, spec, cont_set, years, names){
  
  # INPUT 
  #
  #Note: We could likely write this for dataframes that don't have string identifiers and just have numeric identifiers
  #
  # prep_params[[1]]  = pred: vector containing string with predictor variables. (This becomes X)
  # prep_params[[2]]  = dep: string specifying which one is dependent variable. (This becomes Y and Z)
  # prep_params[[3]]  = u: integer identifying unit variable (which column in data frame specifies index)
  # prep_names[[4]]   = t: integer identifying time variable (which column in data frame specifies time)
  # prep_names[[5]]   = spec: list of special predictors (also becomes X)
  # prep_names[[6]]   = cont_set: vector of the columns for control units.
  # prep_names[[6]]   = subs: vector of the columns for control units.
  # prep_names[[7]]   = years: vector specifying years or another time variable as follows  
  #                             [1]: start for predictor aggregation, [2]: end for predictor aggregation (these determine what years are averaged for for X variables listed in "pred")
  #                             [3]: start of prior period, [4]: end of prior period (should be the treatment year) (these years determine what years are in Z)
  #                             [5]: start time plot, [6]: end time plot  (probably just your whole time period)                                                       
  # prep_names[[8]]   = names: unit names variable (which column in data frame specifies unit names)
  
  # OUTPUT
  #
  # X:                  matrix of non-outcome predictors
  # Y:                  matrix of outcomes
  # Z:                  matrix of prior period outcomes to be used as the main predictors
  
  #run dataprep (from the synth package) to prepare the data
  dataprep.out <-
    dataprep(
      foo = d,
      predictors    = pred,
      dependent     = dep,
      unit.variable = u,
      time.variable = t,
      special.predictors = spec,
      treatment.identifier = ind_treat,
      controls.identifier = cont_set,
      time.predictors.prior = years[1]:years[2],
      time.optimize.ssr = years[3]:years[4],
      unit.names.variable = names,
      time.plot = years[5]:years[6]
    )

  #Extract the treatment and control units from the dataprep object  
  X0 <- dataprep.out$X0
  X1 <- dataprep.out$X1
  
  Z1 <- dataprep.out$Z1
  Z0 <- dataprep.out$Z0
  
  Y1 <- dataprep.out$Y1plot
  Y0 <- dataprep.out$Y0plot
  
  #Re-bind data so that the treated unit is in the first position
  # [X1,X0]
  X <- cbind(X1, X0)
  
  # [Y1,Y0]
  Y <- cbind(Y1, Y0)
  
  # [Z1,Z0]
  Z <- cbind(Z1, Z0)
  
  output <- list( "X" = X, "Y" = Y, "Z" = Z)
}

# Find optimal tuning parameters (alpha and lambda values) for the elastic net fit
tuning_parameters_elastic_net <- function(Y,Z,X,lambda_grid, alpha_grid, ind_treatment=1){
  
  # INPUT:
  #
  # Y:              matrix of outcome variables for treated unit and controls 
  # Z:              matrix of pre-treatment outcomes for treated unit and controls 
  # X:              matrix of covariates for treated unit and controls 
  # lambda_grid:    pre-specified lambda grid over which we are optimizing
  # alpha_grid:     pre-specified alpha grid over which we are optimizing
  # ind_treatment:  indicator which unit is the treatment unit in Y,Z and X matrices
  
  # OUTPUT:
  #
  # alpha:          optimal alpha
  # lambda:         optimal lambda
  
  ################ Find parameters using inputted data ################
  
  N <- dim(Y)[2]  # Number of units
  T <- dim(Y)[1]  # Number of time periods
  T0 <- dim(Z)[1] # Time of intervention
  T1 <- T - T0    # Number of time periods after intervention
  
  nlambda <- length(lambda_grid)  # Length of lambda grid
  na <- length(alpha_grid)        # Length of alpha grid
  
  err_alpha <- matrix(0, nrow = na, ncol = 1)         # Matrix to store error terms associated with each value of alpha
  lambda_opt_alpha <- matrix(0, nrow = na, ncol = 1)  # Matrix to store optimal lambda associated with each alpha
  
  # Rearrange matrices to have treated unit first to make the following procedure easier to generalize
  # In the general code, this will already have happened, but it's good to make sure this happens in case we run this function by itself
  
  Y <- as.matrix(cbind(Y[,ind_treatment], Y[,-ind_treatment]))
  Z <- as.matrix(cbind(Z[,ind_treatment], Z[,-ind_treatment]))
  X <- as.matrix(cbind(X[,ind_treatment], X[,-ind_treatment]))
  
  ################ Loop over alpha and lambda grids find optimal lambda for each alpha ################
  
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
      Z1_tr <- Z1                       # Pre-treatment outcomes for treated unit
      Z0_tr <- Z0                       # Pre-treatment outcomes for control units
      Z1_te <- as.matrix(Y1[-(1:T0),])  # Post-treatment outcomes for treated unit
      Z0_te <- as.matrix(Y0[-(1:T0),])  # Post-treatment outcomes for control units
      
      # Fit elastic net on pre-treatment outcomes
      V1 <- scale(Z1_tr, scale = FALSE) # Scale outcomes to have mean zero
      V0 <- scale(Z0_tr, scale = FALSE) # Scale outcomes to have mean zero
      
      # elastic net function: alpha = a and for all lambda points in lambda grid
      fit <- glmnet(x = V0, y = V1,
                    alpha = a,
                    lambda = lambda_grid,
                    standardize = FALSE,
                    intercept = FALSE) 
      w <- as.matrix(coef(fit, s = lambda_grid))  # Save coefficients of fit for weight
      w <- w[-1,]                                 # Delete intercept since it was forced to be zero anyways
      
      int <- t(as.matrix(apply(Z1_tr[,rep(1, nlambda)] - Z0_tr %*% w, 2, mean)))  # Calculate intercept as mean deviation between fit and true pre-treatment outcomes
      e <- Z1_te[,rep(1, nlambda)] - int[rep(1, T1),] - Z0_te %*% w               # Get post-treatment error for each lambda point
      err[i - 1,] <- colMeans(e ^ 2)                                              # SSR (over post-treatment time periods) for this control unit for each lambda point
    }
    
    # Find optimal lambda
    err <- apply(err, 2, mean)                  # Mean error over all control units
    ind_opt <- which.min(err)                   # Find lambda that minimizes error
    err_alpha[j] <- err[ind_opt]                # Save that error for this alpha
    lambda_opt_alpha[j] <- lambda_grid[ind_opt] # Save the corresponding optimal lambda value
  }
  
  
  ################ Find optimal alpha/lambda pair ################
  
  ind_opt <- which.min(err_alpha)         # Find alpha that minimizes error
  alpha_opt <- alpha_grid[ind_opt]        # Save the corresponding optimal alpha value
  lambda_opt <- lambda_opt_alpha[ind_opt] # Find associated lambda value
  
  # OUTPUT: optimal alpha and lambda
  out <- list("alpha"= alpha_opt, "lambda"= lambda_opt)
}

# Function to find the optimal tuning parameters (i.e. number of control units in set) for the best subset fit
tuning_parameters_best_subset<- function(Y,Z,X,ind_treatment=1){
  
  ## INPUT:
  #
  # Y:              matrix of outcome variables for treated unit and controls 
  # Z:              matrix of pre-treatment outcomes for treated unit and controls 
  # X:              matrix of covariates for treated unit and controls 
  # ind_treatment:  indicator which unit is the treatment unit in Y,Z and X matrices
  
  ## OUTPUT:
  #
  # n_opt:          optimal number of units in subset
  
  ################ Find parameters using inputted data ################
  
  N <- dim(Y)[2]              # Number of units
  T <- dim(Y)[1]              # Number of time periods
  T0 <- dim(Z)[1]             # Time of intervention
  T1 <- T - T0                # Number of time periods after intervention
  T0_tr <- floor(T0 * 2 / 3)  # Trimmed T0 value to prevent overfitting if N ~ T0
  
  # Normalize predictors
  #div <- as.matrix(apply(X, 1, sd)) # Matrix of standard deviations for each predictor
  #X <- X / div[,rep(1, N)] # Standardizes each predictor to have std 1
  # Add back if needed!
  
  # Rearrange matrices to have treated unit first to make the following procedure easier to generalize
  Y <- as.matrix(cbind(Y[,ind_treatment], Y[,-ind_treatment]))
  Z <- as.matrix(cbind(Z[,ind_treatment], Z[,-ind_treatment]))
  X <- as.matrix(cbind(X[,ind_treatment], X[,-ind_treatment]))
  
  ################ Find the optimal number of control units in a subset ################
  n_max <- N - 1 # maximum number of units in subset
  
  # Why N-2 and N-1???? 
  n_grid <- c(0:min(T0_tr - 1, n_max, N - 2)) # Define n_grid: only goes up to min(T0_tr-1, n_max) to prevent overfitting
  nn <- length(n_grid)                        # Number of points in n grid
  err <- matrix(0, nrow = N - 1, ncol = nn)   # Matrix for storage of errors for each unit and n
  c <- matrix(1, nrow = T0, ncol = 1)         # Columns of 1s to add as constant for fit in subset
  
  ################ Iterate over which unit is treated for crossvalidation ################
  
  for (i in 2:N) {
    
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
    
    # Fit the best subset model
    V1 <- Z1_tr
    
    ################ Iterate over subset size to find best subset of each size ################
    
    for (n in 0:(nn - 1)) {                                     # iterate over number of control units in subset
      subs_n <- combn(c(1:(N - 2)), n, simplify = FALSE)        # Generates all possible combinations of 1:(N-2) of size n 
      int <- matrix(0, nrow = 1, ncol = length(subs_n))         # Intercept for each subset combination
      w <- matrix(0, nrow = N - 2, ncol = length(subs_n))       # Weight matrix for each of the N-2 units and subset combination
      err_cur <- matrix(0, nrow = length(subs_n), ncol = 1)     # Matrix for storage of error for each subset
      
      for (j in 1:length(subs_n)) {                             # iterate over all possible combinations for a subset of length n
        sub <- subs_n[[j]]
        V0 <- cbind(c, Z0_tr[,sub])
        w_cur <- solve(t(V0) %*% V0, t(V0) %*% V1)              # Fit subset by OLS
        int[1,j] <- w_cur[1]
        w[sub,j] <- w_cur[-1]
        err_cur[j,1] <- mean((V1 - V0 %*% w_cur) ^ 2)           # Save error term of fit for each possible combination
      }
      
      # Optimal subset of size n
      j_opt <- which.min(err_cur)                               # Find subset combination that minimizes error
      e <- Z1_te - int[rep(1, T1),j_opt] - Z0_te %*% w[,j_opt]  # Find error for that subset combination
      err[i - 1,n + 1] <- mean(e ^ 2)                           # Save error for optimal subset
    }
  }
  
  # Optimal n
  err <- apply(t(scale(t(err))), 2, mean) # Mean error over all units
  ind_opt <- which.min(err) 
  n_opt <- n_grid[ind_opt]                # Find optimal number of units in subset
  
  return(n_opt)
}



#default 1 because when you run it in the normal function it will be in 1. Need option for standard errors
tuning_parameters_synth <- function(Y, Z, X, ind_treat=1){
 
  ## INPUT:
  #
  # Y:              matrix of outcome variables for treated unit and controls (note that these will  be Y_tune in general_estimate)
  # Z:              matrix of pre-treatment outcomes for treated unit and controls 
  # X:              matrix of covariates for treated unit and controls 
  # ind_treatment:  indicator which unit is the treatment unit in Y,Z and X matrices
  
  ## OUTPUT:
  #
  # vw:             the vweights for synthetic control estimation
  
  ############## Fit data to extract vweights ####################
  synth.out <- 
    synth(
      X1 = as.matrix(X[, ind_treat]),
      X0 = as.matrix(X[, -ind_treat]),
      Z1 = as.matrix(Z[, ind_treat]),
      Z0 = as.matrix(Z[, -ind_treat]),
      Margin.ipop=.005,Sigf.ipop=7,Bound.ipop=6
    )
  
  ############# output vweights ####################
  vw <- synth.out$solution.v
  return(vw)
}


# Function to find the weights for the control units using elastic net approach
find_weights_elastic_net <- function(Y, Z, X, alpha, lambda, lambda_grid, ind_treatment=1){
  
  ## INPUT:
  #
  # Y:              matrix of outcome variables for treated unit and controls 
  # Z:              matrix of pre-treatment outcomes for treated unit and controls 
  # X:              matrix of covariates for treated unit and controls 
  # alpha:          optimal value for elastic fit alpha
  # lambda:         optimal value for elastic fit lambda
  # lambda_grid:    lambda_grid used to find optimal lambda -> faster fit if whole sequence of lambdas is fitted instead of single one
  # ind_treatment:  indicator which unit is the treatment unit in Y,Z and X matrices
  
  ## OUTPUT:
  #
  # intercept:      intercept resulting from elastic net fit
  # weights:        weights resulting from elastic net fit
  
  ################ Find parameters using inputted data ################
  N <- dim(Y)[2] # Number of units
  
  # Normalize predictors
  #div <- as.matrix(apply(X, 1, sd)) # Matrix of standard deviations for each predictor
  #X <- X / div[,rep(1, N)] # Standardizes each predictor to have std 1
  ## ADD IF NEEDED
  
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
  
  ################ Fit elastic net using given tuning parameters ################
  fit <- glmnet(x = V0, y = V1,
                alpha = alpha,
                lambda = lambda_grid,
                standardize = FALSE, 
                intercept = FALSE)                # Function to fit elastic net with optimal alpha over lambda grid; no standardization nor intercept
  w <- as.matrix(coef(fit, s = lambda))           # Save coefficients for optimal lambda
  w <- w[-1,]                                     # Delete intercept since it was forced to be zero anyways
  int <- as.matrix(apply(Z1 - Z0 %*% w, 2, mean)) # Calculate intercept as mean deviation between fit and true pre-treatment outcomes
  
  # Output
  out <- list("intercept" = int, "weights" = w)
  
}



find_weights_synth <- function(Y, Z, X, ind_treat=1, vweight){
  
  ## INPUT:
  #
  # Y:              matrix of outcome variables for treated unit and controls 
  # Z:              matrix of pre-treatment outcomes for treated unit and controls 
  # X:              matrix of covariates for treated unit and controls 
  # ind_treatment:  indicator which unit is the treatment unit in Y,Z and X matrices
  # vweight:        list of vweights that we found from tuning_parameters_synth
  
  ## OUTPUT:
  #
  # w:              the estimated synthetic control weights and an intercept (which is zero)
  
  ############## Fit data to extract vweights ####################
  synth.out <- 
    synth(
      X1 = as.matrix(X[, ind_treat]),
      X0 = as.matrix(X[, -ind_treat]),
      Z1 = as.matrix(Z[, ind_treat]),
      Z0 = as.matrix(Z[, -ind_treat]),
      custom.v=as.numeric(vweight)
    )
  
  w<- list("intercept" = 0, "weights" = synth.out$solution.w)
}


find_weights_subset <- function(Y,Z,X,n_opt,ind_treatment=1){
  
  ## INPUT:
  #
  # Y:              matrix of outcome variables for treated unit and controls 
  # Z:              matrix of pre-treatment outcomes for treated unit and controls 
  # X:              matrix of covariates for treated unit and controls 
  # n_opt:          optimal number of control units
  # ind_treatment:  indicator which unit is the treatment unit in Y,Z and X matrices
  
  ## OUTPUT:
  #
  # intercept:      intercept resulting from best subset selection fit
  # weights:        weights resulting from best subset selection fit
  
  ############# Find parameters using the data ####################
  N <- dim(Y)[2]              # Number of units
  T <- dim(Y)[1]              # Number of time periods
  T0 <- dim(Z)[1]             # Time of intervention
  T1 <- T - T0                # Number of time periods after intervention
  T0_tr <- floor(T0 * 2 / 3)  # Trimmed T0 value to prevent overfitting if N ~ T0
  
  
  # Normalize predictors
  #div <- as.matrix(apply(X, 1, sd)) # Matrix of standard deviations for each predictor
  #X <- X / div[,rep(1, N)] # Standardizes each predictor to have std 1
  # Add back if needed!
  
  # Rearrange matrices to have treated unit first to make the following procedure easier to generalize
  Y <- as.matrix(cbind(Y[,ind_treatment], Y[,-ind_treatment]))
  Z <- as.matrix(cbind(Z[,ind_treatment], Z[,-ind_treatment]))
  X <- as.matrix(cbind(X[,ind_treatment], X[,-ind_treatment]))
  
  # Matrices for storage
  int <- matrix(0, nrow = 1, ncol = 1)
  w <- matrix(0, nrow = N - 1, ncol = 1)
  c <- matrix(1, nrow = T0, ncol = 1)
  
  # Fix treatment unit
  i <- 1
  Y1 <- as.matrix(Y[,i])
  Y0 <- as.matrix(Y[,-i])
  Z1 <- as.matrix(Z[,i])
  Z0 <- as.matrix(Z[,-i])
  X1 <- as.matrix(X[,i])
  X0 <- as.matrix(X[,-i])
  
  ############# Find optimal subset of size n ####################
  V1 <- Z1
  subs_n <- combn(c(1:(N - 1)), n_opt, simplify = FALSE) 
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
  
  # Pick out optimal subset
  j_opt <- which.min(err_cur)
  
  ############# Fit Model ####################
  
  # Pick fitted intercept and weights from optimal subset
  int <- as.matrix(int_s[1,j_opt])
  w <- w_s[,j_opt]
  
  # Output
  out <- list("intercept" = int, "weights" =w)
}

#find weights for difference-in-differences
find_weights_did <- function(Y, Z, X, ind_treatment=1){
  
  ## INPUT:
  #
  # Y:              matrix of outcome variables for treated unit and controls 
  # Z:              matrix of pre-treatment outcomes for treated unit and controls 
  # X:              matrix of covariates for treated unit and controls 
  # ind_treatment:  indicator which unit is the treatment unit in Y,Z and X matrices
  
  ## OUTPUT:
  #
  # intercept:      intercept resulting from diff-in-diff fit
  # weights:        weights resulting from diff-in-diff fit
  
  ########### Find Parameters from the data ###################
  N <- dim(Y)[2]
  
  #diff-in-diff weights units equally
  w <- matrix(1 / (N - 1), nrow = N - 1, ncol = 1) 
  int <- as.matrix(mean(Z[,ind_treatment]) - mean(Z[,-ind_treatment]))
  
  #output weights
  out <- list("intercept" = int, "weights" = w)
}

# Function to find the weights for the control units using constrained regression
find_weights_constr_reg <- function(Y,Z,X,ind_treatment=1){
  
  ## INPUT:
  #
  # Y:              matrix of outcome variables for treated unit and controls 
  # Z:              matrix of pre-treatment outcomes for treated unit and controls
  # X:              matrix of covariates for treated unit and controls 
  # ind_treatment:  indicator which unit is the treatment unit in Y,Z and X matrices
  
  ## OUTPUT:
  #
  # intercept:      intercept resulting from constrained regression fit
  # weights:        weights resulting from constrained regression fit
  
  ########### Find Parameters from the data ###################
  
  N <- dim(Y)[2]    # Number of units
  T <- dim(Y)[1]    # Number of time periods
  T0 <- dim(Z)[1]   # Time of intervention
  T1 <- T - T0      # Number of time periods after intervention
  print(N)
  # Matrices for storage
  int <- matrix(0, nrow = 1, ncol = 1)
  w <- matrix(0, nrow = N - 1, ncol = 1)
  
  # Normalize predictors
  #div <- as.matrix(apply(X, 1, sd)) # Matrix of standard deviations for each predictor
  #X <- X / div[,rep(1, N)] # Standardizes each predictor to have std 1
  # Add back if needed!
  
  # Fix treatment unit
  i <- ind_treatment
  Y1 <- as.matrix(Y[,i])
  Y0 <- as.matrix(Y[,-i])
  Z1 <- as.matrix(Z[,i])
  Z0 <- as.matrix(Z[,-i])
  X1 <- as.matrix(X[,i])
  X0 <- as.matrix(X[,-i])
  V1 <- Z1
  V0 <- Z0
  
  ########### Fit Constrained Regression ###################
  
  #should we explain this code so that people understand the steps?
  #H = t(V0)%*%V0
  #f = as.vector(-t(V0)%*%V1)
  
  #Aeq = rep(1,N-1)
  #beq = 1
  #lb = rep(0, N-1)
  #ub = rep(1,N-1)
  
  #w = quadprog(H, f, NULL, NULL, Aeq, beq, lb, ub)
  
  #out <- list("intercept" = 0,"weights"= w$x)
  V <- diag(length(V1))
  H <- t(V0) %*% V %*% (V0)
  a <- V1
  c <- -1 * c(t(a) %*% V %*% (V0))
  A <- t(rep(1, length(c)))
  b <- 1
  l <- rep(0, length(c))
  u <- rep(1, length(c))
  r <- 0
  res <- LowRankQP(Vmat = H, dvec = c, Amat = A, bvec = 1, 
                   uvec = rep(1, length(c)), method = "LU")
  w <- as.matrix(res$alpha)
  
  out <- list("intercept" = 0,"weights"= w)
  
}

#########################################################
#########################################################


general_se <- function(Y, Z, X, method=NULL, se_method="unit", tune_params=NULL,  ind_treatment=1){
  
  # INPUT
  #
  # ## INPUT:
  #
  # Y:              matrix of outcome variables for treated unit and controls 
  # Z:              matrix of pre-treatment outcomes for treated unit and controls 
  # X:              matrix of covariates for treated unit and controls 
  # method:         weights imputation method to be used. Could be "diff_in_diff", "elastic_net", "constr_reg", "synth" or "best_subset"
  # se_method:      chooses whether we're calculating as if the randomness were across units, across time, or both. Should be "unit", "time", or "unit_time". We may end up spitting out all values simulatneously
  # tune_params:    need tuning parameters for synth, best_subset and elastic_net. Currently, elastic net doesn't re-select tuning parameters for each loop, so simply
  #                     supply the optimal values for alpha and lambda found with the original run of tuning_parameters_elastic_net. Similarly, for best_subset, we currently
  #                     just re-use the optimal subset size for the orignal problem. Eventually, we might change tune_parameters for synth so that we supply the initial
  #                     v-weights found. However, currently, tune_params will be a list containing Y_tune, Z_tune, and X_tune named Y, X, Z.
  # ind_treatment:  column indicator for the treated unit
  #
  #
  # OUTPUT
  #
  # se:             the calculated standard error
  
  #must supply a supported method
  if(is.null(method) || (method != "diff_in_diff" && method != "elastic_net" && method != "constr_reg" && method != "synth" && method != "best_subset")){
    stop('Please specify one of the following methods: "diff_in_diff", "elastic_net", "constr_reg", "synth" or "best_subset"!')
  }
  
  #must supply a supported standard error method
  if(se_method!="unit" && se_method!="time" && se_method!="unit_time"){
    stop('Please specify one of the following methods for standard error calculation: "unit", "time" or "unit_time"')
    }else{
      
    ##################### Calculate standard errors ########################
    
    if(se_method=="unit"){
      se <- se_unit(Y, Z, X, method, tune_params, ind_treatment)
    }
    if(se_method=="time"){
      se <- se_time(Y, Z, X, method, tune_params, ind_treatment)
    }
    if(se_method=="unit_time"){
      se <- se_it(Y, Z, X, method, tune_params, ind_treatment)
    }
      
  }
  output <- se
}
  

# Standard error over units for elastic net, best subset, constrained regression and diff-in-diff
se_unit <- function(Y,Z,X, method, tune_params, ind_treatment=1){
  
  # INPUT
  #
  # Y:              matrix of outcome variables for treated unit and controls 
  # Z:              matrix of pre-treatment outcomes for treated unit and controls 
  # X:              matrix of covariates for treated unit and controls 
  # method:         weights imputation method to be used. Could be "diff_in_diff", "elastic_net", "constr_reg", "synth" or "best_subset"
  # tune_params:    need tuning parameters for synth, best_subset and elastic_net. Currently, elastic net doesn't re-select tuning parameters for each loop, so simply
  #                     supply the optimal values for alpha and lambda found with the original run of tuning_parameters_elastic_net. Similarly, for best_subset, we currently
  #                     just re-use the optimal subset size for the orignal problem. Eventually, we might change tune_parameters for synth so that we supply the initial
  #                     v-weights found. However, currently, tune_params will be a list containing Y_tune, Z_tune, and X_tune named Y, X, Z.
  # ind_treatment:  column indicator for the treated unit
  #
  #
  # OUTPUT
  #
  # se:             the calculated standard error
  
  ################# Find parameters from the data ####################
  
  N <- dim(Y)[2]                  # Number of units
  T <- dim(Y)[1]                  # Number of time periods
  T0 <- dim(Z)[1]                 # Time of intervention
  T1 <- T - T0                    # Number of time periods after intervention
  
  std_err_i <- matrix(0, N-1, T1) # Storage matrix for error terms
  
  # Define new Y,Z,X matrices without original treatment unit to feed into find_weights function
  Y <- Y[,-ind_treatment]
  Z <- Z[,-ind_treatment]
  X <- X[,-ind_treatment]
  
  ################# Loop across units to find standard error ####################
  
  for (i in 1:(N-1)) {
    
    # Find weights
    if(method == "elastic_net"){
      # tune_params need to be the alpha you found in the previous part, the lambda you found, then the lambda grid you used
      w <- find_weights_elastic_net(Y, Z, X, tune_params[[1]], tune_params[[2]], tune_params[[3]], i) 
    }
    if(method == "best_subset"){
      # tune_params needs to be the optimal n you found in the last part.
      w <- find_weights_subset(Y,Z,X,tune_params[[1]], i)
    }
    if(method == "constr_reg"){
      w <- find_weights_constr_reg(Y,Z,X, i)
    }
    if(method == "diff_in_diff"){
      w <- find_weights_did(Y,Z,X, i)
    }
    
    if(method == "synth"){
      #tune params needs to be a list with Y_tune, Z_tune, and X_tune
      #redraw vweights using the placebo dataset.
      v <- tuning_parameters_synth(tune_params$Y, tune_params$Z_tune, tune_params$X_tune, i)
      w <- find_weights_synth(Y, Z, X, i, v)
    }
    
    # Get standard error
    std_err_i[i,] <- (Y[-c(1:T0),i] - rep(w$intercept,T1) - Y[-c(1:T0),-i] %*% w$weights) ^ 2
  }
  std_err_i <- as.matrix(sqrt(apply(std_err_i, 2, mean)))
  return(std_err_i)
}

# Standard error over time 
se_time <- function(Y,Z,X, method, tune_params, ind_treatment=1){
  # INPUT
  #
  # Y:              matrix of outcome variables for treated unit and controls 
  # Z:              matrix of pre-treatment outcomes for treated unit and controls 
  # X:              matrix of covariates for treated unit and controls 
  # method:         weights imputation method to be used. Could be "diff_in_diff", "elastic_net", "constr_reg", "synth" or "best_subset"
  # tune_params:    need tuning parameters for synth, best_subset and elastic_net. Currently, elastic net doesn't re-select tuning parameters for each loop, so simply
  #                     supply the optimal values for alpha and lambda found with the original run of tuning_parameters_elastic_net. Similarly, for best_subset, we currently
  #                     just re-use the optimal subset size for the orignal problem. Eventually, we might change tune_parameters for synth so that we supply the initial
  #                     v-weights found. However, currently, tune_params will be a list containing Y_tune, Z_tune, and X_tune named Y, X, Z.
  # ind_treatment:  column indicator for the treated unit
  #
  #
  # OUTPUT
  #
  # se:             the calculated standard error
  
  T0 <- dim(Z)[1] # Time of intervention
  s <- floor(T0 / 2)
  std_err_t <- matrix(0, s, 1) # Storage matrix for error terms
  
  for (t in 1:s) {
    
    # Fix matrix to be according to time period
    Z <- Z[c(1:(T0 - t)),]
    
    # Find weights
    if(method == "elastic_net"){
      w <- find_weights_elastic_net(Y, Z, X, tune_params[[1]], tune_params[[2]], tune_params[[3]], ind_treatment) 
    }
    if(method == "best_subset"){
      w <- find_weights_subset(Y,Z,X, tune_params[[1]], ind_treatment)
    }
    if(method == "constr_reg"){
      w <- find_weights_constr_reg(Y,Z,X, ind_treatment)
    }
    if(method == "diff_in_diff"){
      w <- find_weights_did(Y,Z,X, ind_treatment)
    }
    if(method == "synth"){
      v <- tuning_parameters_synth(tune_params$Y, tune_params$Z_tune, tune_params$X_tune, ind_treatment)
      w <- find_weights_synth(Y, Z, X, ind_treatment, v)
    }
    
    # Get standard error
    std_err_t[t,1] <- (Y[T0 - t + 1,ind_treatment] - w$intercept - Y[T0 - t + 1,-ind_treatment] %*% w$weights) ^ 2
  }
  std_err_t <- as.matrix(sqrt(apply(std_err_t, 2, mean)))
  
  return(std_err_t)
}

# Standard error over units and time 
se_it <- function(Y,Z,X, method, tune_params, ind_treatment=1){
  
  # INPUT
  #
  # Y:              matrix of outcome variables for treated unit and controls 
  # Z:              matrix of pre-treatment outcomes for treated unit and controls 
  # X:              matrix of covariates for treated unit and controls 
  # method:         weights imputation method to be used. Could be "diff_in_diff", "elastic_net", "constr_reg", "synth" or "best_subset"
  # tune_params:    need tuning parameters for synth, best_subset and elastic_net. Currently, elastic net doesn't re-select tuning parameters for each loop, so simply
  #                     supply the optimal values for alpha and lambda found with the original run of tuning_parameters_elastic_net. Similarly, for best_subset, we currently
  #                     just re-use the optimal subset size for the orignal problem. Eventually, we might change tune_parameters for synth so that we supply the initial
  #                     v-weights found. However, currently, tune_params will be a list containing Y_tune, Z_tune, and X_tune named Y, X, Z.
  # ind_treatment:  column indicator for the treated unit
  #
  #
  # OUTPUT
  #
  # se:             the calculated standard error
  
  N <- dim(Y)[2] # Number of units
  T0 <- dim(Z)[1] # Time of intervention
  s <- floor(T0 / 2)
  std_err_it <- matrix(0, N - 1, 1) # Storage matrix for error
  
  # Define new Y,Z,X matrices without original treatment unit to feed into find_weights function
  Y <- Y[,-ind_treatment]
  Z <- Z[,-ind_treatment]
  X <- X[,-ind_treatment]
  
  for (i in 1:(N-1)) {
    
    std_err_temp <- matrix(0, s, 1) # Storage matrix for temporary standard error for each time period
    
    for (t in 1:s) {
      
      # Fix matrix to be according to time period
      Z_temp <- Z[c(1:(T0 - t)),]
      
      # Find weights
      if(method == "elastic_net"){
        w <- find_weights_elastic_net(Y, Z_temp, X, tune_params[[1]], tune_params[[2]], tune_params[[3]], i) 
      }
      if(method == "best_subset"){
        w <- find_weights_subset(Y,Z_temp,X,tune_params[[1]], i)
      }
      if(method == "constr_reg"){
        w <- find_weights_constr_reg(Y,Z_temp,X, i)
      }
      if(method == "diff_in_diff"){
        w <- find_weights_did(Y,Z_temp,X, i)
      }
      if(method == "synth"){
        v <- tuning_parameters_synth(tune_params$Y, tune_params$Z_tune, tune_params$X_tune, i)
        w <- find_weights_synth(Y, Z, X, i, v)
      }
      
      # Get standard error
      std_err_temp[t,1] <- (Y[T0 - t + 1,i] - w$intercept - Y[T0 - t + 1,-i] %*% w$weights) ^ 2
    }
    std_err_temp <- as.matrix(apply(std_err_temp, 2, mean))
    std_err_it[i,1] <- std_err_temp
  }
  std_err_it <- as.matrix(sqrt(apply(std_err_it, 2, mean)))
  
  return(std_err_it)
}

general_treatment_plot <- function(fit_obj, legend = NULL, time_vec = NULL){
  
  # INPUT
  #
  # fit_obj: Object returned by the function general_estimate
  # legend: vector supplied if legend is wished for (e.g. c("Actual data", "Fitted data"))
  # time_vec: vector of start and end point for x-axis of the plot (e.g. c("1960","1990"))
  #
  #
  # OUTPUT
  #
  # Plot of the true data and estimated data
  
  if(is.null(time_vec)){
    t = length(fit_obj$Y_true)
    plot(1:t, fit_obj$Y_true, type = "l", lty = 2, xlim = c(1,t), col = "black",  xlab = " ", ylab = "", las = 1, bty = 'L')
    lines(1:t, fit_obj$Y_est, lty = 1, col= "red")
    abline(v = (fit_obj$T_0-1), col = "black")
    if(!is.null(legend)){
      legend("topleft",legend= legend, col = c("black", "red"),lty=c(2,1), ncol=1, bty = 'n', cex = 0.7)
    }
  }else{
    plot(time_vec[1]:time_vec[2], fit_obj$Y_true, type = "l", lty = 2, xlim = c(time_vec[1],time_vec[2]), col = "black",  xlab = " ", ylab = "", las = 1, bty = 'L')
    lines(time_vec[1]:time_vec[2], fit_obj$Y_est, lty = 1, col= "red")
    abline(v = (time_vec[1]+fit_obj$T_0-1), col = "black")
    if(!is.null(legend)){
      legend("topleft",legend= legend, col = c("black", "red"),lty=c(2,1), ncol=1, bty = 'n', cex = 0.7)
    }
  }
}

general_std_error_plot <- function(fit_obj, legend = NULL, time_vec = NULL){
  
  # INPUT
  #
  # fit_obj: Object returned by the function general_estimate
  # legend: vector supplied if legend is wished for (e.g. c("Actual data", "Fitted data"))
  # time_vec: vector of start and end point for x-axis of the plot (e.g. c("1960","1990"))
  #
  #
  # OUTPUT
  #
  # Plot of the estimated treatment effect and confidence intervals
  
  start = fit_obj$T_0
  end = length(fit_obj$Y_true)
  tau <- fit_obj$Y_true[start:end]-fit_obj$Y_est[start:end] 
  std_err <- fit_obj$std_err_i
  
  if(is.null(time_vec)){
    t = end
    plot(1:t, tau, type = "l", lty = 1, xlim = c(1,t), col = "black",  xlab = " ", ylab = "", las = 1, bty = 'L')
    lines(1:t, tau+1.96*std_err, lty = 3, col= "red")
    lines(1:t, tau-1.96*std_err, lty = 3, col= "red")
    abline(h = 0, col= "black")
    if(!is.null(legend)){
      legend("topleft",legend= legend, col = c("black", "red"),lty=c(1,3), ncol=1, bty = 'n', cex = 0.7)
    }
  }else{
    start = time_vec[1]+fit_obj$T_0
    end = time_vec[2]
    plot(start:end, tau, type = "l", lty = 1, xlim = c(start,end), col = "black",  xlab = " ", ylab = "", las = 1, bty = 'L')
    lines(start:end, tau+1.96*std_err, lty = 3, col= "red")
    lines(start:end, tau-1.96*std_err, lty = 3, col= "red")
    abline(h = 0, col= "black")
    if(!is.null(legend)){
      legend("topleft",legend= legend, col = c("black", "red"),lty=c(1,3), ncol=1, bty = 'n', cex = 0.7)
    }
  }
}

general_weights_plot <- function(fit_obj, control_names, method){
  
  # INPUT
  #
  # fit_obj: Object returned by the function general_estimate
  # control_names: vector containing names of the control units
  # method: string indicating which method was used
  #
  #
  # OUTPUT
  #
  # Plot of the regression weights
  
  weights <- as.data.frame(cbind(control_names,fit_obj$w))
  colnames(weights) <- c("controls","w")
  weights$w <- as.numeric(as.character(weights$w))
  
  p <- ggplot(weights, aes(x=controls, y=w))+geom_bar(stat="identity", fill = "blue",color ="black", show.legend = FALSE)+labs(title="",x="", y = method)+scale_fill_manual(values = c("royalblue"))+ylim(-1, 1)+coord_flip()
  p
  
  
}