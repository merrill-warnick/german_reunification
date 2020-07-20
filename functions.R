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

general_estimate <- function(data_df, method = NULL, prep_params, special_params = NULL){
  
  ## INPUT:
  #
  # data: specify whether just frame 
  # method: "diff_in_diff", "elastic_net", "constr_reg", "synth","best_subset"; 
  # special_params: list of lambda grid and alpha grid over which we are optimizing; special_params = list(lambda_grid,alpha_grid)

  #prep parameters is the vector for data prep
  
  # I think I'll change it so that we have a special_parameters input?
  #then it will have lambda_grid/alpha_grid and all of the parameters for synth
  #but it's null by default?
  # Yeah that's fine
  
  # Are we going to make it so that the user doesn't have to plug in inputs for the methods they don't use?
  # Yes! So special_params should be NULL by default indeed if we don't pick elastic net as a default method
  # But that's the only thing that should be by default, right?
  
  ## OUTPUT:
  #
  # int: estimated intercept
  # w: estimated weights
  # Y_est: fitted Y values
  # Y_true: true Y-values
  # alpha_opt: optimal alpha (if method is elastic net)
  # lambda_opt: optimal lambda (if method is elastic net)
  # n_opt: optimal number of control units (if method is best subset)
  # std_err_i: standard error over units
  # std_err_t: standard error over time
  # std_err_it: standard error over unit and time
  # T_0: time point of intervention
  # T_1: time points after intervention
  
  # Make sure that one of the methods that is feasible is specified
  if(is.null(method)){# | method != "diff_in_diff" | method != "elastic_net" | method != "constr_reg" | method != "synth" | method != "best_subset"){
    stop('Please specify one of the following methods: "diff_in_diff", "elastic_net", "constr_reg", "synth" or "best_subset"!')
  }else{
    ####
    # Dataprep

    # prep_params[[1]] = pred: vector containing string with predictor variables 
    # prep_params[[2]] = dep: string specifying which one is dependent variable
    # prep_params[[3]] = u: integer identifying unit variable (which column in data frame specifies index!)
    # prep_names[[4] = t: integer identifying time variable (which column in data frame specifies time!)
    # prep_names[[5]] = spec: list of special predictors
    # prep_names[[6]] = i: treatment identifier -> should be just ind_treatment, no? then can delete as separate input in function I guess
    # prep_names[[7]] = j: 
    # prep_names[[8]] = subs: together with j identifies control identifiers -> should we just give it as a single input?
    # prep_names[[9]] = years: vector specifying [1]: start predictor priors, [2]: end predictor priors
    #                                            [3]: start time optimize ssr, [4]: end time optimize ssr
    #                                            [5]: start time plot, [6]: end time plot
    # prep_names[[10]] = names: unit names variable (which column in data frame specifies unit names!)
      
    data <- prep_data(data_df, prep_params[[1]], prep_params[[2]], prep_params[[3]], prep_params[[4]], prep_params[[5]], prep_params[[6]], prep_params[[7]], prep_params[[8]], prep_params[[9]], prep_params[[10]])
    Y<- data$Y
    Z<- data$Z
    X<- data$X
    
    ####
    # Get weights
    
    # Elastic Net: Find tuning parameters and weights
    if(method == "elastic_net"){
      params <- tuning_parameters_elastic_net(Y,Z,X, special_params[[1]], special_params[[2]])
      w <- find_weights_elastic_net(Y, Z, X, params$alpha, params$lambda, special_params[[1]]) 
    }
    
    #synthetic control: find tuning parameters and weights
    if(method == "synth"){
      #probably need to pick out parameters from a parameter list or something like that.
      #hmmm it's also important to know that for this part we need to plug in the raw dataframe and not the dataprepped dataframe....
      #yeah synth just has a lot of nasty extra parameters because we have to reuse the dataprep function several times.
      #Here are some ideas:
      #   rewrite it so that we pick out the parameter values inside the function. This might be less clear to readers what to plug in tho
      #   write up a principled way to automatically figure out what the tuning years and the other years should be. Plausibly we could do this
      #       if we knew which column was the time column and what year the treatment was in. Yeah actually, now that I think about it, we could
      #       probably rewrite this to take in general Y Z X matrices without a time column and get fewer frills. I can play around with dataprep on 
      #       my branch I guess to figure this out.
      
      #I should rewrite these a bit/think harder about it.
      #or at least write it so that it can take a general parameter set or something like that?
      
      #realistically, it seems very annoying to use the special parameter vector...maybe I should do the picking out of the vector inside the function?
      #maybe actually I should pass in like "shared parameters", "v parameters" and "w parameters".
      # Another thing to think about is that we should probably think about a way to automatically figure out what the years should be from T, T1, T0 or whatever
      #maybe we split it up into index_parameters/time_parameters/other_parameters? idk. and i guess that some of these are actually names for columns, right? hmmmm
      
      v <- tuning_parameters_synth(data_df, pred, y, u, t, spec, ind_treatment, ind_treatment, cont_set, predyear0, predyear1, optyear0, optyear1,  year0, year1, names)
      w <- find_weights_synth(data_df, pred, y, u, t, spec, ind_treament, ind_treatment, cont_set, predyear0, predyear1, optyear0, optyear1, year0, year1, names, v, FALSE)
    }
    
    # Best subset: Find tuning parameter and weights
    if(method == "best_subset"){
      n_opt <- tuning_parameters_best_subset(Y,Z,X)
      w <- find_weights_subset(Y,Z,X,n_opt)
    }
    
    # Diff-in-diff: Find weights
    if (method == "diff_in_diff"){
      w <- find_weights_did(Y, Z, X)
    }
    
    # Constrained regression: Find weights
    if(method == "constr_reg"){
      w <- find_weights_constr_reg(Y,Z,X)
    }
    
    ###
    # Get estimate
    Y_est = rep(w$intercept,nrow(Y)) + Y[,-1]%*%w$weights
    Y_true = Y[,1]
    
    ###
    # Get standard error
    std_err_i = 0
    std_err_t = 0
    std_err_it = 0
    
    ###
    # Output
    if(method == "elastic_net"){
      out <- list("int" = w$intercept, "w" = w$weights, "Y_est" = Y_est, "Y_true" = Y_true, "alpha_opt" = params$alpha, "lambda_opt" = params$lambda,"std_err_i" = std_err_i, "std_err_t" = std_err_t, "std_err_it" = std_err_it, "T_0"= dim(Z)[1],"T_1"=dim(Y)[1]-dim(Z)[1])
    }else{
      if(method == "best_subset"){
        out <- list("int" = w$intercept, "w" = w$weights, "Y_est" = Y_est, "Y_true" = Y_true,"n_opt" = n_opt, "std_err_i" = std_err_i, "std_err_t" = std_err_t, "std_err_it" = std_err_it, "T_0"= dim(Z)[1],"T_1"=dim(Y)[1]-dim(Z)[1])
      }else{
        out <- list("int" = w$intercept, "w" = w$weights, "Y_est" = Y_est, "Y_true" = Y_true,"std_err_i" = std_err_i, "std_err_t" = std_err_t, "std_err_it" = std_err_it, "T_0"= dim(Z)[1],"T_1"= dim(Y)[1]-dim(Z)[1])
      }
    }
  }
}




###############################
##### Auxiliary functions #####
###############################

#I guess with the way that we're likely writing the prep_data function, it would probably make sense to 
#just specify the treatment column as one for all subsequent functions. If prep_data is inside the function, this makes sense
#but if it's outside the function we don't have to do it like that I think.
#I guess for now I'll just put synth inside the function and leave this out until Lea and I meet.
#I think that we should have like a "special parameters" input that lets you put in the method-specific parameters.
#I thought about this more up above by the way.
prep_data <- function(d, pred, dep, u, t, spec,i, j, subs, years, names){
  
  dataprep.out <-
    dataprep(
      foo = d,
      predictors    = pred,
      dependent     = dep,
      unit.variable = u,
      time.variable = t,
      special.predictors = spec,
      treatment.identifier = i,
      controls.identifier = subs[-j],
      time.predictors.prior = years[1]:years[2],
      time.optimize.ssr = years[3]:years[4],
      unit.names.variable = names,
      time.plot = years[5]:years[6]
    )
  
  #######################################
    #prepping the pieces
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
    
    datlist <- list( x <- X, y <- Y, z <- Z)
    names(datlist) <- c( "X", "Y", "Z")
    output<-datlist
}

tuning_parameters_elastic_net <- function(Y,Z,X,lambda_grid, alpha_grid, ind_treatment=1){
  
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

tuning_parameters_best_subset<- function(Y,Z,X,ind_treatment=1){
  
  ## INPUT:
  #
  # Y: matrix of outcome variables for treated unit and controls 
  # Z: matrix of pre-treatment outcomes for treated unit and controls 
  # X: matrix of covariates for treated unit and controls 
  # ind_treatment: indicator which unit is the treatment unit
  
  ## OUTPUT:
  #
  # n_opt: optimal number of units in subset
  
  # Parameters
  N <- dim(Y)[2] # Number of units
  T <- dim(Y)[1] # Number of time periods
  T0 <- dim(Z)[1] # Time of intervention
  T1 <- T - T0 # Number of time periods after intervention
  T0_tr <- floor(T0 * 2 / 3) # what is this?? - here it is actually used
  
  
  # Normalize predictors
  #div <- as.matrix(apply(X, 1, sd)) # Matrix of standard deviations for each predictor
  #X <- X / div[,rep(1, N)] # Standardizes each predictor to have std 1
  # Add back if needed!
  
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
  
  return(n_opt)
}

tuning_parameters_synth <- function(d, pred, y, u, t, spec, i,j,cont_set, predyear0, predyear1, optyear0, optyear1, year0, year1, names){
  #d is the dataframe of the panel data
  #pred is a string of predictor variables
  #y is the string name of the dependent variable
  #u is the value of the unit identifier column
  #t is the value of the time identifier column 
  #spec is a list of special predictors that should be what you're plugging into special.predictors
  #  for the dataprep function
  #i is the index of the treatment identifier column
  #j helps with the index of the control identifier columns in case you need it
  #cont_set is whatever the set of controls that you're picking out it
  #predyearX are the first and last years that you want to use as predictors
  #optyearX are the first and last years that you want to poptimize over for crossvalidation.
  #names is the column of name identifiers for the units
  #yearX is the first and last year that you want for time.plot
  
  dataprep.out <-
    dataprep(
      
      foo = d,
      
      predictors    = pred,
      
      dependent     = y,
      
      unit.variable = u,
      
      time.variable = t,
      
      special.predictors = spec,
      
      treatment.identifier = i,
      
      controls.identifier = cont_set[-j],
      
      time.predictors.prior = predyear0:predyear1,
      
      time.optimize.ssr = optyear0:optyear1,
      
      unit.names.variable = names,
      
      time.plot = year0:year1
    )
  
  #fit training model to pull out vweights
  synth.out <- 
    synth(
      data.prep.obj=dataprep.out,
      
      Margin.ipop=.005,Sigf.ipop=7,Bound.ipop=6
    )
  output <- synth.out$solution.v
}

find_weights_elastic_net <- function(Y, Z, X, alpha, lambda, lambda_grid, ind_treatment=1){
  
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

#find_weights_synth <- function(d, pred, y, u, t, spec, i, j,cont_set, predyear0, predyear1, optyear0, optyear1, names, year0, year1, names, vweight, yinclude){
#  
#  dataprep.out <-
#    dataprep(
#      foo = d,
#      predictors    = pred,
#      dependent     = y,
#      unit.variable = u,
#      time.variable = t,
#      special.predictors = spec,
#      treatment.identifier = i,
#     controls.identifier = cont_set[-j],
 #     
#      time.predictors.prior = predyear0:predyear1,
#      time.optimize.ssr = optyear0:optyear1,
#      unit.names.variable = names,
#     time.plot = year0:year1
#    )
#  
#  
#  synth.out <- synth(
#    data.prep.obj=dataprep.out,
#    custom.v=as.numeric(vweight)
#  )
#  
#  w <- synth.out$solution.w
#  out <- list("intercept" = 0, "weights" = w)
#  if(yinclude==TRUE){
#    out<- list("intercept" = 0, "weights" = w, "Y1"<-dataprep.out$Y1, "Y0"<-dataprep.out$Y0)
#  }
#}

find_weights_subset <- function(Y,Z,X,n_opt,ind_treatment=1){
  
  ## INPUT:
  #
  # Y: matrix of outcome variables for treated unit and controls 
  # Z: matrix of pre-treatment outcomes for treated unit and controls 
  # X: matrix of covariates for treated unit and controls 
  # n_opt: optimal number of control units
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
  #div <- as.matrix(apply(X, 1, sd)) # Matrix of standard deviations for each predictor
  #X <- X / div[,rep(1, N)] # Standardizes each predictor to have std 1
  # Add back if needed!
  
  # Rearrange matrices to have treated unit first
  Y <- as.matrix(cbind(Y[,ind_treatment], Y[,-ind_treatment]))
  Z <- as.matrix(cbind(Z[,ind_treatment], Z[,-ind_treatment]))
  X <- as.matrix(cbind(X[,ind_treatment], X[,-ind_treatment]))
  
  
  ################################
  ####### OPTIMAL SUBSET #########
  ################################
  
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

find_weights_did <- function(Y, Z, X, ind_treatment=1){
  N <- dim(Y)[2]
  w <- matrix(1 / (N - 1), nrow = N - 1, ncol = 1) 
  int <- as.matrix(mean(Z[,ind_treatment]) - mean(Z[,-ind_treatment]))
  out <- list("intercept" = int, "weights" = w)
}

find_weights_constr_reg <- function(Y,Z,X,ind_treatment=1){
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

#########################################################
#########################################################
# Std Error to be put into general function for elastic net, best_subset and constrained regression

# Not sure if this works out...gotta try once done! Need to think about it still quickly if correct application with weights function
# (e.g. matrices and indicator for treatment variable)

#I added diff in diff to the standard error function

#to put synth in the general function it could look something like this 

#standard_errors <- function(data, params){
#run prep_data to pick out parameters
#if(method==synth){
#se <- specific synth funciton
#}
#else{
#
#
#put in the code from below
#}
#}



# Over units
se_unit <- function(Y,Z,X, method, alpha= NULL, lambda= NULL, n_opt = NULL, ind_treatment=1){
  # Over Units 
  N <- dim(Y)[2]
  T <- dim(Y)[1]
  T0 <- dim(Z)[1]
  T1 <- T-T0
  std_err_i <- matrix(0, N-1, T1)
  
  # Define new Y,Z,X matrices without original treatment unit to feed into find_weights function
  Y <- Y[,-ind_treatment]
  Z <- Z[,-ind_treatment]
  X <- X[,-ind_treatment]
  
  for (i in 1:(N-1)) {
    
    # Find weights
    if(method == "elastic_net"){
      w <- find_weights_elastic_net(Y, Z, X, alpha, lambda, lambda_grid, i) 
    }
    if(method == "best_subset"){
      w <- find_weights_subset(Y,Z,X,n_opt, i)
    }
    if(method == "constr_reg"){
      w <- find_weights_constr_reg(Y,Z,X, i)
    }
    if(method == "diff_in_diff"){
      w <- find_weights_did(Y,Z,X, i)
    }
    
    # Get standard error
    std_err_i[i,] <- (Y[-c(1:T0),i] - rep(w$intercept,T1) - Y[-c(1:T0),-i] %*% w$weights) ^ 2
  }
  std_err_i <- as.matrix(sqrt(apply(std_err_i, 2, mean)))
  return(std_err_i)
}

# Over Time
se_time <- function(Y,Z,X, method, alpha= NULL, lambda= NULL,n_opt = NULL, ind_treatment=1){
  T0 <- dim(Z)[1]
  s <- floor(T0 / 2)
  std_err_t <- matrix(0, s, 1)
  
  for (t in 1:s) {
    
    # Fix matrix to be according to time period
    Z <- Z[c(1:(T0 - t)),]
    
    # Find weights
    if(method == "elastic_net"){
      w <- find_weights_elastic_net(Y, Z, X, alpha, lambda, lambda_grid, ind_treatment) 
    }
    if(method == "best_subset"){
      w <- find_weights_subset(Y,Z,X,n_opt, ind_treatment)
    }
    if(method == "constr_reg"){
      w <- find_weights_constr_reg(Y,Z,X, ind_treatment)
    }
    if(method == "diff_in_diff"){
      w <- find_weights_did(Y,Z,X, i)
    }
    
    std_err_t[t,1] <- (Y[T0 - t + 1,ind_treatment] - w$intercept - Y[T0 - t + 1,-ind_treatment] %*% w$weights) ^ 2
  }
  std_err_t <- as.matrix(sqrt(apply(std_err_t, 2, mean)))
  
  return(std_err_t)
}

# Over Units and Time
se_it <- function(Y,Z,X, method, alpha= NULL, lambda= NULL,n_opt = NULL, ind_treatment=1){
  # Over Units and Time
  N <- dim(Y)[2]
  T0 <- dim(Z)[1]
  s <- floor(T0 / 2)
  std_err_it <- matrix(0, N - 1, 1)
  # Define new Y,Z,X matrices without original treatment unit to feed into find_weights function
  Y <- Y[,-ind_treatment]
  Z <- Z[,-ind_treatment]
  X <- X[,-ind_treatment]
  
  for (i in 1:(N-1)) {
    
    std_err_temp <- matrix(0, s, 1)
    for (t in 1:s) {
      # Fix matrix to be according to time period
      Z_temp <- Z[c(1:(T0 - t)),]
      
      # Find weights
      if(method == "elastic_net"){
        w <- find_weights_elastic_net(Y, Z_temp, X, alpha, lambda, lambda_grid, i) 
      }
      if(method == "best_subset"){
        w <- find_weights_subset(Y,Z_temp,X,n_opt, i)
      }
      if(method == "constr_reg"){
        w <- find_weights_constr_reg(Y,Z_temp,X, i)
      }
      if(method == "diff_in_diff"){
        w <- find_weights_did(Y,Z_temp,X, i)
      }
      
      std_err_temp[t,1] <- (Y[T0 - t + 1,i] - w$intercept - Y[T0 - t + 1,-i] %*% w$weights) ^ 2
    }
    std_err_temp <- as.matrix(apply(std_err_temp, 2, mean))
    std_err_it[i,1] <- std_err_temp
  }
  std_err_it <- as.matrix(sqrt(apply(std_err_it, 2, mean)))
  
  return(std_err_it)
}


#############Synth standard error functions#######################

N <- dim(Y)[2]
T <- dim(Y)[1]
T0 <- dim(Z)[1]
T1 <- T-T0

#need to think about inputting the right parameters, or where we should be doing that rather.
#probably just have a parameter vector? We can worry aobut this a bit later.

se_unit_synth <- function(data, N, T, T0, pred, y, u, t, cspec, spec, cont_set, cyears, years, names){
  
  T1<- T - T0
  
  std_err_i <- matrix(0, N - 1, T1)
  for (j in 1:(N - 1)) {
    i <- cont_set[j]
    vw <- tuning_parameters_synth(data, pred, y, u, t, cspec, i,j,cont_set, cyears[1], cyears[2], cyears[3], cyears[4], names, cyears[5], cyears[6])
   
    w <- find_weights_synth(data, pred, y, u, t, spec, i, j, cont_set, years[1], years[2], years[3], years[4], years[5], years[6], names, vw, TRUE)
    
    std_err_i[j,] <- (w$Y1[-c(1:T0),] - w$intercept - w$Y0[-c(1:T0),] %*% w$weights) ^ 2
  }
  
  std_err_i <- as.matrix(sqrt(apply(std_err_i, 2, mean)))
}



se_time <- function(data, N, T, T0,  pred, y, u, t, cspec,i,j, spec, cont_set, cyears, years, names){
  
  s <- floor(T0 / 2)
  std_err_t <- matrix(0, s, 1)
  vweights <- matrix(0,s,6)
  for (k in 1:s) {
    cat('(Std. Error) Over Time t =', toString(k), '\n')
    
    vw <- find_vweights(data, pred, y, u, t, cspec, i, j, cont_set, cyears[1], cyears[2], cyears[3], cyears[4], names, cyears[5], cyears[6],names)
    
    w <- find_weights_synth(data, pred, y, u, t, spec, i,j, cont_set, years[1], years[2], years[3], years[4]-k, years[5], years[6],names, vw, TRUE)
    
    std_err_t[k,1] <- (w$Y1[T0 - k + 1,] - w$intercept - w$Y0[T0 - k + 1,] %*% w$weights) ^ 2
  }
  std_err_t <- as.matrix(sqrt(apply(std_err_t, 2, mean)))
}



se_unit_time <- function(d, N, T, T0, pred, y, u, t, cspec, spec, cont_set, cyears, years, names){
  
  s <- floor(T0 / 2)
  
  std_err_it <- matrix(0, N - 1, 1)
  for (j in 1:(N - 1)) {
    i <- cont_set[j]
    std_err_temp <- matrix(0, s, 1)
    for (k in 1:s) {
      
      vw <- find_vweights(d, pred, y, u, t, cspec, i,j, cont_set, cyears[1], cyears[2], cyears[3], cyears[4], cyears[5], cyears[6], names)
      
      w <- find_weights_synth(d, pred, y, u, t, spec, i,j, cont_set, years[1], years[2], years[3], years[4] - k, years[5], years[6], names, vw, TRUE)
      
      std_err_temp[k,1] <- (w$Y1[T0 - k + 1,] - w$intercept - w$Y0[T0 - k + 1,] %*% w$weights) ^ 2
    }
    std_err_temp <- as.matrix(apply(std_err_temp, 2, mean))
    std_err_it[j,1] <- std_err_temp
  }
  
  std_err_it <-  as.matrix(sqrt(apply(std_err_it, 2, mean)))
}
