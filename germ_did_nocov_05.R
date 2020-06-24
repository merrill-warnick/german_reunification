## Clear
rm(list = ls())

## Options

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
#a <- 1/3 # What does a stand for??? - never used
K <- dim(X)[1] # Number of predictors
N <- dim(Y)[2] # Number of units (controls+ treatment unit)
T <- dim(Y)[1] # Number of time periods
T0 <- dim(Z)[1] # Time of intervention
T1 <- T - T0 # Number of time periods after intervention
#T0_tr <- floor(T0 * 2 / 3) # What is this???
#T0_te <- T0 - T0_tr # What is this??? - never used

# Normalize predictors
div <- as.matrix(apply(X, 1, sd)) # Matrix of standard deviations for each predictor
X <- X / div[,rep(1, N)] # Standardizes each predictor to have std 1


###########################
####### Fit Model #########
###########################

## Fit the model

# Initialize storage matrices
int_did <- matrix(0, nrow = 1, ncol = 1) # Intercept for DiD
w_did <- matrix(1 / (N - 1), nrow = N - 1, ncol = 1) # Weights are all equal 
Y_did <- matrix(0, nrow = T, ncol = 1)
Y_did <- matrix(0, nrow = T, ncol = 1) # Why twice???

###############
# Define each Z0,Z1,etc correct -> think about what to do here..not necessary if do dataprep in beginning

# West Germany
i <- 1
Y1 <- as.matrix(Y[,i])
Y0 <- as.matrix(Y[,-i])
Z1 <- as.matrix(Z[,i])
Z0 <- as.matrix(Z[,-i])
X1 <- as.matrix(X[,i])
X0 <- as.matrix(X[,-i])

###########
# Fit DiD
w <- w_did
int_did <- as.matrix(mean(Z1) - mean(Z0))
Y_did <- int_did[rep(1, T),] + Y0 %*% w # Estimated Y (no treatment)
Y_true <- Y1

#################################
####### Standard Errors #########
#################################

# Over units
std_err_i <- matrix(0, T1, N - 1)
for (i in 2:N) {
  
  # Define matrices appropriately; take treatment out of control!
  Y1 <- as.matrix(Y[,i])
  Y0 <- as.matrix(Y[,-c(1,i)])
  Z1 <- as.matrix(Z[,i])
  Z0 <- as.matrix(Z[,-c(1,i)])
  X1 <- as.matrix(X[,i])
  X0 <- as.matrix(X[,-c(1,i)])
  
  # Fit DiD
  w <- matrix(1 / (N - 2), N - 2, 1) # N-2 instead of N-1 because treatment unit also out of control group here
  int <- as.matrix(mean(Z1) - mean(Z0))
  std_err_i[,i - 1] <- (Y1[-c(1:T0),] - int[rep(1, T1),] - Y0[-c(1:T0),] %*% w) ^ 2
}
std_err_i <- as.matrix(sqrt(apply(std_err_i, 1, mean))) # Std Er is sqrt(SSR) over units

# Over time
s <- floor(T0 / 2) # Number of periods for standard errors
std_err_t <- matrix(0, s, 1) # Storage matrix

# Fix matrices for Y and X for time-varying standard errors
Y1 <- as.matrix(Y[,1])
Y0 <- as.matrix(Y[,-1])
X1 <- as.matrix(X[,1])
X0 <- as.matrix(X[,-1])

for (t in 1:s) {
  # Have varying Z matrices: Z = pre-treatment outcome over which MSPE to be minimized
  Z1 <- as.matrix(Z[c(1:(T0 - t)),1]) # Pretend intervention happens earlier
  Z0 <- as.matrix(Z[c(1:(T0 - t)),-1])
  
  # Fit DiD
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
  # Do same time-varying as before but now for each unit
  for (t in 1:s) {
    Z1 <- as.matrix(Z[c(1:(T0 - t)),i])
    Z0 <- as.matrix(Z[c(1:(T0 - t)),-c(1,i)])
    
    # Fit DiD
    w <- matrix(1 / (N - 2), N - 2, 1)
    int <- as.matrix(mean(Z1) - mean(Z0))
    std_err_temp[t,1] <- (Y1[T0 - t + 1,] - int - Y0[T0 - t + 1,] %*% w) ^ 2
  }
  std_err_temp <- as.matrix(apply(std_err_temp, 2, mean))
  std_err_it[i - 1,1] <- std_err_temp # save time-varying std err for each unit
}
std_err_it <- as.matrix(sqrt(apply(std_err_it, 2, mean))) # then do sqrt of mean for each unit that is not treatment


##############################
####### Save Results #########
##############################

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
