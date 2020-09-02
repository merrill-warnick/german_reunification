# Other data sets

## Clear
rm(list = ls())

library(R.matlab)
library(shape)
library(ggplot2)

source('functions.R')

# Load data

####################
#### California ####
####################

data <- read.table("MLAB_data.txt")
## Built Indices
# California is state no 3, stored in the last column no 39
index_tr <- c(39)
# 38 Control states are 1,2 & 4,5,...,38, stored in columns 1 to 38
index_co <- c(1:38)

# Predcitors are stored in rows 2 to 8
index_predict <- c(2:8)
# Outcome Data is stored in rows 9 to 39; for 1970, 1971,...,2000
index_Y <- c(9:39)

# The pre-treatment time periods are 1970--1980
index_pre <- c(1:19)

## Define Matrices for Predictors
# X0 : 7 X 38 matrix (7 smoking predictors for 38 control states)
X0 <- as.matrix(data[index_predict,index_co])

# X1 : 10 X 1 matrix (10 crime predictors for 1 treated states)
X1 <- as.matrix(data[index_predict,index_tr])

# [X1,X0]
X <- cbind(X1, X0)

## Define Matrices for Outcome Data
# Y0 : 31 X 38 matrix (31 years of smoking data for 38 control states)
Y0 <- as.matrix(data[index_Y,index_co])
# Y1 : 31 X 1 matrix (31 years of smoking data for 1 treated state)
Y1 <- as.matrix(data[index_Y,index_tr])

# [Y1,Y0]
Y <- cbind(Y1, Y0)

# Now pick Z matrices, i.e. the pretreatment period
# over which the loss function should be minmized
# Here we pick Z to go from 1970 to 1988 
# Z0 : 19 X 38 matrix (19 years of pre-treatment smoking data for 38 control states)
Z0 <- as.matrix(Y0[index_pre,])
# Z1 : 19 X 1 matrix (19 years of pre-treatment smoking data for 1 treated state)
Z1 <- as.matrix(Y1[index_pre,1])

# [Z1,Z0]
Z <- cbind(Z1, Z0)




W <- matrix(0, 1, dim(Y)[2])

W[1] <- 1


tune_params_synth <- list('Y' = Y, 'Z' = Z, 'X' = X)

#run synth so we can put it in the tables.
#Later I'll add these to general code.

smok_synth <- general_estimate(Y, Z, X, W, method = 'synth', tune_params = tune_params_synth)

smok_synth_tau <- smok_synth$Y_true - smok_synth$Y_est

#1990 should be 29?
#1995 should be 34
#2000 should be 39

##################
#### Boatlift ####
##################

dat <- readMat("boat_data.mat")
X <- dat$X
Y <- dat$Y
Z <- dat$Z

W <- matrix(0, 1, dim(Y)[2])

W(1) <- 1


#run synth so we can put it in the tables.
#Later I'll add these to general code.

boat_synth <- general_estimate(Y, Z, X, W, method = 'synth', tune_params = list("Y" = Y, "Z" = Z, "X" = X))


boat_synth_tau <- boat_synth$Y_true - boat_synth$Y_est
