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


source('functions.R')

vw_old <- tuning_parameters_synth(Y, Z, X, 1)

#vw <- c(1, 0.472278279279498, 1.42088591771919, 12501611.8714282, 0.0990703167049116)


vw_id <- c(1,1,1,1,1,1,1)

w <- find_weights_synth(Y, Z, X, vweight = vw_id)


ind_treat <- 1

synth.out <- 
  synth(
    X1 = as.matrix(X[, ind_treat]),
    X0 = as.matrix(X[, -ind_treat]),
    Z1 = as.matrix(Z[, ind_treat]),
    Z0 = as.matrix(Z[, -ind_treat]),
    
  )





##################
#### Boatlift ####
##################

dat <- readMat("boat_data.mat")
X <- dat$X
Y <- dat$Y
Z <- dat$Z

fit_elastic_net <- general_estimate(X = X, Y = Y, Z = Z,method = "elastic_net", 
                                    tune_params = list(c(seq(from = 1e-04, to = 1e-03, by = 1e-04),
                                                         seq(from = 2e-03, to = 1e-02, by = 1e-03),
                                                         seq(from = 2e-02, to = 1e-01, by = 1e-02),
                                                         seq(from = 2e-01, to = 100, by = 1e-01), 
                                                         seq(from = 200, to = 50000, by = 100)), seq(from = 0.1, to = 0.9, by = 0.1)),
                                    ind_treatment = 1)


source('functions.R')


boat_synth <- 
  synth(
    X1 = as.matrix(X[, ind_treat]),
    X0 = as.matrix(X[, -ind_treat]),
    Z1 = as.matrix(Z[, ind_treat]),
    Z0 = as.matrix(Z[, -ind_treat]),
    quadopt = "LowRankQP"
  )


vw_old <- tuning_parameters_synth(Y, Z, X, 1)

vw <- c(1, 0.472278279279498, 1.42088591771919, 12501611.8714282, 0.0990703167049116)


#the big difference is that there are way more non-zero weights

#when we just use our way of doing weights, the results are muuuuch closer. the zeros are not as small, however.

#yeah when we use our procedure with our weights

w_check <- synth(
  X1 = as.matrix(X[, ind_treat]),
  X0 = as.matrix(X[, -ind_treat]),
  Z1 = as.matrix(Z[, ind_treat]),
  Z0 = as.matrix(Z[, -ind_treat]),
  vweight = boat_synth$solution.v,
)

w <- find_weights_synth(Y, Z, X, vweight = boat_synth$solution.v)


#goal: look at synthetic control in Matlab and figure out what is going on