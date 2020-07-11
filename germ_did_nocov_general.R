#this will contain general functions for diff-in-diff




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
library(Synth)



#hmmm I might need to change how we do the other ones a little bit? for diff in diff it makes sense to do dataprep
#beforehand but for synth it makes sense to do it within...I'll talk to Lea about it?

prep_data <- function(d, pred, dep, u, t, spec,i, j, subs, year1, year2, year3, year4, year5, year6, names ){
  
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
      time.predictors.prior = year1:year2,
      time.optimize.ssr = year3:year4,
      unit.names.variable = names,
      time.plot = year5:year6
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
  names(datlist) <- c( "x", "y", "z")
  output<-datlist
  }

find_weights_did <- function(data){
  N <- dim(data$y)[2]
  w_did <- matrix(1 / (N - 1), nrow = N - 1, ncol = 1)  
}


find_y_did <- function(data, weights, i ,T){
 
  int_did <- as.matrix(mean(data$z[,i]) - mean(data$z[,-i]))
  
  Y_did <- int_did[rep(1, T),] + data$y[,-i] %*% weights
  output <- list(int <- int_did, Yhat <- Y_did)
  names(output) <- c("intercept", "Y_did")
  output<- output
}

find_treatment <- function(data, weights, T){
  
  treats <- data$Y1 - find_y_did(data, weights, T)
}




se_units_did <- function(data, T0, T, j){
  N <- dim(data$y)[2]

  data_trunc <- list(y<-data$y[,-j], z<-data$z[,-j], x<-data$x[,-j])
  names(data_trunc) <- c("y", "z", "x")
  
  w <- find_weights_did(data_trunc)

  std_err_i <- matrix(0, T1, N - 1)

  for (i in 1:(N-1)) {
    Y1 <- as.matrix(data_trunc$y[,i])
    y_did <- find_y_did(data_trunc, w, i, T)
  
    std_err_i[,i] <- (Y1[-c(1:T0),] - y_did$Y_did[-c(1:T0),]) ^ 2
  
  }

  std_err_i <- as.matrix(sqrt(apply(std_err_i, 1, mean)))
}

#I can finish this at Chloey's I guess
se_time_did <- function(data, T0, T, j){
  
  

}

## Load the data
#load("germany_data.RData")
load("repgermany.Rdata")
d <- x

predictors <- c("gdp","trade","infrate")
special = list(
  list("industry" ,1981:1990, c("mean")),
  list("schooling",c(1980,1985), c("mean")),
  list("invest80" ,1980, c("mean"))
)

year <- c(1981, 1990, 1960, 1989, 1960, 2003)
  
data <- prep_data(d, predictors, "gdp", 1, 3, special, 7, 7, unique(d$index), year[1], year[2], year[3], year[4], year[5], year[6], 2)

#wait why can't i use dollar signs anymore?????

#hmmmmm do we want to just write a function that finds parameters to plug into other functions? that might be useful.

#number of covariates
K <- dim(data$x)[1]

#number of units
N <- dim(data$y)[2]

#number of periods
T <- dim(data$y)[1]

#number of pretreatment periods
T0 <- dim(data$z)[1]
#number of post-treatment peiords
T1 <- T - T0
#not sure what these are. They don't ever get used as far as I can tell.
T0_tr <- floor(T0 * 2 / 3)
T0_te <- T0 - T0_tr

w <- find_weights_did(data)

#okay so the main difference between our y0 (from the dataprep i guess?) and their y0 is this one only goes back to 1989.
#I guess that means we do ahve to use y instead. good to know.

#this bit is very different right now and we need tow work on it.

#ask Lea why naming directly works sometimes but not other times.
output <- find_y_did(data, w, 1, T)

Y_did <- output$Y_did
intercept <- output$intercept

treatment <- find_treatment(data, w, T)

std_err_i <- se_units_did(data, T0, T, 1)




  
  
  
  
  
  
  
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

