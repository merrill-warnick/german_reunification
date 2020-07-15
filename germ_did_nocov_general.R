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

#prep data needs tochane anyway? I guess I can just deal with that right now as part of the restructure.

#so prep data will just take in the data and spit out either the full structure to be used for synth or just Y X Z for other stuff
#I gotta delve deeper and figure out what turns into Y and what turns into Z.
#We might need to change this when we start adding covariates I guess?

#actually now that I think about it, the output for prep data will be the same unless full = 1.

prep_data <- function(d, pred, dep, u, t, spec,i, j, subs, year1, year2, year3, year4, year5, year6, names, full ){
  
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
  
  #so I think that if full==FALSE, it willjust output the dataprep.out object.
  #I guess that this will make it so that in our data, the true dataset will be in the first column.
  #we should still write component functions to take in the input of which column is treated because that will be important to change
  #when calculating standard errors.
  if (full==FALSE){
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
  
  }

find_weights_did <- function(Y, Z, X){
  N <- dim(Y)[2]
  w_did <- matrix(1 / (N - 1), nrow = N - 1, ncol = 1)  
}


find_y_did <- function(Y, Z, X, weights, i ,T){
  #cat(mean(Z[,i]))
  int_did <- as.matrix(mean(Z[,i]) - mean(Z[,-i]))
  
  Y_did <- int_did[rep(1, T),] + Y[,-i] %*% weights
  output <- list(int <- int_did, Yhat <- Y_did)
  names(output) <- c("intercept", "Y_did")
  output<- output
}

find_treatment <- function(Y, Z, X, weights, i, T){
  #I guess the treatment period is T0+1, but on second thought we just give thema ll the treatment effects
  treats <- Y[,i] - find_y_did(Y, Z, X, weights, i, T)$Y_did
}

se_units_did <- function(Y, Z, X, i, T0, T){
  N <- dim(data$y)[2]
  
  Y_trunc <- Y[,-i]
  Z_trunc <- Z[,-i]
  X_trunc <- X[,-i]
  
  w <- find_weights_did(Y_trunc, Z_trunc, X_trunc)

  std_err_i <- matrix(0, T1, N - 1)

  for (j in 1:(N-1)) {
    
    y_did <- find_y_did(Y_trunc, Z_trunc, X_trunc, w, j ,T)
    
    std_err_i[,j] <- (Y_trunc[-c(1:T0),j] - y_did$Y_did[-c(1:T0),]) ^ 2
  
  }

  std_err_i <- as.matrix(sqrt(apply(std_err_i, 1, mean)))
}


se_time_did <- function(Y, Z, X, i, T0, T){
  
  N <- dim(Y)[2]
  w <- find_weights_did(Y, Z, X)
  
  # Over time
  s <- floor(T0 / 2)
  std_err_t <- matrix(0, s, 1)
  
  #cat(s)
  for (t in 1:s) {
    #
    #okay so we should pass in truncated versions I think? I'm not sure actually, it seems a little weird, need to investigate.
    #After looking at this again, the prep_data part makes it so that the Zs are only for pretreatment periods. That means that I guess that we 
    #*do* want to just pass in the truncated Zs. Passing in the full Ys should be fine, evaluate Xs later.
    # Fit did
    
    y_did <- find_y_did(Y, Z[c(1:(T0-t)),], X, w, i ,T)
    
    #I think that this will work
    std_err_t[t,1] <- (Y[T0 - t + 1, i] - y_did$Y_did[T0 - t + 1,]) ^ 2
  }
  std_err_t <- as.matrix(sqrt(apply(std_err_t, 2, mean)))

}


se_units_time_did <- function(Y, Z, X, i, T0, T){
  N <- dim(Y)[2]
  
  s <- floor(T0 / 2)
  std_err_it <- matrix(0, N - 1, 1)
  
  Y_trunc <- Y[,-i]
  Z_trunc <- Z[,-i]
  X_trunc <- X[,-i]
  
  w <- find_weights_did(Y_trunc, Z_trunc, X_trunc)
  
  for (j in 1:(N-1)){
    std_err_temp <- matrix(0, s, 1)
    for (t in 1:s){
      y_did <- find_y_did(Y_trunc, Z_trunc[c(1:(T0-t)),], X_trunc, w, j ,T)
      
      std_err_temp[t,1] <- (Y_trunc[T0 - t + 1, j] - y_did$Y_did[T0 - t + 1,]) ^ 2
    }
    std_err_temp <- as.matrix(apply(std_err_temp, 2, mean))
    std_err_it[j,1] <- std_err_temp
  }
  std_err_it <- as.matrix(sqrt(apply(std_err_it, 2, mean)))
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
  
data <- prep_data(d, predictors, "gdp", 1, 3, special, 7, 7, unique(d$index), year[1], year[2], year[3], year[4], year[5], year[6], 2, FALSE)

Y <- data$y
Z <- data$z
X <- data$x


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

w <- find_weights_did(Y, Z, X)

#okay so the main difference between our y0 (from the dataprep i guess?) and their y0 is this one only goes back to 1989.
#I guess that means we do ahve to use y instead. good to know.

#this bit is very different right now and we need tow work on it.

#ask Lea why naming directly works sometimes but not other times.
output <- find_y_did(Y, Z, X, w, 1, T)

Y_did <- output$Y_did


intercept <- output$intercept

treatment <- find_treatment(Y, Z, X, w, 1, T)

std_err_i <- se_units_did(Y, Z, X, 1, T0, T)

std_err_t <- se_time_did(Y, Z, X, 1, T0, T)

std_err_it <- se_units_time_did(Y, Z, X, 1, T0, T)

  
#okayt this all waorks. Need to figure out the next part I guess.

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

