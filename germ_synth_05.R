## Replication Code for
# A. Abadie, A. Diamond, and J. Hainmueller. 2014.
# Comparative Politics and the Synthetic Control Method
# American Journal of Political Science.

rm(list=ls())
library(foreign)
library(Synth)
library(xtable)
library(R.matlab)

# Load Data 
#setwd("/Users/nikolayd/Dropbox/Research/Synth/Germany/")
d <- read.dta("repgermany.dta")

## Table 1 & 2, Figure 1, 2, & 3

## pick v by cross-validation
# data setup for training model
cat('*** Main ***\n')
dataprep.out <-
  dataprep(
           foo = d,
           predictors    = c("gdp","trade","infrate"),
           dependent     = "gdp",
           unit.variable = 1,
           time.variable = 3,
           special.predictors = list(
            list("industry", 1971:1980, c("mean")),
            list("schooling",c(1970,1975), c("mean")),
            list("invest70" ,1980, c("mean"))
           ),
           treatment.identifier = 7,
           controls.identifier = unique(d$index)[-7],
           time.predictors.prior = 1971:1980,
           time.optimize.ssr = 1981:1990,
           unit.names.variable = 2,
           time.plot = 1960:2003
         )

# fit training model
synth.out <- 
  synth(
        data.prep.obj=dataprep.out,
        Margin.ipop=.005,Sigf.ipop=7,Bound.ipop=6
        )

# data prep for main model
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

# fit main model with v from training model
synth.out <- synth(
  data.prep.obj=dataprep.out,
  custom.v=as.numeric(synth.out$solution.v)
  )

#### Main results
w_synth <- synth.out$solution.w
Y_synth <- dataprep.out$Y0 %*% w_synth

## Compute the standard errors
N <- 17
T <- 44
T0 <- 30
T1 <- T - T0
units_co <- c(1, 2, 3, 4, 5, 6, 8, 9, 10, 12, 14, 16, 18, 19, 20, 21)
# Over units
std_err_i <- matrix(0, N - 1, T1)
for (j in 1:(N - 1)) {
  i <- units_co[j]
  cat('(Std. Error) Over Unit i =', toString(i), '\n')
  ## pick v by cross-validation
  # data setup for training model
  dataprep.out <-
    dataprep(
      foo = d,
      predictors    = c("gdp","trade","infrate"),
      dependent     = "gdp",
      unit.variable = 1,
      time.variable = 3,
      special.predictors = list(
        list("industry", 1971:1980, c("mean")),
        list("schooling",c(1970,1975), c("mean")),
        list("invest70" ,1980, c("mean"))
      ),
      treatment.identifier = i,
      controls.identifier = units_co[-j],
      time.predictors.prior = 1971:1980,
      time.optimize.ssr = 1981:1990,
      unit.names.variable = 2,
      time.plot = 1960:2003
    )
  
  # fit training model
  synth.out <- 
    synth(
      data.prep.obj=dataprep.out,
      Margin.ipop=.005,Sigf.ipop=7,Bound.ipop=6
    )
  
  # data prep for main model
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
      treatment.identifier = i,
      controls.identifier = units_co[-j],
      time.predictors.prior = 1981:1990,
      time.optimize.ssr = 1960:1989,
      unit.names.variable = 2,
      time.plot = 1960:2003
    )
  
  # fit main model with v from training model
  synth.out <- synth(
    data.prep.obj=dataprep.out,
    custom.v=as.numeric(synth.out$solution.v)
  )
  
  # Solutions
  w <- synth.out$solution.w
  Y_pred <- dataprep.out$Y0 %*% w
  Y1 <- dataprep.out$Y1
  std_err_i[j,] <- (Y1[-c(1:T0),] - Y_pred[-c(1:T0),]) ^ 2
}
std_err_i <- as.matrix(sqrt(apply(std_err_i, 2, mean)))

# Over time
s <- floor(T0 / 2)
std_err_t <- matrix(0, s, 1)
for (t in 1:s) {
  cat('(Std. Error) Over Time t =', toString(t), '\n')
  ## pick v by cross-validation
  # data setup for training model
  dataprep.out <-
    dataprep(
      foo = d,
      predictors    = c("gdp","trade","infrate"),
      dependent     = "gdp",
      unit.variable = 1,
      time.variable = 3,
      special.predictors = list(
        list("industry", 1971:1980, c("mean")),
        list("schooling",c(1970,1975), c("mean")),
        list("invest70" ,1980, c("mean"))
      ),
      treatment.identifier = 7,
      controls.identifier = unique(d$index)[-7],
      time.predictors.prior = 1971:1980,
      time.optimize.ssr = 1981:1990,
      unit.names.variable = 2,
      time.plot = 1960:2003
    )
  
  # fit training model
  synth.out <- 
    synth(
      data.prep.obj=dataprep.out,
      Margin.ipop=.005,Sigf.ipop=7,Bound.ipop=6
    )
  
  # data prep for main model
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
      time.optimize.ssr = 1960:(1989 - t),
      unit.names.variable = 2,
      time.plot = 1960:2003
    )
  
  # fit main model with v from training model
  synth.out <- synth(
    data.prep.obj=dataprep.out,
    custom.v=as.numeric(synth.out$solution.v)
  )
  
  # Solutions
  w <- synth.out$solution.w
  Y_pred <- dataprep.out$Y0 %*% w
  Y1 <- dataprep.out$Y1
  std_err_t[t,1] <- (Y1[T0 - t + 1,] - Y_pred[T0 - t + 1,]) ^ 2
}
std_err_t <- as.matrix(sqrt(apply(std_err_t, 2, mean)))

# Over units and time
std_err_it <- matrix(0, N - 1, 1)
for (j in 1:(N - 1)) {
  i <- units_co[j]
  std_err_temp <- matrix(0, s, 1)
  for (t in 1:s) {
    cat('(Std. Error) Over Unit and Time ( i , t ) = (',toString(i), ',', 
        toString(t), ')\n')
    ## pick v by cross-validation
    # data setup for training model
    dataprep.out <-
      dataprep(
        foo = d,
        predictors    = c("gdp","trade","infrate"),
        dependent     = "gdp",
        unit.variable = 1,
        time.variable = 3,
        special.predictors = list(
          list("industry", 1971:1980, c("mean")),
          list("schooling",c(1970,1975), c("mean")),
          list("invest70" ,1980, c("mean"))
        ),
        treatment.identifier = i,
        controls.identifier = units_co[-j],
        time.predictors.prior = 1971:1980,
        time.optimize.ssr = 1981:1990,
        unit.names.variable = 2,
        time.plot = 1960:2003
      )
    
    # fit training model
    synth.out <- 
      synth(
        data.prep.obj=dataprep.out,
        Margin.ipop=.005,Sigf.ipop=7,Bound.ipop=6
      )
    
    # data prep for main model
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
        treatment.identifier = i,
        controls.identifier = units_co[-j],
        time.predictors.prior = 1981:1990,
        time.optimize.ssr = 1960:(1989 - t),
        unit.names.variable = 2,
        time.plot = 1960:2003
      )
    
    # fit main model with v from training model
    synth.out <- synth(
      data.prep.obj=dataprep.out,
      custom.v=as.numeric(synth.out$solution.v)
    )
    
    # Solutions
    w <- synth.out$solution.w
    Y_pred <- dataprep.out$Y0 %*% w
    Y1 <- dataprep.out$Y1
    std_err_temp[t,1] <- (Y1[T0 - t + 1,] - Y_pred[T0 - t + 1,]) ^ 2
  }
  std_err_temp <- as.matrix(apply(std_err_temp, 2, mean))
  std_err_it[j,1] <- std_err_temp
}
std_err_it <- as.matrix(sqrt(apply(std_err_it, 2, mean)))

# Copy the standard errors
std_err_synth_i <- std_err_i
std_err_synth_t <- std_err_t
std_err_synth_it <- std_err_it

## pick v by cross-validation
# data setup for training model (counterfactual)
cat('*** Counterfactual ***\n')
dataprep.out <-
  dataprep(
    foo = d,
    predictors    = c("gdp","trade","infrate"),
    dependent     = "gdp",
    unit.variable = 1,
    time.variable = 3,
    special.predictors = list(
      list("industry", 1971:1980, c("mean")),
      list("schooling",c(1970,1975), c("mean")),
      list("invest70" ,1980, c("mean"))
    ),
    treatment.identifier = 7,
    controls.identifier = unique(d$index)[-7],
    time.predictors.prior = 1971:1980,
    time.optimize.ssr = 1981:1990,
    unit.names.variable = 2,
    time.plot = 1960:2003
  )

# fit training model
synth.out <- 
  synth(
    data.prep.obj=dataprep.out,
    Margin.ipop=.005,Sigf.ipop=7,Bound.ipop=6
  )

# data prep for main model
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
    time.optimize.ssr = 1960:1980,
    unit.names.variable = 2,
    time.plot = 1960:1989
  )

# fit main model with v from training model
synth.out <- synth(
  data.prep.obj=dataprep.out,
  custom.v=as.numeric(synth.out$solution.v)
)

# Save the counterfactual results
w_synth_co <- synth.out$solution.w
Y_synth_co <- dataprep.out$Y0 %*% w_synth_co
Y_true_co <- dataprep.out$Y1

## Compute the counterfactual standard errors
N <- 17
T0_co <- 21
T1_co <- T0 - T0_co
units_co <- c(1, 2, 3, 4, 5, 6, 8, 9, 10, 12, 14, 16, 18, 19, 20, 21)
# Over units
std_err_i_co <- matrix(0, N - 1, T1_co)
for (j in 1:(N - 1)) {
  i <- units_co[j]
  cat('(Std. Error) Over Unit i =', toString(i), '\n')
  ## pick v by cross-validation
  # data setup for training model
  dataprep.out <-
    dataprep(
      foo = d,
      predictors    = c("gdp","trade","infrate"),
      dependent     = "gdp",
      unit.variable = 1,
      time.variable = 3,
      special.predictors = list(
        list("industry", 1971:1980, c("mean")),
        list("schooling",c(1970,1975), c("mean")),
        list("invest70" ,1980, c("mean"))
      ),
      treatment.identifier = i,
      controls.identifier = units_co[-j],
      time.predictors.prior = 1971:1980,
      time.optimize.ssr = 1981:1990,
      unit.names.variable = 2,
      time.plot = 1960:2003
    )
  
  # fit training model
  synth.out <- 
    synth(
      data.prep.obj=dataprep.out,
      Margin.ipop=.005,Sigf.ipop=7,Bound.ipop=6
    )
  
  # data prep for main model
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
      treatment.identifier = i,
      controls.identifier = units_co[-j],
      time.predictors.prior = 1981:1990,
      time.optimize.ssr = 1960:1980,
      unit.names.variable = 2,
      time.plot = 1960:1989
    )
  
  # fit main model with v from training model
  synth.out <- synth(
    data.prep.obj=dataprep.out,
    custom.v=as.numeric(synth.out$solution.v)
  )
  
  # Solutions
  w <- synth.out$solution.w
  Y_pred <- dataprep.out$Y0 %*% w
  Y1 <- dataprep.out$Y1
  std_err_i_co[j,] <- (Y1[-c(1:T0_co),] - Y_pred[-c(1:T0_co),]) ^ 2
}
std_err_i_co <- as.matrix(sqrt(apply(std_err_i_co, 2, mean)))

# Copy the counterfactual standard errors
std_err_synth_i_co <- std_err_i_co

## Save the results
save(list = c("w_synth", "Y_synth", 
              "std_err_synth_i", "std_err_synth_t", "std_err_synth_it", 
              "Y_synth_co", "Y_true_co", "std_err_synth_i_co"), 
     file = "germ_synth_05.RData")
writeMat("germ_synth_05.mat", 
         w_synth = w_synth, 
         Y_synth = Y_synth,
         std_err_synth_i = std_err_i, 
         std_err_synth_t = std_err_t, 
         std_err_synth_it = std_err_it,
         Y_synth_co = Y_synth_co,
         Y_true_co = Y_true_co,
         std_err_synth_i_co = std_err_i_co)
