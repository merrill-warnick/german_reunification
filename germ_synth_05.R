## Replication Code for
# A. Abadie, A. Diamond, and J. Hainmueller. 2014.
# Comparative Politics and the Synthetic Control Method
# American Journal of Political Science.

rm(list=ls())
library(foreign)
library(Synth)
library(xtable)
library(R.matlab)


######################
####### DATA #########
######################

# Load Data 
d <- read.dta("repgermany.dta")

## Table 1 & 2, Figure 1, 2, & 3


########################
####### PICK V #########
########################

# v determines the relative importance of each predictor

## pick v by cross-validation

#########
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
            list("industry", 1971:1980, c("mean")), # difference here
            list("schooling",c(1970,1975), c("mean")), # difference here
            list("invest70" ,1980, c("mean"))
           ),
           treatment.identifier = 7,
           controls.identifier = unique(d$index)[-7],
           time.predictors.prior = 1971:1980,
           time.optimize.ssr = 1981:1990, # this is where difference is!
           unit.names.variable = 2,
           time.plot = 1960:2003
         )


####################################
####### FIT TRAINING MODEL #########
####################################

# fit training model
synth.out <- 
  synth(
        data.prep.obj=dataprep.out,
        Margin.ipop=.005,Sigf.ipop=7,Bound.ipop=6
        ) # use pre-specified function for synthetic control

############################
####### MAIN MODEL #########
############################

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

################################
####### FIT MAIN MODEL #########
################################

# Use optimal V from training model!
# fit main model with v from training model
synth.out <- synth(
  data.prep.obj=dataprep.out,
  custom.v=as.numeric(synth.out$solution.v)
  ) # Built-in Synth function -> custom v applied

#########################
####### RESULTS #########
#########################

w_synth <- synth.out$solution.w
Y_synth <- dataprep.out$Y0 %*% w_synth # Estimated Y (no treatment)

#################################
####### Standard Errors #########
#################################

## Define parameters first (manually....why? To get it in framework we could do X,Y,Z again)
N <- 17
T <- 44
T0 <- 30
T1 <- T - T0
units_co <- c(1, 2, 3, 4, 5, 6, 8, 9, 10, 12, 14, 16, 18, 19, 20, 21)

# Over units
std_err_i <- matrix(0, N - 1, T1) # Storage matrix

for (j in 1:(N - 1)) {
  i <- units_co[j]
  cat('(Std. Error) Over Unit i =', toString(i), '\n')
  
  ## STEP 1: Pick V each time (why???? In other cases also never pick new alpha e.g.)
  
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
      treatment.identifier = i, # This is changed in comparison to above
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
  
  ## STEP 2: Fit main model
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
  
  ## STEP 3: Save solutions
  # Solutions
  w <- synth.out$solution.w
  Y_pred <- dataprep.out$Y0 %*% w
  Y1 <- dataprep.out$Y1
  std_err_i[j,] <- (Y1[-c(1:T0),] - Y_pred[-c(1:T0),]) ^ 2 # SST for each unit
}
std_err_i <- as.matrix(sqrt(apply(std_err_i, 2, mean))) # Sqrt of mean of SSR over each unit gives std err for each time point

# Over time
s <- floor(T0 / 2) # Number of time periods
std_err_t <- matrix(0, s, 1) # Storage matrix

for (t in 1:s) {
  cat('(Std. Error) Over Time t =', toString(t), '\n')
  
  ## STEP 1: Pick V
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
  
  ## STEP 2: Fit Main Model
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
      time.optimize.ssr = 1960:(1989 - t), # HERE is the change to general set-up!
      unit.names.variable = 2,
      time.plot = 1960:2003
    )
  
  # fit main model with v from training model
  synth.out <- synth(
    data.prep.obj=dataprep.out,
    custom.v=as.numeric(synth.out$solution.v)
  )
  
  ## STEP 3: Save solutions
  # Solutions
  w <- synth.out$solution.w
  Y_pred <- dataprep.out$Y0 %*% w
  Y1 <- dataprep.out$Y1
  std_err_t[t,1] <- (Y1[T0 - t + 1,] - Y_pred[T0 - t + 1,]) ^ 2 # SSR over time
}
std_err_t <- as.matrix(sqrt(apply(std_err_t, 2, mean))) # Sqrt of mean over time gives std err

# Over units and time
std_err_it <- matrix(0, N - 1, 1) # Storage matrix

for (j in 1:(N - 1)) {
  i <- units_co[j]
  std_err_temp <- matrix(0, s, 1) # Storage for temporary std err over time for specific unit
  
  for (t in 1:s) {
    cat('(Std. Error) Over Unit and Time ( i , t ) = (',toString(i), ',', 
        toString(t), ')\n')
    
    # STEP 1: Pick V via cross validation
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
        treatment.identifier = i, # Fix this here from unit loop
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
    
    ## STEP 2: Fit Main Model
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
        treatment.identifier = i, # HERE change
        controls.identifier = units_co[-j],
        time.predictors.prior = 1981:1990,
        time.optimize.ssr = 1960:(1989 - t), # And HERE change
        unit.names.variable = 2,
        time.plot = 1960:2003
      )
    
    # fit main model with v from training model
    synth.out <- synth(
      data.prep.obj=dataprep.out,
      custom.v=as.numeric(synth.out$solution.v)
    )
    
    ## STEP 3: Save solutions
    # Solutions
    w <- synth.out$solution.w
    Y_pred <- dataprep.out$Y0 %*% w
    Y1 <- dataprep.out$Y1
    std_err_temp[t,1] <- (Y1[T0 - t + 1,] - Y_pred[T0 - t + 1,]) ^ 2
  }
  std_err_temp <- as.matrix(apply(std_err_temp, 2, mean))
  std_err_it[j,1] <- std_err_temp # Get std err over time for each unit
}
std_err_it <- as.matrix(sqrt(apply(std_err_it, 2, mean))) # Sqrt of mean over units to get std err

# Copy the standard errors
std_err_synth_i <- std_err_i
std_err_synth_t <- std_err_t
std_err_synth_it <- std_err_it

################################
####### Counterfactual #########
################################

# Pretend intervention happened in T0 = 21 and not T0 = 30

#######
# Optimal V

## pick v by cross-validation
# data setup for training model (counterfactual)
cat('*** Counterfactual ***\n')

## No difference here to earlier??

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

#######
# fit training model
synth.out <- 
  synth(
    data.prep.obj=dataprep.out,
    Margin.ipop=.005,Sigf.ipop=7,Bound.ipop=6
  )
#######
# MAIN MODEL


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
    time.optimize.ssr = 1960:1980, # difference here!!! Only until 1980!!
    unit.names.variable = 2,
    time.plot = 1960:1989
  )

######
# fit main model with v from training model

synth.out <- synth(
  data.prep.obj=dataprep.out,
  custom.v=as.numeric(synth.out$solution.v)
)

#########################
####### RESULTS #########
#########################

w_synth_co <- synth.out$solution.w
Y_synth_co <- dataprep.out$Y0 %*% w_synth_co # Estimated fit counterfactual (no treatment)
Y_true_co <- dataprep.out$Y1


#################################
####### Standard Errors #########
#################################

######
## Compute the counterfactual standard errors
N <- 17
T0_co <- 21
T1_co <- T0 - T0_co
units_co <- c(1, 2, 3, 4, 5, 6, 8, 9, 10, 12, 14, 16, 18, 19, 20, 21)

# Over units
std_err_i_co <- matrix(0, N - 1, T1_co) # Storage matrix

for (j in 1:(N - 1)) {
  i <- units_co[j]
  cat('(Std. Error) Over Unit i =', toString(i), '\n')
  
  ## STEP 1: Pick V by cross-validation
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
  
  ## STEP 2: Fit Main Model
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
      time.optimize.ssr = 1960:1980, # Still only until 1980 here!
      unit.names.variable = 2,
      time.plot = 1960:1989
    )
  
  # fit main model with v from training model
  synth.out <- synth(
    data.prep.obj=dataprep.out,
    custom.v=as.numeric(synth.out$solution.v)
  )
  
  ## STEP 3: Save results
  # Solutions
  w <- synth.out$solution.w
  Y_pred <- dataprep.out$Y0 %*% w
  Y1 <- dataprep.out$Y1
  std_err_i_co[j,] <- (Y1[-c(1:T0_co),] - Y_pred[-c(1:T0_co),]) ^ 2 # SSR for each unit and time point
}
std_err_i_co <- as.matrix(sqrt(apply(std_err_i_co, 2, mean))) # Sqrt of mean over time to get unit std err

# Copy the counterfactual standard errors
std_err_synth_i_co <- std_err_i_co

##############################
####### Save Results #########
##############################

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
