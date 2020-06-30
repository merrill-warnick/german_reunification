#The plan in this file is to do some work to prep cross-validation experimentation for vweights in the
#synthetic control part of the project. My plan is to first write a function that finds the vweights generally?
#now I guess the question is, do I want to try to make it whole hog general now? Yeah I guess so, that makes sense, esp since
#during variance calculations I need to assign different treatment units/years.


#this is all the code I need for it, so I just need to go through and make the argument general


rm(list=ls())
library(foreign)
library(Synth)
library(xtable)
library(R.matlab)


find_vweights <- function(d, pred, y, u, t, spec, i,j, predyear0, predyear1, optyear0, optyear1, names, year0, year1){
  #d is the dataframe of the panel data
  #pred is a string of predictor variables
  #y is the string name of the dependent variable
  #u is the value of the unit identifier column
  #t is the value of the time identifier column 
  #spec is a list of special predictors that should be what you're plugging into special.predictors
  #  for the dataprep function
  #i is the index of the treatment identifier column
  #j helps with the index of the control identifier columns in case you need it
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
    
    #I'm not sure what this does but I think it's just supposed to be -i
    #controls.identifier = unique(d$index)[-7],
    controls.identifier = unique(d$index)[-j],
    
    time.predictors.prior = predyear0:predyear1,
    
    time.optimize.ssr = optyear0:optyear1,
    
    unit.names.variable = names,
    
    time.plot = year0:year1
  )

#fit training model to pull out vweights
synth.out <- 
  synth(
    data.prep.obj=dataprep.out,
    
    #these are just optimization things, for now I'm just going to leave them as is, but I guess if I look at later stuff and see that we make different
    # choices, I'll let these be inputs.
    Margin.ipop=.005,Sigf.ipop=7,Bound.ipop=6
  )
output <- synth.out$solution.v
}









## Replication Code for
# A. Abadie, A. Diamond, and J. Hainmueller. 2014.
# Comparative Politics and the Synthetic Control Method
# American Journal of Political Science.


#first, we load in the German data
d <- read.dta("repgermany.dta")

## Table 1 & 2, Figure 1, 2, & 3

cat('*** Main ***\n')


special<-list(
  list("industry" ,1981:1990, c("mean")),
  list("schooling",c(1980,1985), c("mean")),
  list("invest80" ,1980, c("mean"))
)

#some change that we'll want to assign these differently later
#I thin kthat the first thing to do before I start investigating is to make sure that it runs and gives the same results under
#the new function

#while it's running, I might as well work on a small function that generalizes the whole process I guess.

#vweights
vw <- find_vweights(d, c("gdp","trade","infrate"), "gdp", 1, 3, special, 7,7, 1971, 1980, 1981, 1990, 2, 1960, 2003)

# data prep for main model
# why do we redo dataprep?
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
    
    #hmmmm why are the ten years before the only predictors but we optimize over all forty years?
    time.predictors.prior = 1981:1990,
    time.optimize.ssr = 1960:1989,
    unit.names.variable = 2,
    time.plot = 1960:2003
  )

# fit main model with v from training model
#this gives us what our synthetic control is I think.
synth.out <- synth(
  data.prep.obj=dataprep.out,
  #this is where we specify the vs that we got from the last piece
  custom.v=as.numeric(vw)
)

#### Main results
#this gives us our synthetic weights and our synthetic outcoms.
w_synth <- synth.out$solution.w
Y_synth <- dataprep.out$Y0 %*% w_synth

## Compute the standard errors
#I think we're doing this by doing the counterfactual exercise over units--pretend that the wrong unit is the treated unit etc
N <- 17
T <- 44
T0 <- 30
T1 <- T - T0
units_co <- c(1, 2, 3, 4, 5, 6, 8, 9, 10, 12, 14, 16, 18, 19, 20, 21)

#need to save vweight results to think about them. Actually we only really need to do it for the unitsXtime one since it
# will capture everything going on in the other parts.

# Over units
std_err_i <- matrix(0, N - 1, T1)
for (j in 1:(N - 1)) {
  i <- units_co[j]
  cat('(Std. Error) Over Unit i =', toString(i), '\n')
  ## pick v by cross-validation
  # data setup for training model
  
  vw <- find_vweights(d, c("gdp","trade","infrate"), "gdp", 1, 3, special, i,j, 1971, 1980, 1981, 1990, 2, 1960, 2003)
  
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
    custom.v=as.numeric(vw)
  )
  
  # Solutions
  w <- synth.out$solution.w
  Y_pred <- dataprep.out$Y0 %*% w
  Y1 <- dataprep.out$Y1
  #I don't quite understand the indexing but I see what is going on here, I think.
  std_err_i[j,] <- (Y1[-c(1:T0),] - Y_pred[-c(1:T0),]) ^ 2
}
std_err_i <- as.matrix(sqrt(apply(std_err_i, 2, mean)))

# Over time
#calculating standard errors again, but we're assigning treatment to the wrong year now.
s <- floor(T0 / 2)
std_err_t <- matrix(0, s, 1)
for (t in 1:s) {
  cat('(Std. Error) Over Time t =', toString(t), '\n')
  ## pick v by cross-validation
  # data setup for training model
  
  #question: why don't we change how we do crossvalidation when we change which year we're treating as treated? We need to change that I think.
  vw <- find_vweights(d, c("gdp","trade","infrate"), "gdp", 1, 3, special, 7,7, 1971, 1980, 1981, 1990, 2, 1960, 2003)
  
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
      #waaaait I see that we're optimizing over the right set of years, but why is the set of time predictors unchanged?
      #I guess after reading the documentation I'm not 100% sure what time.predictors.prior is doing, and I'll need to talk to somebody about it I think.
      time.predictors.prior = 1981:1990,
      time.optimize.ssr = 1960:(1989 - t),
      unit.names.variable = 2,
      time.plot = 1960:2003
    )
  
  # fit main model with v from training model
  synth.out <- synth(
    data.prep.obj=dataprep.out,
    custom.v=as.numeric(vw)
  )
  
  # Solutions
  w <- synth.out$solution.w
  Y_pred <- dataprep.out$Y0 %*% w
  Y1 <- dataprep.out$Y1
  std_err_t[t,1] <- (Y1[T0 - t + 1,] - Y_pred[T0 - t + 1,]) ^ 2
}
std_err_t <- as.matrix(sqrt(apply(std_err_t, 2, mean)))

# Over units and time

#need to save vweights

#not sure where we assign s
#also not sure how big vw is oging to be
vweights_it <-array(dim=c(N-1, s, length(vw)))

std_err_it <- matrix(0, N - 1, 1)
for (j in 1:(N - 1)) {
  i <- units_co[j]
  std_err_temp <- matrix(0, s, 1)
  for (t in 1:s) {
    cat('(Std. Error) Over Unit and Time ( i , t ) = (',toString(i), ',', 
        toString(t), ')\n')
    ## pick v by cross-validation
    
    #ohhhh I guess I need to think about i and j better real quick
    vw <- find_vweights(d, c("gdp","trade","infrate"), "gdp", 1, 3, special, i, j, 1971, 1980, 1981, 1990, 2, 1960, 2003)
    
    #save vweights so we can look at them
    vweights_it(i,t,) <- vw
    
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
      custom.v=as.numeric(vw)
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
#gotta find where these standard errors show up in the paper. Hmmm it looks like they don't, I think that the
#paper only talks about p-values. well, where do these show up? maybe in the original abadie cali smoking paper? Yeah I don't see it.
#I guess we should talk about this. Maybe just guido calculates this? yeah okay it looks like that these standard errors are not an
#abadie thing but a guido thing
#yeah I reread the paper and this is definitely from our paper
std_err_synth_i <- std_err_i
std_err_synth_t <- std_err_t
std_err_synth_it <- std_err_it

## pick v by cross-validation
# data setup for training model (counterfactual)

#okay, what is this section doing?
#tbh I'm not totally clear on what the counterfactual exercise is doing in the paper
#or rather I'm not sure how it works, maybe I can read the code and figure it out.

#this counterfactual section seems redunant, I think it's recalculating stuff we did in previous loops--the across units and time bit.

cat('*** Counterfactual ***\n')

vw <- find_vweights(d, c("gdp","trade","infrate"), "gdp", 1, 3, special, 7,7, 1971, 1980, 1981, 1990, 2, 1960, 2003)

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
    #is it just that we're doing a placebo treatment for the treatment being at 1980???
    time.optimize.ssr = 1960:1980,
    unit.names.variable = 2,
    time.plot = 1960:1989
  )

# fit main model with v from training model
synth.out <- synth(
  data.prep.obj=dataprep.out,
  custom.v=as.numeric(vw)
)

# Save the counterfactual results
w_synth_co <- synth.out$solution.w
Y_synth_co <- dataprep.out$Y0 %*% w_synth_co
Y_true_co <- dataprep.out$Y1

## Compute the counterfactual standard errors
#only do this across units and not across time
#my main question is why is this separated out from the stuff we did previously? Couldn't we just save it when we get to 1980 in the previous loops??
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
  
  vw <- find_vweights(d, c("gdp","trade","infrate"), "gdp", 1, 3, special, i,j, 1971, 1980, 1981, 1990, 2, 1960, 2003)
  
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
      #yeah I guess that the counterfactual exercise is just saying "treatment happens in 1980". I'll compare to other methods to see if I still think that
      time.optimize.ssr = 1960:1980,
      unit.names.variable = 2,
      time.plot = 1960:1989
    )
  
  # fit main model with v from training model
  synth.out <- synth(
    data.prep.obj=dataprep.out,
    custom.v=as.numeric(vw)
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


############################################################################


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
