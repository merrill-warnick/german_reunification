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

#first, we load in the German data
d <- read.dta("repgermany.dta")

## Table 1 & 2, Figure 1, 2, & 3

#choose values for the v
#I don't see where we do that yet

# data setup for training model
cat('*** Main ***\n')

#dataprep is a function that prepares your panel data for using synth
dataprep.out <-
  dataprep(
    #name of the dataframe
           foo = d,
           #names of predictor variables
           predictors    = c("gdp","trade","infrate"),
           #name of dependent variable
           dependent     = "gdp",
           #give the value of the unit identifier column
           unit.variable = 1,
           #give the value of the time identifier column
           time.variable = 3,
           #special predictors, I think this is mostly for control variables
           special.predictors = list(
            list("industry", 1971:1980, c("mean")),
            list("schooling",c(1970,1975), c("mean")),
            list("invest70" ,1980, c("mean"))
           ),
           #give the value of the treatment identifier column
           treatment.identifier = 7,
           
           controls.identifier = unique(d$index)[-7],
           #identify prior periods
           
           #ohhhhh so if you look at the time periods on the time.predictors and time.optimize, they're different from the 
           #one that's lower down. 
           
           #yes look closesr and you can see this does the cross-validation
           #I guess one thing to do is to compare this cross-validation method to the one that we use in the general case, maybe they use
           #years in the same way to choose tuning parameters? That could be useful when we generalize the code.
           
           #no, so the general cross-validation exercize is different and is more like the placebo test.
           
           #anyway, what's going on here is that it optimizes and finds v-weights for us. We use the first half of the prior period
           #as the predictors and the second half of the prior period as the out-of sample fitting set.
           
           time.predictors.prior = 1971:1980,
           #which periods we want to optimize over
           time.optimize.ssr = 1981:1990,
           
           unit.names.variable = 2,
           #time period that we want to make plots over for a different command we'll use later
           time.plot = 1960:2003
         )

# fit training model
# need to plug in the structure we got out of the dataprep function into synth

#running synth here lets us grab the vs that we're going to use later
synth.out <- 
  synth(
        data.prep.obj=dataprep.out,
        Margin.ipop=.005,Sigf.ipop=7,Bound.ipop=6
        #when we don't specify vs it optimizes for us, letting us find the vs that we're actually going to use
        )

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
  custom.v=as.numeric(synth.out$solution.v)
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
# Over units
std_err_i <- matrix(0, N - 1, T1)
for (j in 1:(N - 1)) {
  i <- units_co[j]
  cat('(Std. Error) Over Unit i =', toString(i,j), '\n')
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
      #yeah, see, this is where we tell it that the wrong guy is the treated unit.
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
      #okay my one quibbble with this bit is that when we're changing the treatment year, shouldn't it change what years
      #we're using for the cross-validation?
      #yeah, in the paper say say "we divide the pretreatment years into a training period from blah to blah but if you're chanign the treatment 
      #year, I wonder if that kinda messes this up. Idk if it will matter that much but that's something to think about.
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
cat(toString(std_err_it))
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



#okay so this works now and gives same results as the original, prolem in the counterfactual step i guess?

#okay, what is this section doing?
#tbh I'm not totally clear on what the counterfactual exercise is doing in the paper
#or rather I'm not sure how it works, maybe I can read the code and figure it out.

#this counterfactual section seems redunant, I think it's recalculating stuff we did in previous loops--the across units and time bit.

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
    #is it just that we're doing a placebo treatment for the treatment being at 1980???
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
      #yeah I guess that the counterfactual exercise is just saying "treatment happens in 1980". I'll compare to other methods to see if I still think that
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
