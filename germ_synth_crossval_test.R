#The plan in this file is to do some work to prep cross-validation experimentation for vweights in the
#synthetic control part of the project. My plan is to first write a function that finds the vweights generally?
#now I guess the question is, do I want to try to make it whole hog general now? Yeah I guess so, that makes sense, esp since
#during variance calculations I need to assign different treatment units/years.


#this is all the code I need for it, so I just need to go through and make the argument general


#there is something wrong with my code....so I need to work a bit harder on it I guess. I'm going to rerun the OG code once to make sure
#that it gives the same stuff every time, but I guess that what I can do is redo this program, start from the ground up.
#maybe I didn't notice small differences hanging around or something like that.

#now I need to compartmentalize the rest of this

rm(list=ls())
library(foreign)
library(Synth)
library(xtable)
library(R.matlab)


#okay wehn I'm doing standard error over units I'm getting some kind of missing data error for control unit 7. Then when we get to unit 8, it says that the 
#treated unit is among controls.

#idk if the default is gonna work

#I think that we're going to need a bunch of different years parameters so in the root function you just input a years vector and then 
#unpack within the if statements.

find_vweights <- function(d, pred, y, u, t, spec, i,j,cont_set, predyear0, predyear1, optyear0, optyear1, names, year0, year1){
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
    
    #I'm not sure what this does but I think it's just supposed to be -i
    #oh I see what the problem is. We're taking a subset that is not the same as unique(d$index)...okay I need to fix that then
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
    
    #these are just optimization things, for now I'm just going to leave them as is, but I guess if I look at later stuff and see that we make different
    # choices, I'll let these be inputs.
    Margin.ipop=.005,Sigf.ipop=7,Bound.ipop=6
  )
output <- synth.out$solution.v
}



find_ysynth <- function(d, pred, y, u, t, spec, i,j,cont_set, predyear0, predyear1, optyear0, optyear1, names, year0, year1, vweight){

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


synth.out <- synth(
  data.prep.obj=dataprep.out,
  custom.v=as.numeric(vweight)
)

w <- synth.out$solution.w

#okay so this doesn't quite concatenate like I want it to. It just lines things up onto on eline, but
#I actually want it to have separate pieces. Do I put them in a list or something?
#maybe I can use the year variables to make it work out? Easiest would be if I could store them like I want to.
#I can try putting them in a list.
output <- list(y1 = dataprep.out$Y1, ysynth = dataprep.out$Y0 %*% w, w = w)

}

#I think we just feed in ysynth and y1 and tell it the treatment year
#won't work with multiple treated units, I guess we could figure that out later. I don't know what you're supposed to do with multiple
#treated units anyway.
find_treatment <- function(Y1, Ysynth, tyear){
  
  output <-Y1[tyear] - Ysynth[tyear]
}


#we could either tell it the number of units or just feed it in I'll just feed in for now.
se_unit <- function(N, T, T0, d, pred, y, u, t, cspec, spec, cont_set, cyears, years, names){
  
  T1<- T - T0
  
std_err_i <- matrix(0, N - 1, T1)
for (j in 1:(N - 1)) {
  i <- cont_set[j]
  cat('(Std. Error) Over Unit i =', toString(i), '\n')
  
  vw <- find_vweights(d, pred, y, u, t, cspec, i,j,cont_set, cyears[1], cyears[2], cyears[3], cyears[4], names, cyears[5], cyears[6])
  
  y_both <- find_ysynth(d, pred, y, u, t, spec, i,j, cont_set, years[1], years[2], years[3], years[4], names, years[5], years[6], vw)
    
  Y1 <- y_both$y1
  Y_synth<-y_both$ysynth
  
  #I don't quite understand the indexing but I see what is going on here, I think.
  std_err_i[j,] <- (Y1[-c(1:T0),] - Y_synth[-c(1:T0),]) ^ 2
}

std_err_i <- as.matrix(sqrt(apply(std_err_i, 2, mean)))
}

se_time <- function(N, T, T0, d, pred, y, u, t, cspec,i,j, spec, cont_set, cyears, years, names){
  
  s <- floor(T0 / 2)
  std_err_t <- matrix(0, s, 1)
  for (k in 1:s) {
    cat('(Std. Error) Over Time t =', toString(k), '\n')
    
    vw <- find_vweights(d, pred, y, u, t, cspec, i,j,cont_set, cyears[1], cyears[2], cyears[3], cyears[4], names, cyears[5], cyears[6])
    
    y_both <- find_ysynth(d, pred, y, u, t, spec, i,j, cont_set, years[1], years[2], years[3], years[4]-k, names, years[5], years[6], vw)
    
    Y1 <- y_both$y1
    Y_synth<-y_both$ysynth
    
    
    std_err_t[k,1] <- (Y1[T0 - k + 1,] - Y_synth[T0 - k + 1,]) ^ 2
  }
  std_err_t <- as.matrix(sqrt(apply(std_err_t, 2, mean)))
}


#hmmm should saving vweights be an option or something like that? We only need to do it for now....
#i'll just keep it now then delete it later I guess.

#something is going wrong in this one, we got the wrong values and the wrong vweights I think, meaning the vweights didn't even come out right


se_unit_time <- function(N, T, T0, d, pred, y, u, t, cspec, spec, cont_set, cyears, years, names){
 
  s <- floor(T0 / 2)
  #this is a magic number rn
  vweights_it<-array(0, dim=c(N-1,s,6))
  
  std_err_it <- matrix(0, N - 1, 1)
  for (j in 1:(N - 1)) {
    i <- cont_set[j]
    std_err_temp <- matrix(0, s, 1)
    for (k in 1:s) {
      cat('(Std. Error) Over Unit and Time ( i , t ) = (',toString(i), ',', 
          toString(k), ')\n')
      
      vw <- find_vweights(d, pred, y, u, t, cspec, i,j,cont_set, cyears[1], cyears[2], cyears[3], cyears[4], names, cyears[5], cyears[6])
      
      #save vweights so we can look at them
      vweights_it[j,t,] <- as.matrix(vw)
      
      y_both <- find_ysynth(d, pred, y, u, t, spec, i,j, cont_set, years[1], years[2], years[3], years[4] - k, names, years[5], years[6], vw)
      
     
      
      # Solutions
      Y1 <- y_both$y1
      Y_synth<-y_both$ysynth
      
      std_err_temp[t,1] <- (Y1[T0 - k + 1,] - Y_synth[T0 - k + 1,]) ^ 2
    }
    std_err_temp <- as.matrix(apply(std_err_temp, 2, mean))
    std_err_it[j,1] <- std_err_temp
  }
  cat(toString(std_err_it))
  std_err_it <- list(se = as.matrix(sqrt(apply(std_err_it, 2, mean))), vw = vweights_it)
}

#std_err_it is wrong I think

## Replication Code for
# A. Abadie, A. Diamond, and J. Hainmueller. 2014.
# Comparative Politics and the Synthetic Control Method
# American Journal of Political Science.


#first, we load in the German data
d <- read.dta("repgermany.dta")

## Table 1 & 2, Figure 1, 2, & 3

cat('*** Main ***\n')


#I think that special is correct now, I think that we just plug it in to all the find_vweights.
#it is different for main-mmodel stuff.
cspecial<-list(
  list("industry" ,1971:1980, c("mean")),
  list("schooling",c(1970,1975), c("mean")),
  list("invest70" ,1980, c("mean"))
)

special <- list(
  list("industry" ,1981:1990, c("mean")),
  list("schooling",c(1980,1985), c("mean")),
  list("invest80" ,1980, c("mean"))
)

cyears <-c(1971, 1980, 1981, 1990, 1960, 2003)

years <- c(1981, 1990, 1960, 1989, 1960, 2003)

predict<- c("gdp","trade","infrate")


#vweights
vw <- find_vweights(d, predict, "gdp", 1, 3, cspecial, 7,7,unique(d$index), cyears[1],cyears[2],cyears[3], cyears[4], 2, cyears[5], cyears[6])

y_out <- find_ysynth(d, predict, "gdp", 1, 3, special, 7,7,unique(d$index), years[1],years[2],years[3], years[4], 2, years[5], years[6], vw)


w_synth <- y_out$w
Y_synth <- y_out$ysynth
Y1 <- y_out$y1

teffect <- find_treatment(Y1, Y_synth, 6)


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
std_err_i <- se_unit(N, T, T0, d, predict, "gdp", 1, 3, cspecial, special, units_co, cyears, years, 2)

std_err_t <- se_time(N, T, T0, d, predict, "gdp", 1, 3, cspecial, 7, 7,  special, units_co, cyears, years, 2)

se_ut_output <- se_unit_time(N, T, T0, d, predict, "gdp", 1, 3, cspecial, special, units_co, cyears, years, 2)

#yeah after printing them, these are still much too small. i'll need to run small examples to test and figure out what the problem is.
#It hinki'll do like three units with like the same nubmer of years each or something.
std_err_it <- se_ut_output$se
vweights <- se_ut_output$vw


#I'll rewrite the counterfactual part later after this runs.
#I think you just change the years or somethign




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

vw <- find_vweights(d, c("gdp","trade","infrate"), "gdp", 1, 3, special, 7,7,unique(d$index), 1971, 1980, 1981, 1990, 2, 1960, 2003)

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
  
  vw <- find_vweights(d, c("gdp","trade","infrate"), "gdp", 1, 3, special, i,j,units_co, 1971, 1980, 1981, 1990, 2, 1960, 2003)
  
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
