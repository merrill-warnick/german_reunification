## Clear
rm(list = ls())

library(R.matlab)
library(shape)
library(ggplot2)

source('functions.R')


################################
########## Parameters ##########
################################

#we could add a parameters section later to make things easier, but it depends on what we want to do.


################################
########## Data Prep ###########
################################

# Load data
d <- read.dta("repgermany.dta")


#first, prepare the data to be used for choosing v-weights in synth

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

#Extract the treatment and control units from the dataprep object  
X0 <- dataprep.out$X0
X1 <- dataprep.out$X1

Z1 <- dataprep.out$Z1
Z0 <- dataprep.out$Z0

Y1 <- dataprep.out$Y1plot
Y0 <- dataprep.out$Y0plot

#Re-bind data so that the treated unit is in the first position
# [X1,X0]
X <- cbind(X1, X0)

# [Y1,Y0]
Y <- cbind(Y1, Y0)

# [Z1,Z0]
Z <- cbind(Z1, Z0)

tune_params_synth <- list( "X" = X, "Y" = Y, "Z" = Z)



#prepare the actual estimation data

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
  
#Extract the treatment and control units from the dataprep object  
X0 <- dataprep.out$X0
X1 <- dataprep.out$X1
  
Z1 <- dataprep.out$Z1
Z0 <- dataprep.out$Z0
  
Y1 <- dataprep.out$Y1plot
Y0 <- dataprep.out$Y0plot
  
#Re-bind data so that the treated unit is in the first position
# [X1,X0]
X <- cbind(X1, X0)
  
# [Y1,Y0]
Y <- cbind(Y1, Y0)
  
# [Z1,Z0]
Z <- cbind(Z1, Z0)
  
output <- list( "X" = X, "Y" = Y, "Z" = Z)


#"spec" and "years" are the parameters that are different here.
tune_params_synth <- prep_data(d = d, 
                               ind_treat = 7,
                               pred = c("gdp","trade","infrate"), 
                               dep = "gdp", 
                               u = 1,
                               t = 3, 
                               spec = list(
                                 list("industry", 1971:1980, c("mean")),
                                 list("schooling",c(1970,1975), c("mean")),
                                 list("invest70" ,1980, c("mean"))
                               ),
                               cont_set = unique(d$index)[-7],
                               years = c(1971, 1980, 1981, 1990, 1960, 2003),
                               names = 2)

#prepare W:

W <- matrix(0, 1, dim(data$Y)[2])

W[7] <- 1


###################################################
########## Estimates and Standard Errors ##########
###################################################


fit_elastic_net <- general_estimate(data$Y, data$Z, data$X, W, method = "elastic_net", 
                                    tune_params = list(c(seq(from = 1e-02, to = 1e-01, by = 1e-02),
                                                         seq(from = 2e-01, to = 100, by = 1e-01), 
                                                         seq(from = 200, to = 50000, by = 100)), seq(from = 0.1, to = 0.9, by = 0.1)))

fit_constr_reg <- general_estimate(data$Y, data$Z, data$X, W, method = "constr_reg")

fit_subs <- general_estimate(data$Y, data$Z, data$X, W, method = "best_subset")

fit_synth <- general_estimate(data$Y, data$Z, data$X, W, method = "synth", tune_params = tune_params_synth)

fit_diff_in_diff <- general_estimate(data$Y, data$Z, data$X, W, method = "diff_in_diff")



#################################
########## Save Values ##########
#################################

# Save matrices for future reference and plots

# Elastic Net 
#save(list = c("w", "int", "Y_est", "Y_true", 
#              "std_err_i", "std_err_t", "std_err_it"), 
#     file = "germ_elast_nocov.RData")
writeMat("germ_elast_nocov.mat", 
         w = fit_elastic_net$w, int = fit_elastic_net$int, 
         Y_est = fit_elastic_net$Y_est, Y_true = fit_elastic_net$Y_true, 
         std_err_i = fit_elastic_net$std_err_i, 
         std_err_t = fit_elastic_net$std_err_t, 
         std_err_it = fit_elastic_net$std_err_it)

# Best subset 
#save(list = c("w", "int", "Y_est", "Y_true", 
#              "std_err_i", "std_err_t", "std_err_it"), 
#     file = "germ_subs_nocov.RData")
writeMat("germ_subs_nocov.mat", 
         w = fit_subs$w, int = fit_subs$int, 
         Y_est = fit_subs$Y_est, Y_true = fit_subs$Y_true, 
         std_err_i = fit_subs$std_err_i, 
         std_err_t = fit_subs$std_err_t, 
         std_err_it = fit_subs$std_err_it)

# Constrained regression
#save(list = c("w", "int", "Y_est", "Y_true", 
#              "std_err_i", "std_err_t", "std_err_it"), 
#     file = "germ_constr_reg_nocov.RData")
writeMat("germ_constr_reg_nocov.mat", 
         w = fit_constr_reg$w, int = fit_constr_reg$int, 
         Y_est = fit_constr_reg$Y_est, Y_true = fit_constr_reg$Y_true, 
         std_err_i = fit_constr_reg$std_err_i, 
         std_err_t = fit_constr_reg$std_err_t, 
         std_err_it = fit_constr_reg$std_err_it)

###########################
########## Plots ##########
###########################

## Load data -> or can also simply call from before
data_did <- readMat('germ_did_nocov.mat')
data_elast <- readMat('germ_elast_nocov.mat')
data_subset <- readMat('germ_subs_nocov.mat')
data_synth <- readMat('germ_synth.mat')
data_constr <- readMat('germ_constr_reg_nocov.mat')

### Treatment figure
plot(1960:2003, data_did$Y.true, type = "l", lty = 2, ylim = c(0, 35000), xlim = c(1960,2003), col = "red", main = "West Germany: per capita GDP", xlab = "Year", ylab = "", las = 1, bty = 'L')
lines(1960:2003, data_did$Y.did, lty = 1, col= "yellow")
lines(1960:2003, data_elast$Y.elast, lty = 1, col= "purple4")
lines(1960:2003, data_subset$Y.subs, lty = 1, col= "orange")
lines(1960:2003, data_synth$Y.synth, lty = 1, col= "blue")
abline(v = 1989, col="black")
abline(v = 1960, col = "grey96")
abline(v = 1970, col = "grey96")
abline(v = 1980, col = "grey96")
abline(v = 1990, col = "grey96")
abline(v = 2000, col = "grey96")
legend("topleft",legend=c("Actual data","Difference-in-Differences", expression(paste("Elastic net (opt. ", lambda," and ",alpha,")" )),"Best subset (opt. k)", "Original synth."), col=c("red","yellow","purple4","orange","blue"),lty=c(2,1,1,1,1), ncol=1, bty = 'n', cex = 0.7)
arrows(x0=1987, y0=32500,x1=1988, y1=32499, col=c("black"), lwd=1 , length = 0.05,xpd=TRUE)
text(x=1981,y=32500,pos=4,label = "Reunification", cex = 0.5)

### Standard Errors
tau <- cbind(data_did$Y.true[31:44]-data_synth$Y.synth[31:44],data_did$Y.true[31:44]-data_elast$Y.elast[31:44]) # cbind for each method
std_err <- cbind(data_synth$std.err.synth.i, data_elast$std.err.elast.i)

plot(1990:2003, tau[,1], type = "l", lty = 1, ylim = c(-12500, 12500), xlim = c(1990,2003), col = "blue", main = "West Germany: Standard Errors", xlab = "Year", ylab = "", las = 1, bty = "L")
lines(1990:2003, tau[,1]+1.96*std_err[,1], lty = 3, col= "blue")
lines(1990:2003, tau[,1]-1.96*std_err[,1], lty = 3, col= "blue")
lines(1990:2003, tau[,2], lty = 1, col= "plum2")
lines(1990:2003, tau[,2]+1.96*std_err[,2], lty = 2, col= "plum2")
lines(1990:2003, tau[,2]-1.96*std_err[,2], lty = 2, col= "plum2")
abline(h = 0, col= "black")
abline(v = 1990, col = "grey96")
abline(v = 1992, col = "grey96")
abline(v = 1994, col = "grey96")
abline(v = 1996, col = "grey96")
abline(v = 1998, col = "grey96")
abline(v = 2000, col = "grey96")
abline(v = 2002, col = "grey96")
legend("topright",legend=c("ADH synth. treatment","ADH treatment +/-1.96*std.err.",expression(paste("Elastic net treatment (opt. ", lambda,"and ",alpha,")" )),"Elastic net treatment +/-1.96*std.err."), col=c("blue","blue","plum2","plum2"),lty=c(1,2,1,2), ncol=1, bty = 'n', cex = 0.65)

## Counterfactual
#tau_co <- cbind(data_elast$Y.true.co[(data_elast$T0.co+1):(data_elast$T0.co+data_elast$T1.co)]-data_synth$Y.synth.co[(data_elast$T0.co+1):(data_elast$T0.co+data_elast$T1.co)],data_synth$Y.true.co[(data_elast$T0.co+1):(data_elast$T0.co+data_elast$T1.co)]-data_elast$Y.elast.co[(data_elast$T0.co+1):(data_elast$T0.co+data_elast$T1.co)]) # cbind for each method
#std_err_co <- cbind(data_synth$std.err.synth.i.co, data_elast$std.err.elast.i.co)

#plot(1981:1989, tau_co[,1], type = "l", lty = 1, ylim = c(-12500, 12500), xlim = c(1981,1989), col = "blue", main = "West Germany: Counterfactual", xlab = "Year", ylab = "", las = 1, bty = "L")
#lines(1981:1989, tau_co[,1]+1.96*std_err_co[,1], lty = 3, col= "blue")
#lines(1981:1989, tau_co[,1]-1.96*std_err_co[,1], lty = 3, col= "blue")
#lines(1981:1989, tau_co[,2], lty = 1, col= "plum2")
#lines(1981:1989, tau_co[,2]+1.96*std_err_co[,2], lty = 2, col= "plum2")
#lines(1981:1989, tau_co[,2]-1.96*std_err_co[,2], lty = 2, col= "plum2")
#abline(h = 0, col= "black")
#abline(v = 1982, col = "grey96")
#abline(v = 1983, col = "grey96")
#abline(v = 1984, col = "grey96")
#abline(v = 1985, col = "grey96")
#abline(v = 1986, col = "grey96")
#abline(v = 1987, col = "grey96")
#abline(v = 1988, col = "grey96")
#abline(v = 1989, col = "grey96")
#legend("topright",legend=c("ADH synth. treatment","ADH treatment +/-1.96*std.err.",expression(paste("Elastic net treatment (opt. ", lambda,"and ",alpha,")" )),"Elastic net treatment +/-1.96*std.err."), col=c("blue","blue","plum2","plum2"),lty=c(1,2,1,2), ncol=1, bty = 'n', cex = 0.65)

## Weights
weights <- cbind(data_synth$w.synth, data_elast$w.elast, data_subset$w.subs)
theme_set(theme_bw())

control_names <- c("USA", "GBR", "AUT", "BEL", "DNK", "FRA", "ITA", "NLD", "NOR", "CHE", "JPN", "GRC", "PRT", "ESP","AUS","NZL")
weights_synth <- as.data.frame(cbind(control_names,weights[,3]))
colnames(weights_synth) <- c("controls","w")
weights_synth$w <- as.numeric(as.character(weights_synth$w))
p <- ggplot(weights_synth, aes(x=controls, y=w))+geom_bar(stat="identity", fill = "blue",color ="black", show.legend = FALSE)+labs(title="",x="", y = "Original synth.")+scale_fill_manual(values = c("royalblue"))+ylim(-1, 1)+coord_flip()
p