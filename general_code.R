## Clear
rm(list = ls())

library(R.matlab)
library(shape)
library(ggplot2) 

source('functions.R')

# Load data
d <- read.dta("repgermany.dta")

################################
########## Parameters ##########
################################

# Maybe needed, maybe not?

###################################################
########## Estimates and Standard Errors ##########
###################################################



# prep_params[[1]] = pred: vector containing string with predictor variables 
# prep_params[[2]] = dep: string specifying which one is dependent variable
# prep_params[[3]] = u: integer identifying unit variable (which column in data frame specifies index!)
# prep_names[[4] = t: integer identifying time variable (which column in data frame specifies time!)
# prep_names[[5]] = spec: list of special predictors

# prep_names[[6]] = subs: together with j identifies control identifiers -> should we just give it as a single input?
# prep_names[[7]] = years: vector specifying [1]: start predictor priors, [2]: end predictor priors
#                                            [3]: start time optimize ssr, [4]: end time optimize ssr
#                                            [5]: start time plot, [6]: end time plot
# prep_names[[8]] = names: unit names variable (which column in data frame specifies unit names!)

prep_params_list <- list(c("gdp","trade","infrate"),
                         "gdp",1,3,list(
                           list("industry" ,1981:1990, c("mean")),
                           list("schooling",c(1980,1985), c("mean")),
                           list("invest80" ,1980, c("mean"))
                         ),unique(d$index)[-7], c(1981,1990,1960,1989,1960,2003),2)


# Get estimates and standard errors
fit_elastic_net <- general_estimate(d, method = "elastic_net", 
                                    prep_params= prep_params_list, 
                                    tune_params = list(c(seq(from = 1e-02, to = 1e-01, by = 1e-02),
                                                            seq(from = 2e-01, to = 100, by = 1e-01), 
                                                            seq(from = 200, to = 50000, by = 100)), seq(from = 0.1, to = 0.9, by = 0.1)),
                                                      ind_treatment=7)
fit_constr_reg <- general_estimate(d, method = "constr_reg", 
                                    prep_params= prep_params_list, ind_treatment=7)



tune_spec_predict <- list(
  list("industry" ,1971:1980, c("mean")),
  list("schooling",c(1970,1975), c("mean")),
  list("invest70" ,1980, c("mean"))
  )

tune_years <- c(1971, 1980, 1981, 1990, 1960, 2003)


#getting the estimate works. Still need to check standard errors. Is our plan to run the standard error code _inside_ the general_estimate funciton?
fit_synth <- general_estimate(d, method = "synth", 
                                   prep_params= prep_params_list,
                                    tune_params = list(tune_spec_predict, tune_years), 7)


#we can put it inside the other function but for now I'm just going to do this.

source('functions.R')
se_synth <- general_se(d, method="synth", se_method="units_time", prep_params=prep_params_list, tune_params=list(tune_spec_predict, tune_years),  ind_treatment=7)

#################################
########## Save Values ##########
#################################

# Save matrices for future reference and plots

# Elastic Net 
save(list = c("w", "int", "Y_est", "Y_true", 
              "std_err_i", "std_err_t", "std_err_it"), 
     file = "germ_elast_nocov.RData")
writeMat("germ_elast_nocov.mat", 
         w = w_elast, int = int_elast, 
         Y_est = Y_est_elast, Y_true = Y_true, 
         std_err_i = std_err_i_elast, 
         std_err_t = std_err_t_elast, 
         std_err_it = std_err_it_elast)

# Best subset 
save(list = c("w", "int", "Y_est", "Y_true", 
              "std_err_i", "std_err_t", "std_err_it"), 
     file = "germ_subs_nocov.RData")
writeMat("germ_subs_nocov.mat", 
         w = w_subs, int = int_subs, 
         Y_est = Y_est_subs, Y_true = Y_true, 
         std_err_i = std_err_i_subs, 
         std_err_t = std_err_t_subs, 
         std_err_it = std_err_it_subs)

# Constrained regression
save(list = c("w", "int", "Y_est", "Y_true", 
              "std_err_i", "std_err_t", "std_err_it"), 
     file = "germ_constr_reg_nocov.RData")
writeMat("germ_constr_reg_nocov.mat", 
         w = w_constr_reg, int = int_constr_reg, 
         Y_est = Y_est_constr_reg, Y_true = Y_true, 
         std_err_i = std_err_i_constr_reg, 
         std_err_t = std_err_t_constr_reg, 
         std_err_it = std_err_it_constr_reg)

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
