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


dat <- readMat("boat_data.mat")
X <- dat$X
Y <- dat$Y
Z <- dat$Z

data <- list('Y' <- dat$Y, 'Z' = dat$Z, 'X' = dat$X)

W <- matrix(0, 1, dim(Y)[2])

W(1) <- 1

#in this setting, we just use the same Y, Z, X to choose v-weights in synth
tune_params_synth <- data


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
writeMat("boat_elast_nocov.mat", 
         w = fit_elastic_net$w, int = fit_elastic_net$int, 
         Y_est = fit_elastic_net$Y_est, Y_true = fit_elastic_net$Y_true, 
         std_err_i = fit_elastic_net$std_err_i, 
         std_err_t = fit_elastic_net$std_err_t, 
         std_err_it = fit_elastic_net$std_err_it)

# Best subset 
#save(list = c("w", "int", "Y_est", "Y_true", 
#              "std_err_i", "std_err_t", "std_err_it"), 
#     file = "germ_subs_nocov.RData")
writeMat("boat_subs_nocov.mat", 
         w = fit_subs$w, int = fit_subs$int, 
         Y_est = fit_subs$Y_est, Y_true = fit_subs$Y_true, 
         std_err_i = fit_subs$std_err_i, 
         std_err_t = fit_subs$std_err_t, 
         std_err_it = fit_subs$std_err_it)

# Constrained regression
#save(list = c("w", "int", "Y_est", "Y_true", 
#              "std_err_i", "std_err_t", "std_err_it"), 
#     file = "germ_constr_reg_nocov.RData")
writeMat("boat_constr_reg_nocov.mat", 
         w = fit_constr_reg$w, int = fit_constr_reg$int, 
         Y_est = fit_constr_reg$Y_est, Y_true = fit_constr_reg$Y_true, 
         std_err_i = fit_constr_reg$std_err_i, 
         std_err_t = fit_constr_reg$std_err_t, 
         std_err_it = fit_constr_reg$std_err_it)

###########################
########## Plots ##########
###########################

## Load data -> or can also simply call from before
data_did <- readMat('boat_did_nocov.mat')
data_elast <- readMat('boat_elast_nocov.mat')
data_subset <- readMat('boat_subs_nocov.mat')
data_synth <- readMat('boat_synth.mat')
data_constr <- readMat('boat_constr_reg_nocov.mat')

### Treatment figure
plot(1973:1991, data_did$Y.true, type = "l", lty = 2, ylim = c(0, 35000), xlim = c(1973,1991), col = "red", main = "Mariel Boatlift: Log Weekly Wages", xlab = "Year", ylab = "", las = 1, bty = 'L')
lines(1973:1991, data_did$Y.did, lty = 1, col= "yellow")
lines(1973:1991, data_elast$Y.elast, lty = 1, col= "purple4")
lines(1973:1991, data_subset$Y.subs, lty = 1, col= "orange")
lines(1973:1991, data_synth$Y.synth, lty = 1, col= "blue")
abline(v = 1989, col="black")
abline(v = 1975, col = "grey96")
abline(v = 1980, col = "grey96")
abline(v = 1985, col = "grey96")
abline(v = 1990, col = "grey96")
abline(v = 1995, col = "grey96")
legend("topleft",legend=c("Actual data","Difference-in-Differences", expression(paste("Elastic net (opt. ", lambda," and ",alpha,")" )),"Best subset (opt. k)", "Original synth."), col=c("red","yellow","purple4","orange","blue"),lty=c(2,1,1,1,1), ncol=1, bty = 'n', cex = 0.7)
arrows(x0=1987, y0=32500,x1=1988, y1=32499, col=c("black"), lwd=1 , length = 0.05,xpd=TRUE)
text(x=1981,y=32500,pos=4,label = "Reunification", cex = 0.5)

### Standard Errors

tau <- cbind(data_did$Y.true[8:19]-data_synth$Y.synth[8:19],data_did$Y.true[8:19]-data_elast$Y.elast[8:19]) # cbind for each method
std_err <- cbind(data_synth$std.err.synth.i, data_elast$std.err.elast.i)

plot(1980,1991, tau[,1], type = "l", lty = 1, ylim = c(-100, 100), xlim = c(1980,1991), col = "blue", main = "Mariel Boatlift: Standard Errors", xlab = "Year", ylab = "", las = 1, bty = "L")
lines(1980,1991, tau[,1]+1.96*std_err[,1], lty = 3, col= "blue")
lines(1980,1991, tau[,1]-1.96*std_err[,1], lty = 3, col= "blue")
lines(1980,1991, tau[,2], lty = 1, col= "plum2")
lines(1980,1991, tau[,2]+1.96*std_err[,2], lty = 2, col= "plum2")
lines(1980,1991, tau[,2]-1.96*std_err[,2], lty = 2, col= "plum2")
abline(h = 0, col= "black")
abline(v = 1981, col = "grey96")
abline(v = 1982, col = "grey96")
abline(v = 1983, col = "grey96")
abline(v = 1984, col = "grey96")
abline(v = 1985, col = "grey96")
abline(v = 1986, col = "grey96")
abline(v = 1987, col = "grey96")
abline(v = 1988, col = "grey96")
abline(v = 1989, col = "grey96")
abline(v = 1990, col = "grey96")
abline(v = 1991, col = "grey96")
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



control_names <- c('NYC', 'LA', 'CHG', 'PHL', 'DET', 'SF', 'DC', 'BST', 'SUFF', 'PIT', 'SLO', 'BLT', 'CLV', 'HOU', 'NWK', 'MNN', 'DAL', 'SEA', 'ANH', 'MLWK', 'ATL', 'CIN', 'PAT', 'SD', 'BUF', 'KNS', 'DEN', 'SBRN', 'IND', 'SJC', 'NORL', 'TMP', 'POR', 'COL', 'ROCH', 'SAC', 'FWOR', 'BIR', 'ALB', 'NOR', 'AKR', 'ECHG', 'GRSB')
weights_synth <- as.data.frame(cbind(control_names,weights[,3]))
colnames(weights_synth) <- c("controls","w")
weights_synth$w <- as.numeric(as.character(weights_synth$w))
p <- ggplot(weights_synth, aes(x=controls, y=w))+geom_bar(stat="identity", fill = "blue",color ="black", show.legend = FALSE)+labs(title="",x="", y = "Original synth.")+scale_fill_manual(values = c("royalblue"))+ylim(-1, 1)+coord_flip()
p