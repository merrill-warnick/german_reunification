## Clear
rm(list = ls())

library(R.matlab)
library(shape)
library(ggplot2)
library(egg)

source('functions.R')

################################
########## Parameters ##########
################################

#we could add a parameters section later to make things easier, but it depends on what we want to do.


################################
########## Data Prep ###########
################################


data <- readMat("boat_data.mat")
X <- data$X
Y <- data$Y
Z <- data$Z


W <- matrix(0, 1, dim(Y)[2])

W[1] <- 1

#in this setting, we just use the same Y, Z, X to choose v-weights in synth
tune_params_synth <- data


###################################################
########## Estimates and Standard Errors ##########
###################################################


fit_elastic_net <- general_estimate(data$Y, data$Z, data$X, W, method = "elastic_net", 
                                    tune_params = list(c(seq(from = 1e-04, to = 1e-03, by = 1e-04),
                                                         seq(from = 2e-03, to = 1e-02, by = 1e-03),
                                                         seq(from = 2e-02, to = 1e-01, by = 1e-02),
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
save(fit_elastic_net, file = "boat_elast_nocov.RData")
writeMat("boat_elast_nocov.mat", 
         w = fit_elastic_net$w, int = fit_elastic_net$int, 
         Y_est = fit_elastic_net$Y_est, Y_true = fit_elastic_net$Y_true, 
         std_err_i = fit_elastic_net$std_err_i, 
         std_err_t = fit_elastic_net$std_err_t, 
         std_err_it = fit_elastic_net$std_err_it)

# Constrained regression
save(fit_constr_reg, file = "boat_constr_reg_nocov.RData")
writeMat("boat_constr_reg_nocov.mat", 
         w = fit_constr_reg$w, int = fit_constr_reg$int, 
         Y_est = fit_constr_reg$Y_est, Y_true = fit_constr_reg$Y_true, 
         std_err_i = fit_constr_reg$std_err_i, 
         std_err_t = fit_constr_reg$std_err_t, 
         std_err_it = fit_constr_reg$std_err_it)

# Best subset 
save(fit_subs, file = "boat_subs_nocov.RData")
writeMat("boat_subs_nocov.mat", 
         w = fit_subs$w, int = fit_subs$int, 
         Y_est = fit_subs$Y_est, Y_true = fit_subs$Y_true, 
         std_err_i = fit_subs$std_err_i, 
         std_err_t = fit_subs$std_err_t, 
         std_err_it = fit_subs$std_err_it)

# Synthetic Control
save(fit_synth, file = "boat_synth_nocov.RData")
writeMat("boat_synth_nocov.mat", 
         w = fit_synth$w, int = fit_synth$int, 
         Y_est = fit_synth$Y_est, Y_true = fit_synth$Y_true, 
         std_err_i = fit_synth$std_err_i, 
         std_err_t = fit_synth$std_err_t, 
         std_err_it = fit_synth$std_err_it)

# Diff-in-diff
save(fit_diff_in_diff, file = "boat_did_nocov.RData")
writeMat("boat_did_nocov.mat", 
         w = fit_diff_in_diff$w, int = fit_diff_in_diff$int, 
         Y_est = fit_diff_in_diff$Y_est, Y_true = fit_diff_in_diff$Y_true, 
         std_err_i = fit_diff_in_diff$std_err_i, 
         std_err_t = fit_diff_in_diff$std_err_t, 
         std_err_it = fit_diff_in_diff$std_err_it)



####################################
########## Counterfactual ##########
####################################
T0 <- T0 <- dim(Z)[1]
T0_co <- 4
T1_co <- T0 - T0_co

Y_co = data$Y[1:T0,]
Z_co = data$Z[1:T0_co,]
X_co = data$X

fit_elastic_net_co <- general_estimate(Y_co, Z_co, X_co, W, method = "elastic_net", 
                                       tune_params = list(c(seq(from = 1e-02, to = 1e-01, by = 1e-02),
                                                            seq(from = 2e-01, to = 100, by = 1e-01), 
                                                            seq(from = 200, to = 50000, by = 100)), seq(from = 0.1, to = 0.9, by = 0.1)))
fit_synth_co <- general_estimate(Y_co, Z_co, X_co, W, method = "synth", tune_params = tune_params_synth)

save(fit_elastic_net_co, file = "boat_en_co_nocov.RData")
save(fit_synth_co, file = "boat_synth_co_nocov.RData")

###########################
########## Plots ##########
###########################

### Treatment figure
plot(1973:1991, fit_diff_in_diff$Y_true, type = "l", lty = 2, ylim = c(4.5, 5.5), xlim = c(1973,1991), col = "red", main = "Mariel Boatlift: Log Weekly Wages", xlab = "Year", ylab = "", las = 1, bty = 'L')
lines(1973:1991, fit_synth$Y_est, lty = 1, col= "blue")
lines(1973:1991, fit_constr_reg$Y_est, lty = 1, col= "green")
lines(1973:1991, fit_elastic_net$Y_est, lty = 1, col= "purple4")
lines(1973:1991, fit_subs$Y_est, lty = 1, col= "orange")
lines(1973:1991, fit_diff_in_diff$Y_est, lty = 1, col= "yellow")
abline(v = 1979, col="black")
abline(v = 1975, col = "grey96")
abline(v = 1980, col = "grey96")
abline(v = 1985, col = "grey96")
abline(v = 1990, col = "grey96")
abline(v = 1995, col = "grey96")
legend("topright",legend=c("Actual data", "Original synth.", "Regression w/restrictions",expression(paste("Elastic net (opt. ", lambda," and ",alpha,")" )),"Best subset (opt. k)","Difference-in-Differences"), col=c("red","blue","green","purple4","orange","yellow"),lty=c(2,1,1,1,1,1), ncol=1, bty = 'n', cex = 0.7)
arrows(x0=1978, y0=5.45,x1=1978.5, y1=5.45, col=c("black"), lwd=1 , length = 0.05,xpd=TRUE)
text(x=1976.3,y=5.45,pos=4,label = "Policy", cex = 0.65)


### Standard Errors
tau <- cbind(fit_synth$Y_true[8:19]-fit_synth$Y_est[8:19],fit_elastic_net$Y_true[8:19]-fit_elastic_net$Y_est[8:19]) # cbind for each method
std_err <- cbind(fit_synth$std_err_i, fit_elastic_net$std_err_i)

plot(1980:1991, tau[,1], type = "l", lty = 1, ylim = c(-1, 1), xlim = c(1980,1991), col = "blue", main = "Mariel Boatlift: Standard Errors", xlab = "Year", ylab = "", las = 1, bty = "L")
lines(1980:1991, tau[,1]+1.96*std_err[,1], lty = 3, col= "blue")
lines(1980:1991, tau[,1]-1.96*std_err[,1], lty = 3, col= "blue")
lines(1980:1991, tau[,2], lty = 1, col= "darkmagenta")
lines(1980:1991, tau[,2]+1.96*std_err[,2], lty = 2, col= "darkmagenta")
lines(1980:1991, tau[,2]-1.96*std_err[,2], lty = 2, col= "darkmagenta")
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
legend("topright",legend=c("ADH synth. treatment","ADH treatment +/-1.96*std.err.",expression(paste("Elastic net treatment (opt. ", lambda," and ",alpha,")" )),"Elastic net treatment +/-1.96*std.err."), col=c("blue","blue","darkmagenta","darkmagenta"),lty=c(1,3,1,2), ncol=1, bty = 'n', cex = 0.65)


### Standard Errors Counterfactual
tau <- cbind(fit_synth$Y_true[5:7]-fit_synth$Y_est[5:7],fit_elastic_net$Y_true[5:7]-fit_elastic_net$Y_est[5:7]) # cbind for each method
std_err <- cbind(fit_synth$std_err_i, fit_elastic_net$std_err_i)

plot(1977:1979, tau[,1], type = "l", lty = 1, ylim = c(-1, 1), xlim = c(1977,1979), col = "blue", main = "Mariel Boatlift: Counterfactual Treatment", xlab = "Year", ylab = "", las = 1, bty = "L")
lines(1977:1979, tau[,1]+1.96*std_err[,1], lty = 3, col= "blue")
lines(1977:1979, tau[,1]-1.96*std_err[,1], lty = 3, col= "blue")
lines(1977:1979, tau[,2], lty = 1, col= "darkmagenta")
lines(1977:1979, tau[,2]+1.96*std_err[,2], lty = 2, col= "darkmagenta")
lines(1977:1979, tau[,2]-1.96*std_err[,2], lty = 2, col= "darkmagenta")
abline(h = 0, col= "black")
abline(v = 1978, col = "grey96")
abline(v = 1979, col = "grey96")
legend("topright",legend=c("ADH synth. treatment","ADH treatment +/-1.96*std.err.",expression(paste("Elastic net treatment (opt. ", lambda," and ",alpha,")" )),"Elastic net treatment +/-1.96*std.err."), col=c("blue","blue","darkmagenta","darkmagenta"),lty=c(1,3,1,2), ncol=1, bty = 'n', cex = 0.65)


## Weights
weights <- cbind(fit_synth$w, fit_constr_reg$w, fit_elastic_net$w, fit_subs$w)
theme_set(theme_bw())

control_names <- c('NYC', 'LA', 'CHG', 'PHL', 'DET', 'SF', 'DC', 'BST', 'SUFF', 'PIT', 'SLO', 'BLT', 'CLV', 'HOU', 'NWK', 'MNN', 'DAL', 'SEA', 'ANH', 'MLWK', 'ATL', 'CIN', 'PAT', 'SD', 'BUF', 'KNS', 'DEN', 'SBRN', 'IND', 'SJC', 'NORL', 'TMP', 'POR', 'COL', 'ROCH', 'SAC', 'FWOR', 'BIR', 'ALB', 'NOR', 'AKR', 'ECHG', 'GRSB')
weights_synth <- as.data.frame(cbind(control_names,weights[,1]))
colnames(weights_synth) <- c("controls","w")
weights_synth$w <- as.numeric(as.character(weights_synth$w))
p <- ggplot(weights_synth, aes(x=controls, y=w))+geom_bar(stat="identity", fill = "blue",color ="black", show.legend = FALSE)+labs(title="",x="", y = "Original synth.")+scale_fill_manual(values = c("royalblue"))+ylim(-1, 1)+coord_flip()

weights_constr_reg <- as.data.frame(cbind(control_names,weights[,2]))
colnames(weights_constr_reg) <- c("controls","w")
weights_constr_reg$w <- as.numeric(as.character(weights_constr_reg$w))
p1 <- ggplot(weights_constr_reg, aes(x=controls, y=w))+geom_bar(stat="identity", fill = "forestgreen",color ="black", show.legend = FALSE)+labs(title="",x="", y = "Reg./w.restr.")+scale_fill_manual(values = c("forestgreen"))+ylim(-1, 1)+coord_flip()

weights_elastic_net <- as.data.frame(cbind(control_names,weights[,3]))
colnames(weights_elastic_net) <- c("controls","w")
weights_elastic_net$w <- as.numeric(as.character(weights_elastic_net$w))
p2 <- ggplot(weights_elastic_net, aes(x=controls, y=w))+geom_bar(stat="identity", fill = "purple",color ="black", show.legend = FALSE)+labs(title="",x="", y = "Elastic Net")+scale_fill_manual(values = c("purple"))+ylim(-1, 1)+coord_flip()

weights_subs <- as.data.frame(cbind(control_names,weights[,4]))
colnames(weights_subs) <- c("controls","w")
weights_subs$w <- as.numeric(as.character(weights_subs$w))
p3 <- ggplot(weights_subs, aes(x=controls, y=w))+geom_bar(stat="identity", fill = "orange",color ="black", show.legend = FALSE)+labs(title="",x="", y = "Best subset")+scale_fill_manual(values = c("orange"))+ylim(-1, 1)+coord_flip()

figure <- ggarrange(p, p1, p2,p3,
                    ncol = 4, nrow = 1, top = "Mariel Boatlift: Weights")
figure