## Clear
rm(list = ls())

library(R.matlab)
library(shape)
library(ggplot2)
library(ggpubr)
library(egg)

source('functions.R')

################################
########## Parameters ##########
################################

#we could add a parameters section later to make things easier, but it depends on what we want to do.


################################
########## Data Prep ###########
################################


data <- read.table("MLAB_data.txt")
## Built Indices
# California is state no 3, stored in the last column no 39
index_tr <- c(39)
# 38 Control states are 1,2 & 4,5,...,38, stored in columns 1 to 38
index_co <- c(1:38)

# Predcitors are stored in rows 2 to 8
index_predict <- c(2:8)
# Outcome Data is stored in rows 9 to 39; for 1970, 1971,...,2000
index_Y <- c(9:39)

# The pre-treatment time periods are 1970--1980
index_pre <- c(1:19)

## Define Matrices for Predictors
# X0 : 7 X 38 matrix (7 smoking predictors for 38 control states)
X0 <- as.matrix(data[index_predict,index_co])

# X1 : 10 X 1 matrix (10 crime predictors for 1 treated states)
X1 <- as.matrix(data[index_predict,index_tr])

# [X1,X0]
X <- cbind(X1, X0)

## Define Matrices for Outcome Data
# Y0 : 31 X 38 matrix (31 years of smoking data for 38 control states)
Y0 <- as.matrix(data[index_Y,index_co])
# Y1 : 31 X 1 matrix (31 years of smoking data for 1 treated state)
Y1 <- as.matrix(data[index_Y,index_tr])

# [Y1,Y0]
Y <- cbind(Y1, Y0)

# Now pick Z matrices, i.e. the pretreatment period
# over which the loss function should be minmized
# Here we pick Z to go from 1970 to 1988 
# Z0 : 19 X 38 matrix (19 years of pre-treatment smoking data for 38 control states)
Z0 <- as.matrix(Y0[index_pre,])
# Z1 : 19 X 1 matrix (19 years of pre-treatment smoking data for 1 treated state)
Z1 <- as.matrix(Y1[index_pre,1])

# [Z1,Z0]
Z <- cbind(Z1, Z0)




W <- matrix(0, 1, dim(Y)[2])

W[1] <- 1

data <- list('Y' = Y, 'Z' = Z, 'X' = X)

#in this setting, we just use the same Y, Z, X to choose v-weights in synth
tune_params_synth <- list('Y' = Y, 'Z' = Z, 'X' = X)


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
save(fit_elastic_net, file = "smoke_elast_nocov.RData")
writeMat("smoke_elast_nocov.mat", 
         w = fit_elastic_net$w, int = fit_elastic_net$int, 
         Y_est = fit_elastic_net$Y_est, Y_true = fit_elastic_net$Y_true, 
         std_err_i = fit_elastic_net$std_err_i, 
         std_err_t = fit_elastic_net$std_err_t, 
         std_err_it = fit_elastic_net$std_err_it)

# Constrained regression
save(fit_constr_reg, file = "smoke_constr_reg_nocov.RData")
writeMat("smoke_constr_reg_nocov.mat", 
         w = fit_constr_reg$w, int = fit_constr_reg$int, 
         Y_est = fit_constr_reg$Y_est, Y_true = fit_constr_reg$Y_true, 
         std_err_i = fit_constr_reg$std_err_i, 
         std_err_t = fit_constr_reg$std_err_t, 
         std_err_it = fit_constr_reg$std_err_it)

# Best subset 
save(fit_subs, file = "smoke_subs_nocov.RData")
writeMat("smoke_subs_nocov.mat", 
         w = fit_subs$w, int = fit_subs$int, 
         Y_est = fit_subs$Y_est, Y_true = fit_subs$Y_true, 
         std_err_i = fit_subs$std_err_i, 
         std_err_t = fit_subs$std_err_t, 
         std_err_it = fit_subs$std_err_it)

# Synthetic Control
save(fit_synth, file = "smoke_synth_nocov.RData")
writeMat("smoke_synth_nocov.mat", 
         w = fit_synth$w, int = fit_synth$int, 
         Y_est = fit_synth$Y_est, Y_true = fit_synth$Y_true, 
         std_err_i = fit_synth$std_err_i, 
         std_err_t = fit_synth$std_err_t, 
         std_err_it = fit_synth$std_err_it)

# Diff-in-diff
save(fit_diff_in_diff, file = "smoke_did_nocov.RData")
writeMat("smoke_did_nocov.mat", 
         w = fit_diff_in_diff$w, int = fit_diff_in_diff$int, 
         Y_est = fit_diff_in_diff$Y_est, Y_true = fit_diff_in_diff$Y_true, 
         std_err_i = fit_diff_in_diff$std_err_i, 
         std_err_t = fit_diff_in_diff$std_err_t, 
         std_err_it = fit_diff_in_diff$std_err_it)

####################################
########## Counterfactual ##########
####################################
T0 <- T0 <- dim(Z)[1]
T0_co <- 10
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
plot(1970:2000, fit_diff_in_diff$Y_true, type = "l", lty = 2, ylim = c(0, 150), xlim = c(1970,2000), col = "red", main = "California: Smoking per capita", xlab = "Year", ylab = "", las = 1, bty = 'L')
lines(1970:2000, fit_synth$Y_est, lty = 1, col= "blue")
lines(1970:2000, fit_constr_reg$Y_est, lty = 1, col= "green")
lines(1970:2000, fit_elastic_net$Y_est, lty = 1, col= "purple4")
lines(1970:2000, fit_subs$Y_est, lty = 1, col= "orange")
lines(1970:2000, fit_diff_in_diff$Y_est, lty = 1, col= "yellow")
abline(v = 1989, col="black")
abline(v = 1975, col = "grey96")
abline(v = 1980, col = "grey96")
abline(v = 1985, col = "grey96")
abline(v = 1990, col = "grey96")
abline(v = 1995, col = "grey96")
abline(v = 2000, col = "grey96")
legend("topright",legend=c("Actual data", "Original synth.", "Regression w/restrictions",expression(paste("Elastic net (opt. ", lambda," and ",alpha,")" )),"Best subset (opt. k)","Difference-in-Differences"), col=c("red","blue","green","purple4","orange","yellow"),lty=c(2,1,1,1,1,1), ncol=1, bty = 'n', cex = 0.7)
arrows(x0=1987.5, y0=145,x1=1988.5, y1=145, col=c("black"), lwd=1 , length = 0.05,xpd=TRUE)
text(x=1984.5,y=145,pos=4,label = "Policy", cex = 0.65)


### Standard Errors

tau <- cbind(fit_synth$Y_true[20:31]-fit_synth$Y_est[20:31],fit_elastic_net$Y_true[20:31]-fit_elastic_net$Y_est[20:31]) # cbind for each method
std_err <- cbind(fit_synth$std_err_i, fit_elastic_net$std_err_i)

plot(1989:2000, tau[,1], type = "l", lty = 1, ylim = c(-100, 100), xlim = c(1989,2000), col = "blue", main = "California: Standard Errors", xlab = "Year", ylab = "", las = 1, bty = "L")
lines(1989:2000, tau[,1]+1.96*std_err[,1], lty = 3, col= "blue")
lines(1989:2000, tau[,1]-1.96*std_err[,1], lty = 3, col= "blue")
lines(1989:2000, tau[,2], lty = 1, col= "darkmagenta")
lines(1989:2000, tau[,2]+1.96*std_err[,2], lty = 2, col= "darkmagenta")
lines(1989:2000, tau[,2]-1.96*std_err[,2], lty = 2, col= "darkmagenta")
abline(h = 0, col= "black")
abline(v = 1990, col = "grey96")
abline(v = 1991, col = "grey96")
abline(v = 1992, col = "grey96")
abline(v = 1993, col = "grey96")
abline(v = 1994, col = "grey96")
abline(v = 1995, col = "grey96")
abline(v = 1996, col = "grey96")
abline(v = 1997, col = "grey96")
abline(v = 1998, col = "grey96")
abline(v = 1999, col = "grey96")
abline(v = 2000, col = "grey96")
legend("topright",legend=c("ADH synth. treatment","ADH treatment +/-1.96*std.err.",expression(paste("Elastic net treatment (opt. ", lambda," and ",alpha,")" )),"Elastic net treatment +/-1.96*std.err."), col=c("blue","blue","darkmagenta","darkmagenta"),lty=c(1,3,1,2), ncol=1, bty = 'n', cex = 0.65)

### Standard Errors Counterfactual

tau <- cbind(fit_synth_co$Y_true[11:19]-fit_synth_co$Y_est[11:19],fit_elastic_net_co$Y_true[11:19]-fit_elastic_net_co$Y_est[11:19]) # cbind for each method
std_err <- cbind(fit_synth_co$std_err_i, fit_elastic_net_co$std_err_i)

plot(1980:1988, tau[,1], type = "l", lty = 1, ylim = c(-100, 100), xlim = c(1980,1988), col = "blue", main = "California: Standard Errors", xlab = "Year", ylab = "", las = 1, bty = "L")
lines(1980:1988, tau[,1]+1.96*std_err[,1], lty = 3, col= "blue")
lines(1980:1988, tau[,1]-1.96*std_err[,1], lty = 3, col= "blue")
lines(1980:1988, tau[,2], lty = 1, col= "darkmagenta")
lines(1980:1988, tau[,2]+1.96*std_err[,2], lty = 2, col= "darkmagenta")
lines(1980:1988, tau[,2]-1.96*std_err[,2], lty = 2, col= "darkmagenta")
abline(h = 0, col= "black")
abline(v = 1981, col = "grey96")
abline(v = 1982, col = "grey96")
abline(v = 1983, col = "grey96")
abline(v = 1984, col = "grey96")
abline(v = 1985, col = "grey96")
abline(v = 1986, col = "grey96")
abline(v = 1987, col = "grey96")
abline(v = 1988, col = "grey96")
legend("topright",legend=c("ADH synth. treatment","ADH treatment +/-1.96*std.err.",expression(paste("Elastic net treatment (opt. ", lambda," and ",alpha,")" )),"Elastic net treatment +/-1.96*std.err."), col=c("blue","blue","darkmagenta","darkmagenta"),lty=c(1,3,1,2), ncol=1, bty = 'n', cex = 0.65)


## Weights
weights <- cbind(fit_synth$w, fit_constr_reg$w, fit_elastic_net$w, fit_subs$w)

control_names <- c('AL', 'AR', 'CO', 'CT', 'DE', 'GA', 'ID', 'IL', 'IN', 'IA', 'KS', 'KY', 'LA', 'ME', 'MN', 'MS', 'MO', 'MT', 'NE', 'NV', 'NH', 'NM', 'NC', 'ND', 'OH', 'OK', 'PA', 'RI', 'SC', 'SD', 'TN', 'TX', 'UT', 'VT', 'VA', 'WV', 'WI', 'WY')
weights_synth <- as.data.frame(cbind(control_names,weights[,1]))
colnames(weights_synth) <- c("controls","w")
weights_synth$w <- as.numeric(as.character(weights_synth$w))
p <- ggplot(weights_synth, aes(x=controls, y=w))+geom_bar(stat="identity", fill = "steelblue3",color ="black", show.legend = FALSE)+labs(title="",x="", y = "Original synth.")+scale_fill_manual(values = c("steelblue3"))+ylim(-1, 1)+coord_flip()+ theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major.x = element_blank(),panel.grid.major.y = element_line( size=.1, color="grey96" ),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(),
        panel.background = element_blank()) 

weights_constr_reg <- as.data.frame(cbind(control_names,weights[,2]))
colnames(weights_constr_reg) <- c("controls","w")
weights_constr_reg$w <- as.numeric(as.character(weights_constr_reg$w))
p1 <- ggplot(weights_constr_reg, aes(x=controls, y=w))+geom_bar(stat="identity", fill = "springgreen3",color ="black", show.legend = FALSE)+labs(title="",x="", y = "Reg./w.restr.")+scale_fill_manual(values = c("springgreen3"))+ylim(-1, 1)+coord_flip()+ theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major.x = element_blank(),panel.grid.major.y = element_line( size=.1, color="grey96" ),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(),
        panel.background = element_blank()) 

weights_elastic_net <- as.data.frame(cbind(control_names,weights[,3]))
colnames(weights_elastic_net) <- c("controls","w")
weights_elastic_net$w <- as.numeric(as.character(weights_elastic_net$w))
p2 <- ggplot(weights_elastic_net, aes(x=controls, y=w))+geom_bar(stat="identity", fill = "maroon4",color ="black", show.legend = FALSE)+labs(title="",x="", y = "Elastic Net")+scale_fill_manual(values = c("maroon4"))+ylim(-1, 1)+coord_flip()+ theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major.x = element_blank(),panel.grid.major.y = element_line( size=.1, color="grey96" ),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(),
        panel.background = element_blank()) 

weights_subs <- as.data.frame(cbind(control_names,weights[,4]))
colnames(weights_subs) <- c("controls","w")
weights_subs$w <- as.numeric(as.character(weights_subs$w))
p3 <- ggplot(weights_subs, aes(x=controls, y=w))+geom_bar(stat="identity", fill = "orange2",color ="black", show.legend = FALSE)+labs(title="",x="", y = "Best subset")+scale_fill_manual(values = c("orange2"))+ylim(-1, 1)+coord_flip()+ theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major.x = element_blank(),panel.grid.major.y = element_line( size=.1, color="grey96" ),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(),
        panel.background = element_blank()) 

figure <- ggarrange(p, p1, p2,p3,
                    ncol = 4, nrow = 1)
annotate_figure(
  figure,
  top = text_grob("California: Weights", face = "bold", size = 14)
)
