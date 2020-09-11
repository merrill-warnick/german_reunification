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
  
data <- list( "X" = X, "Y" = Y, "Z" = Z)


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
## ?? why twice tune_params_synth??

#prepare W:

W <- matrix(0, 1, dim(data$Y)[2])

W[1] <- 1


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
save(fit_elastic_net, file = "germ_elast_nocov.RData")
writeMat("germ_elast_nocov.mat", 
         w = fit_elastic_net$w, int = fit_elastic_net$int, 
         Y_est = fit_elastic_net$Y_est, Y_true = fit_elastic_net$Y_true, 
         std_err_i = fit_elastic_net$std_err_i, 
         std_err_t = fit_elastic_net$std_err_t, 
         std_err_it = fit_elastic_net$std_err_it)

# Constrained regression
save(fit_constr_reg, file = "germ_constr_reg_nocov.RData")
writeMat("germ_constr_reg_nocov.mat", 
         w = fit_constr_reg$w, int = fit_constr_reg$int, 
         Y_est = fit_constr_reg$Y_est, Y_true = fit_constr_reg$Y_true, 
         std_err_i = fit_constr_reg$std_err_i, 
         std_err_t = fit_constr_reg$std_err_t, 
         std_err_it = fit_constr_reg$std_err_it)

# Best subset 
save(fit_subs, file = "germ_subs_nocov.RData")
writeMat("germ_subs_nocov.mat", 
         w = fit_subs$w, int = fit_subs$int, 
         Y_est = fit_subs$Y_est, Y_true = fit_subs$Y_true, 
         std_err_i = fit_subs$std_err_i, 
         std_err_t = fit_subs$std_err_t, 
         std_err_it = fit_subs$std_err_it)

# Synthetic Control
save(fit_synth, file = "germ_synth_nocov.RData")
writeMat("germ_synth_nocov.mat", 
         w = fit_synth$w, int = fit_synth$int, 
         Y_est = fit_synth$Y_est, Y_true = fit_synth$Y_true, 
         std_err_i = fit_synth$std_err_i, 
         std_err_t = fit_synth$std_err_t, 
         std_err_it = fit_synth$std_err_it)

# Diff-in-diff
save(fit_diff_in_diff, file = "germ_did_nocov.RData")
writeMat("germ_did_nocov.mat", 
         w = fit_diff_in_diff$w, int = fit_diff_in_diff$int, 
         Y_est = fit_diff_in_diff$Y_est, Y_true = fit_diff_in_diff$Y_true, 
         std_err_i = fit_diff_in_diff$std_err_i, 
         std_err_t = fit_diff_in_diff$std_err_t, 
         std_err_it = fit_diff_in_diff$std_err_it)

####################################
########## Counterfactual ##########
####################################
T0 <- T0 <- dim(Z)[1]
T0_co <- 21
T1_co <- T0 - T0_co

Y_co = data$Y[1:T0,]
Z_co = data$Z[1:T0_co,]
X_co = data$X

fit_elastic_net_co <- general_estimate(Y_co, Z_co, X_co, W, method = "elastic_net", 
                                    tune_params = list(c(seq(from = 1e-02, to = 1e-01, by = 1e-02),
                                                         seq(from = 2e-01, to = 100, by = 1e-01), 
                                                         seq(from = 200, to = 50000, by = 100)), seq(from = 0.1, to = 0.9, by = 0.1)))

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

# fit training model
synth.out <- 
  synth(
    data.prep.obj=dataprep.out,
    Margin.ipop=.005,Sigf.ipop=7,Bound.ipop=6
  )

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

#Extract the treatment and control units from the dataprep object  
X0_co_synth <- dataprep.out$X0
X1_co_synth <- dataprep.out$X1

Z1_co_synth <- dataprep.out$Z1
Z0_co_synth <- dataprep.out$Z0

Y1_co_synth <- dataprep.out$Y1plot
Y0_co_synth <- dataprep.out$Y0plot

#Re-bind data so that the treated unit is in the first position
# [X1,X0]
X_co_synth <- cbind(X1_co_synth, X0_co_synth)

# [Y1,Y0]
Y_co_synth <- cbind(Y1_co_synth, Y0_co_synth)

# [Z1,Z0]
Z_co_synth <- cbind(Z1_co_synth, Z0_co_synth)

fit_synth_co <- general_estimate(Y_co_synth, Z_co_synth, X_co_synth, W, method = "synth", tune_params = tune_params_synth)

save(fit_elastic_net_co, file = "germ_en_co_nocov.RData")
save(fit_synth_co, file = "germ_synth_co_nocov.RData")

###########################
########## Plots ##########
###########################

### Treatment figure
plot(1960:2003, fit_diff_in_diff$Y_true, type = "l", lty = 2, ylim = c(0, 35000), xlim = c(1960,2003), col = "red", main = "West Germany: per capita GDP", xlab = "Year", ylab = "", las = 1, bty = 'L')
lines(1960:2003, fit_synth$Y_est, lty = 1, col= "blue")
lines(1960:2003, fit_constr_reg$Y_est, lty = 1, col= "green")
lines(1960:2003, fit_elastic_net$Y_est, lty = 1, col= "purple4")
lines(1960:2003, fit_subs$Y_est, lty = 1, col= "orange")
lines(1960:2003, fit_diff_in_diff$Y_est, lty = 1, col= "yellow")
abline(v = 1989, col="black")
abline(v = 1960, col = "grey96")
abline(v = 1970, col = "grey96")
abline(v = 1980, col = "grey96")
abline(v = 1990, col = "grey96")
abline(v = 2000, col = "grey96")
legend("topleft",legend=c("Actual data", "Original synth.", "Regression w/restrictions",expression(paste("Elastic net (opt. ", lambda," and ",alpha,")" )),"Best subset (opt. k)","Difference-in-Differences"), col=c("red","blue","green","purple4","orange","yellow"),lty=c(2,1,1,1,1,1), ncol=1, bty = 'n', cex = 0.6)
arrows(x0=1987, y0=32500,x1=1988.5, y1=32500, col=c("black"), lwd=1 , length = 0.05,xpd=TRUE)
text(x=1979,y=32500,pos=4,label = "Reunification", cex = 0.55)

### Standard Errors
tau <- cbind(fit_synth$Y_true[31:44]-fit_synth$Y_est[31:44],fit_elastic_net$Y_true[31:44]-fit_elastic_net$Y_est[31:44]) # cbind for each method
std_err <- cbind(fit_synth$std_err_i, fit_elastic_net$std_err_i)

plot(1990:2003, tau[,1], type = "l", lty = 1, ylim = c(-12500, 12500), xlim = c(1990,2003), col = "blue", main = "West Germany: Standard Errors", xlab = "Year", ylab = "", las = 1, bty = "L")
lines(1990:2003, tau[,1]+1.96*std_err[,1], lty = 3, col= "blue")
lines(1990:2003, tau[,1]-1.96*std_err[,1], lty = 3, col= "blue")
lines(1990:2003, tau[,2], lty = 1, col= "darkmagenta")
lines(1990:2003, tau[,2]+1.96*std_err[,2], lty = 2, col= "darkmagenta")
lines(1990:2003, tau[,2]-1.96*std_err[,2], lty = 2, col= "darkmagenta")
abline(h = 0, col= "black")
abline(v = 1990, col = "grey96")
abline(v = 1992, col = "grey96")
abline(v = 1994, col = "grey96")
abline(v = 1996, col = "grey96")
abline(v = 1998, col = "grey96")
abline(v = 2000, col = "grey96")
abline(v = 2002, col = "grey96")
legend("topright",legend=c("ADH synth. treatment","ADH treatment +/-1.96*std.err.",expression(paste("Elastic net treatment (opt. ", lambda," and ",alpha,")" )),"Elastic net treatment +/-1.96*std.err."), col=c("blue","blue","darkmagenta","darkmagenta"),lty=c(1,3,1,2), ncol=1, bty = 'n', cex = 0.5)

### Standard Errors Counterfactual
tau <- cbind(fit_synth_co$Y_true[22:30]-fit_synth_co$Y_est[22:30],fit_elastic_net_co$Y_true[22:30]-fit_elastic_net_co$Y_est[22:30]) # cbind for each method
std_err <- cbind(fit_synth_co$std_err_i, fit_elastic_net_co$std_err_i)

plot(1981:1989, tau[,1], type = "l", lty = 1, ylim = c(-12500, 12500), xlim = c(1981,1989), col = "blue", main = "West Germany: Counterfactual Treatment", xlab = "Year", ylab = "", las = 1, bty = "L")
lines(1981:1989, tau[,1]+1.96*std_err[,1], lty = 3, col= "blue")
lines(1981:1989, tau[,1]-1.96*std_err[,1], lty = 3, col= "blue")
lines(1981:1989, tau[,2], lty = 1, col= "darkmagenta")
lines(1981:1989, tau[,2]+1.96*std_err[,2], lty = 2, col= "darkmagenta")
lines(1981:1989, tau[,2]-1.96*std_err[,2], lty = 2, col= "darkmagenta")
abline(h = 0, col= "black")
abline(v = 1982, col = "grey96")
abline(v = 1983, col = "grey96")
abline(v = 1984, col = "grey96")
abline(v = 1985, col = "grey96")
abline(v = 1986, col = "grey96")
abline(v = 1987, col = "grey96")
abline(v = 1988, col = "grey96")
abline(v = 1989, col = "grey96")
legend("topright",legend=c("ADH synth. treatment","ADH treatment +/-1.96*std.err.",expression(paste("Elastic net treatment (opt. ", lambda," and ",alpha,")" )),"Elastic net treatment +/-1.96*std.err."), col=c("blue","blue","darkmagenta","darkmagenta"),lty=c(1,3,1,2), ncol=1, bty = 'n', cex = 0.5)

# Weights
weights <- cbind(fit_synth$w, fit_constr_reg$w, fit_elastic_net$w, fit_subs$w)
#theme_set(theme_bw())

control_names <- c("USA", "GBR", "AUT", "BEL", "DNK", "FRA", "ITA", "NLD", "NOR", "CHE", "JPN", "GRC", "PRT", "ESP","AUS","NZL")
weights_synth <- as.data.frame(cbind(control_names,weights[,1]))
colnames(weights_synth) <- c("controls","w")
weights_synth$w <- as.numeric(as.character(weights_synth$w))
p <- ggplot(weights_synth, aes(x=controls, y=w))+geom_bar(stat="identity", fill = "blue",color ="black", show.legend = FALSE)+labs(title="",x="", y = "Original synth.")+scale_fill_manual(values = c("royalblue"))+ylim(-1, 1)+coord_flip()+ theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major.x = element_blank(),panel.grid.major.y = element_line( size=.1, color="grey96" ),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(),
        panel.background = element_blank())

weights_constr_reg <- as.data.frame(cbind(control_names,weights[,2]))
colnames(weights_constr_reg) <- c("controls","w")
weights_constr_reg$w <- as.numeric(as.character(weights_constr_reg$w))
p1 <- ggplot(weights_constr_reg, aes(x=controls, y=w))+geom_bar(stat="identity", fill = "forestgreen",color ="black", show.legend = FALSE)+labs(title="",x="", y = "Reg./w.restr.")+scale_fill_manual(values = c("forestgreen"))+ylim(-1, 1)+coord_flip()+ theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major.x = element_blank(),panel.grid.major.y = element_line( size=.1, color="grey96" ),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(),
        panel.background = element_blank())

weights_elastic_net <- as.data.frame(cbind(control_names,weights[,3]))
colnames(weights_elastic_net) <- c("controls","w")
weights_elastic_net$w <- as.numeric(as.character(weights_elastic_net$w))
p2 <- ggplot(weights_elastic_net, aes(x=controls, y=w))+geom_bar(stat="identity", fill = "purple",color ="black", show.legend = FALSE)+labs(title="",x="", y = "Elastic Net")+scale_fill_manual(values = c("purple"))+ylim(-1, 1)+coord_flip()+ theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major.x = element_blank(),panel.grid.major.y = element_line( size=.1, color="grey96" ),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(),
        panel.background = element_blank())

weights_subs <- as.data.frame(cbind(control_names,weights[,4]))
colnames(weights_subs) <- c("controls","w")
weights_subs$w <- as.numeric(as.character(weights_subs$w))
p3 <- ggplot(weights_subs, aes(x=controls, y=w))+geom_bar(stat="identity", fill = "orange",color ="black", show.legend = FALSE)+labs(title="",x="", y = "Best subset")+scale_fill_manual(values = c("orange"))+ylim(-1, 1)+coord_flip()+ theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major.x = element_blank(),panel.grid.major.y = element_line( size=.1, color="grey96" ),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(),
        panel.background = element_blank()) 

figure <- ggarrange(p, p1, p2,p3,
                    ncol = 4, nrow = 1)
annotate_figure(
  figure,
  top = text_grob("West Germany: Weigths", color = "black", face = "bold", size = 14),
)
figure

