## Clear
rm(list = ls())

library(R.matlab)
library(shape)

source('functions.R')
##############################
########## Dataprep ##########
##############################
# Here or in general_estimate function??


###################################################
########## Estimates and Standard Errors ##########
###################################################

# Get estimates and standard errors

# Save matrices for future reference and plots

###########################
########## Plots ##########
###########################

## Load data -> or can also simply call from before
data_did <- readMat('germ_did_nocov_05.mat')
data_elast <- readMat('germ_elast_nocov_05.mat')
data_subset <- readMat('germ_subs_nocov_05.mat')
data_synth <- readMat('germ_synth_05.mat')

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
tau_co <- cbind(data_elast$Y.true.co[(data_elast$T0.co+1):(data_elast$T0.co+data_elast$T1.co)]-data_synth$Y.synth.co[(data_elast$T0.co+1):(data_elast$T0.co+data_elast$T1.co)],data_synth$Y.true.co[(data_elast$T0.co+1):(data_elast$T0.co+data_elast$T1.co)]-data_elast$Y.elast.co[(data_elast$T0.co+1):(data_elast$T0.co+data_elast$T1.co)]) # cbind for each method
std_err_co <- cbind(data_synth$std.err.synth.i.co, data_elast$std.err.elast.i.co)

plot(1981:1989, tau_co[,1], type = "l", lty = 1, ylim = c(-12500, 12500), xlim = c(1981,1989), col = "blue", main = "West Germany: Counterfactual", xlab = "Year", ylab = "", las = 1, bty = "L")
lines(1981:1989, tau_co[,1]+1.96*std_err_co[,1], lty = 3, col= "blue")
lines(1981:1989, tau_co[,1]-1.96*std_err_co[,1], lty = 3, col= "blue")
lines(1981:1989, tau_co[,2], lty = 1, col= "plum2")
lines(1981:1989, tau_co[,2]+1.96*std_err_co[,2], lty = 2, col= "plum2")
lines(1981:1989, tau_co[,2]-1.96*std_err_co[,2], lty = 2, col= "plum2")
abline(h = 0, col= "black")
abline(v = 1982, col = "grey96")
abline(v = 1983, col = "grey96")
abline(v = 1984, col = "grey96")
abline(v = 1985, col = "grey96")
abline(v = 1986, col = "grey96")
abline(v = 1987, col = "grey96")
abline(v = 1988, col = "grey96")
abline(v = 1989, col = "grey96")
legend("topright",legend=c("ADH synth. treatment","ADH treatment +/-1.96*std.err.",expression(paste("Elastic net treatment (opt. ", lambda,"and ",alpha,")" )),"Elastic net treatment +/-1.96*std.err."), col=c("blue","blue","plum2","plum2"),lty=c(1,2,1,2), ncol=1, bty = 'n', cex = 0.65)

## Weights
weights <- cbind(data_synth$w.synth, data_elast$w.elast, data_subset$w.subs)
barplot(weights[,1], horiz = T, xlim = c(-1,1), xlab = "Original Synthetic Control", names.arg = c("USA", "GBR", "AUT", "BEL", "DNK", "FRA", "ITA", "NLD", "NOR", "CHE", "JPN", "GRC", "PRT", "ESP","AUS","NZL"), col = "royalblue", las = 1, cex.axis = 0.8, cex.names = 0.8)
box(which = "plot", lty = "solid",tck = -0.5)
at_tick <- seq(0,16, by = 2)

# Try 2
library(ggplot2) 
theme_set(theme_bw())

control_names <- c("USA", "GBR", "AUT", "BEL", "DNK", "FRA", "ITA", "NLD", "NOR", "CHE", "JPN", "GRC", "PRT", "ESP","AUS","NZL")
weights_synth <- as.data.frame(cbind(control_names,weights[,3]))
colnames(weights_synth) <- c("controls","w")
weights_synth$w <- as.numeric(as.character(weights_synth$w))
p <- ggplot(weights_synth, aes(x=controls, y=w))+geom_bar(stat="identity", fill = "blue",color ="black", show.legend = FALSE)+labs(title="",x="", y = "Original synth.")+scale_fill_manual(values = c("royalblue"))+ylim(-1, 1)+coord_flip()
p
