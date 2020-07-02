## Clear
rm(list = ls())

library(R.matlab)
library(shape)
data <- readMat('germ_did_nocov_05.mat')

### Treatment figure
plot(1960:2003, data$Y.true, type = "l", lty = 2, ylim = c(0, 35000), xlim = c(1960,2003), col = "red", main = "West Germany: per capita GDP", xlab = "Year", ylab = "", las = 1, bty = 'L')
lines(1960:2003, data$Y.did, lty = 1, col= "yellow")
abline(v = 1989, col="black")
abline(v = 1960, col = "grey96")
abline(v = 1970, col = "grey96")
abline(v = 1980, col = "grey96")
abline(v = 1990, col = "grey96")
abline(v = 2000, col = "grey96")
legend("topleft",legend=c("Actual data","Difference-in-Differences"), col=c("red","yellow"),lty=c(1,1), ncol=1, bty = 'n', cex = 0.8)
arrows(x0=1987, y0=32500,x1=1988, y1=32499, col=c("black"), lwd=1 , length = 0.05,xpd=TRUE)
text(x=1981.5,y=32500,pos=4,label = "Reunification", cex = 0.5)

### Standard Errors
tau <- cbind(data$Y.true[31:44]-data$Y.did[31:44]) # cbind for each method
std_err <- cbind(data$std.err.did.i)

plot(1990:2003, tau, type = "l", lty = 1, ylim = c(-12500, 12500), xlim = c(1990,2003), col = "blue", main = "West Germany: Standard Errors", xlab = "Year", ylab = "", las = 1, bty = "L")
lines(1990:2003, tau+1.96*std_err, lty = 2, col= "blue")
lines(1990:2003, tau-1.96*std_err, lty = 2, col= "blue")
abline(h = 0, col= "black")
abline(v = 1990, col = "grey96")
abline(v = 1992, col = "grey96")
abline(v = 1994, col = "grey96")
abline(v = 1996, col = "grey96")
abline(v = 1998, col = "grey96")
abline(v = 2000, col = "grey96")
abline(v = 2002, col = "grey96")
legend("topright",legend=c("ADH synth. treatment","ADH treatment +/-1.96*std.err."), col=c("blue","blue"),lty=c(1,2), ncol=1, bty = 'n')

# For grey lines for years add abline(v = 1960,col = "grey96")
