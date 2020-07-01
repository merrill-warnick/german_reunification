## Clear
rm(list = ls())

library(R.matlab)
data <- readMat('germ_did_nocov_05.mat')

### Treatment figure
plot(1960:2003, data$Y.true, type = "l", lty = 2, ylim = c(0, 35000), xlim = c(1960,2003), col = "red", main = "West Germany: per capita GDP", xlab = "Year", ylab = "", las = 1)
lines(1960:2003, data$Y.did, lty = 1, col= "yellow")
abline(v = 1989, col="black")
legend("topleft",legend=c("Actual data","Difference-in-Differences"), col=c("red","yellow"),lty=c(1,1), ncol=1, bty = 'n')
arrows(x0=1988, y0=32500, col=c("black"), lwd=2, xpd=TRUE)
text(x=1986,y=32500,pos=4,label = "Reunification")
box(which = "plot", bty = "l")