library(survival)

library(rms)

data(package="survival")
dd<-datadist(lung)

options(datadist="dd")

f <- cph(Surv(time, status) ~age+sex +ph.karno, data = lung, x=T, y = T,surv=T)

survival <- Survival(f)
survival1 <- function(x) survival(365, x)
survival2 <- function(x) survival(730, x)


nom <- nomogram(f, fun = list(survival1, survival2),
                fun.at =c(0.05,seq(0.1,0.9, by =0.05), 0.95),
                funlabel = c('1 year survival', '2 year survival'))

plot(nom)


cal <- calibrate(f, cmethod = "KM", method = "boot", u = 365, m = 38, B = 228)
cal2 <- calibrate(f, cmethod = "KM", method = "boot", u = 730, m = 38, B = 228)
## Using Cox survival estimates at 120 Days
cal
## calibrate.cph(fit = cph, cmethod = "KM", method = "boot", u = 365, 
##     m = 38, B = 228)
## 
## 
## n=226  B=228  u=365 Day
## 
##      index.orig     training         test mean.optimism mean.corrected   n
## [1,] -0.4836310 -0.002478164 -0.003134677  0.0006565131     -0.4842875 228
## [2,] -0.5180526  0.010037703 -0.058465514  0.0685032169     -0.5865558  28
## [3,] -0.4576156  0.075069139  0.053964406  0.0211047329     -0.4787204  11
## [4,] -0.3499039          NaN          NaN           NaN            NaN   0
## [5,] -0.3336164          NaN          NaN           NaN            NaN   0
##      mean.predicted        KM KM.corrected   std.err
## [1,]      0.7422109 0.2585799    0.2579234 0.2736880
## [2,]      0.8120250 0.2939724    0.2254692 0.2582661
## [3,]      0.8530646 0.3954489    0.3743442 0.1908749
## [4,]      0.8882504 0.5383464          NaN 0.1517891
## [5,]      0.9229901 0.5893736          NaN 0.1466856
plot(cal, lwd = 2, lty = 1, #errbar.col = c(rgb(0, 118, 192, maxColorValue = 255)),
     #xlab = "Nomogram-Predicted Probability of 1-Year OS", ylab = "Actual 1-Year DFS (proportion)",
     #col = c(rgb(192, 98, 83, maxColorValue = 255)), subtitles = FALSE, 
     #xlim = c(0,1), ylim = c(0, 1), main = "Calibrate plot"
     )
lines(cal[, c("mean.predicted", "KM")], type = "l", lwd = 2, col = c(rgb(192, 98,
                                                                         83, maxColorValue = 255)), pch = 16)
abline(0, 1, lty = 3, lwd = 2, col = c(rgb(0, 118, 192, maxColorValue = 255)))


plot(cal, lwd = 2, lty = 1, errbar.col = c(rgb(0, 118, 192, maxColorValue = 255)),
     xlab = "Nomogram-Predicted Probability of 1-Year OS", ylab = "Actual 1-Year DFS (proportion)",
     col = c(rgb(192, 98, 83, maxColorValue = 255)), subtitles = FALSE, xlim = c(0,
                                                                                 1), ylim = c(0, 1), main = "Calibrate plot")
