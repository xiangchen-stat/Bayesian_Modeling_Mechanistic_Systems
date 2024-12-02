library(invgamma)
library(LaplacesDemon)

tt <- 10000
a <- rgamma(tt, shape = 2, rate = 3)
b <- invgamma::rinvgamma(tt, shape = 2, rate = 3) # error in package: rate should be scale
c <- LaplacesDemon::rinvgamma(tt, shape = 2, scale = 3)
plot(x = sort(1/a), y = sort(b))
abline(a=0, b=1, col = "red")
plot(x = sort(1/a), y = sort(c))
abline(a=0, b=1, col = "red")