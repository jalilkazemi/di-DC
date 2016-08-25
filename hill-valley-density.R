library(ggplot2)
library(stats4)

HVDensity <- function(x, hill, valley, delta) {
  if (abs(delta) >= 1/pi)
    return(NA)
  
  Mod <- function(x, m) x - floor(x/m)*m 
  alpha <- Mod(valley, 2*pi)
  beta <- Mod(hill-valley, 2*pi)
  print(c(alpha=alpha, beta=beta, delta=delta))
  
  b1 <- function(x) (x-beta)^2 * (x+0.5*beta) / (0.5*beta^3)
  b2 <- function(x) x^2 * (x-1.5*beta) / (-0.5*beta^3)

  b3 <- function(x) (x-beta)^2 * (x-3*pi+0.5*beta) / (-4*(pi-0.5*beta)^3)
  b4 <- function(x) (x-2*pi)^2 * (x+pi-1.5*beta) / (4*(pi-0.5*beta)^3)

  eps <- 0.5/pi-0.5*delta

  sft_x <- ifelse(x>alpha, x-alpha, x-alpha+2*pi)
  ifelse(sft_x<beta, eps * b1(sft_x) + (eps+delta) * b2(sft_x),
         eps * b3(sft_x) + (eps+delta) * b4(sft_x))
}

HVLogl <- function(t, hill, valley, delta) -sum(log(HVDensity(t, hill, valley, delta)))

FitHV <- function(t) {
  ### t from 0 to <2*pi
  fit <- mle(HVLogl, start = list(hill=2*pi, valley=pi, delta=0.5/pi), fixed = list(t = t), control = list())
  list(coef = fit@coef, vcov = fit@vcov, logl = fit@min)
}


### Test examples
xx <- seq(0, 2*pi, by=0.01)
yy <- HVDensity(xx, valley=0.25*2*pi, hill=(0.666+0.25)*2*pi, 0.5/pi)
ggplot(data.frame(x=xx, y=yy), aes(x=x, y=y)) + geom_path()
yy <- HVDensity(xx, valley=0.125*2*pi, hill=(0.2+0.125)*2*pi, 0.9/pi)
ggplot(data.frame(x=xx, y=yy), aes(x=x, y=y)) + geom_path()
yy <- HVDensity(xx, valley=0.125*2*pi, hill=(0.2+0.125)*2*pi, 0)
ggplot(data.frame(x=xx, y=yy), aes(x=x, y=y)) + geom_path()
yy <- HVDensity(xx, valley=0, hill=(1+0.1)*2*pi, 0.5/pi)
ggplot(data.frame(x=xx, y=yy), aes(x=x, y=y)) + geom_path()

HVLogl(xx, valley=0.25*2*pi, hill=(0.666+0.25)*2*pi, 0.5/pi)

fit <- FitHV(xx)
yy <- HVDensity(xx, fit$coef[1], fit$coef[2], fit$coef[3])
ggplot(data.frame(x=xx, y=yy), aes(x=x, y=y)) + geom_path()

