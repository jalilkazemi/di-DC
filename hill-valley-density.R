library(ggplot2)
library(stats4)

HVDensity <- function(x, alpha, beta, delta) {
  print(c(alpha=alpha, beta=beta, delta=delta))
  if (alpha < 0 || alpha >= 2*pi || beta <= 0 || beta >= 2*pi || abs(delta) >= 1/pi)
    return(NA)
  
  b1 <- function(x) (x-beta)^2 * (x+0.5*beta) / (0.5*beta^3)
  b2 <- function(x) x^2 * (x-1.5*beta) / (-0.5*beta^3)

  b3 <- function(x) (x-beta)^2 * (x-3*pi+0.5*beta) / (-4*(pi-0.5*beta)^3)
  b4 <- function(x) (x-2*pi)^2 * (x+pi-1.5*beta) / (4*(pi-0.5*beta)^3)

  eps <- 0.5/pi-0.5*delta

  sft_x <- ifelse(x>alpha, x-alpha, x-alpha+2*pi)
  ifelse(sft_x<beta, eps * b1(sft_x) + (eps+delta) * b2(sft_x),
         eps * b3(sft_x) + (eps+delta) * b4(sft_x))
}

HVLogl <- function(t, alpha, beta, delta) -sum(log(HVDensity(t, alpha, beta, delta)))

FitHV <- function(t) {
  ### t from 0 to <2*pi
  fit <- mle(HVLogl, start = list(alpha=pi, beta=pi, delta=0.5/pi), fixed = list(t = t), control = list())
  list(coef = fit@coef, vcov = fit@vcov, logl = fit@min)
}


### Test examples
xx <- seq(0, 2*pi, by=0.01)
yy <- HVDensity(xx, 0.25*2*pi, 0.666*2*pi, 0.5/pi)
ggplot(data.frame(x=xx, y=yy), aes(x=x, y=y)) + geom_path()
yy <- HVDensity(xx, 0.125*2*pi, 0.2*2*pi, 1/pi)
ggplot(data.frame(x=xx, y=yy), aes(x=x, y=y)) + geom_path()
yy <- HVDensity(xx, 0.125*2*pi, 0.2*2*pi, 0)
ggplot(data.frame(x=xx, y=yy), aes(x=x, y=y)) + geom_path()

HVLogl(xx, 0.25*2*pi, 0.666*2*pi, 0.5/pi)

fit <- FitHV(xx)
yy <- HVDensity(xx, fit$coef[1], fit$coef[2], fit$coef[3])
ggplot(data.frame(x=xx, y=yy), aes(x=x, y=y)) + geom_path()

# eps * (x-beta)^2 * (x+0.5*beta) / (0.5*beta^3) +
#   (eps+delta) * x^2 * (x-1.5*beta) / (-0.5*beta^3) +
#   eps * (x-beta)^2 * (x-3*pi+0.5*beta) / (-4*(pi-0.5*beta)^3) +
#   (eps+delta) * (x-2*pi)^2 * (x+pi-1.5*beta) / (4*(pi-0.5*beta)^3)
# 
# 
# eps * (x-beta)^2 * (x+0.5*beta) / (0.5*beta^3) + 0.5*(eps+delta) * (x-2*pi)^2 * (x+pi-1.5*beta) / (4*(pi-0.5*beta)^3)
# 0.5*(eps+delta) * x^2 * (x-1.5*beta) / (-0.5*beta^3) + eps * (x-beta)^2 * (x-3*pi+0.5*beta) / (-4*(pi-0.5*beta)^3)
# 
# 8*pi^4*eps/beta^3 - 8*pi^3*eps/beta^2 + 4*pi^3*(pi-beta)*(delta+eps)/(2*pi-beta)^3+2*pi*eps
# 2*pi*(-beta^3 + 6*pi*beta^2 - 8*pi^2*beta + 4*pi^3)*eps/(2*pi-beta)^3 - 4*pi^3*(pi-beta)*(delta+eps)/beta^3
# 
# 4*pi^3*(beta*(delta+eps)-pi*(delta-eps))/beta^3 + (4*pi-8*pi^3/beta^2)*eps + 4*pi^3*(beta-pi)*(delta-eps)/(beta-2*pi)^3
