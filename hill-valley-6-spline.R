HVDensity <- function(x, tau, kappa, lambda) {
  # print(c(tau=tau, kappa=kappa, lambda=lambda))
  lambada <- lambda^2
  
  Mod <- function(x, m) x - floor(x/m)*m 
  x <- Mod(x, 2*pi)
  x <- x/(2*pi)
  gamma <- 1 - (  lambda/30 + (kappa *lambda)/7 + 0.25 *lambda *tau)
lambda* (kappa *(x^6 - 5* x^3 + 4.5 *x^2 - 0.5* x) + 
            tau *(3 *x^5 - 10* x^3 + 7.5* x^2 - 0.5 *x) + (x^4 - 2* x^3 + 
                                                        x^2)) + gamma
}

HVLogl <- function(t, tau, kappa, lambda) -sum(log(HVDensity(t, tau, kappa, lambda)))

FitHV <- function(t) {
  ### t from 0 to <2*pi
  fit <- mle(HVLogl, start = list(tau=1, kappa=1, lambda=1), fixed = list(t = t), control = list())
  list(coef = fit@coef, vcov = fit@vcov, logl = fit@min)
}

xx <- seq(0, 2*pi, by=0.01)
yy <- HVDensity(xx, 0.2*2*pi, 0, 0.5/pi)
plot(xx, yy, type='l')

HVLogl(xx, 0.25*2*pi, (0.666+0.25)*2*pi, 0.5/pi)

fit <- FitHV(xx)
yy <- HVDensity(xx, fit$coef[1], fit$coef[2], fit$coef[3])
ggplot(data.frame(x=xx, y=yy), aes(x=x, y=y)) + geom_path()
