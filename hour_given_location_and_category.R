library(stats4)
library(lubridate)
source('~/r-workspace/data_incubator/di-DC/feature_selection.R')

HVDensity <- function(x, tau, kappa, lambda) {
  lambada <- lambda^2
  
  Mod <- function(x, m) x - floor(x/m)*m 
  x <- Mod(x, 2*pi)
  x <- x/(2*pi)
  gamma <- 1 - (  lambda/30 + (kappa *lambda)/7 + 0.25 *lambda *tau)
  lambda* (kappa *(x^6 - 5* x^3 + 4.5 *x^2 - 0.5* x) + 
             tau *(3 *x^5 - 10* x^3 + 7.5* x^2 - 0.5 *x) + (x^4 - 2* x^3 + 
                                                              x^2)) + gamma
}

HVLogl <- function(tl, theta) {
  t <- tl[[1]]
  k <- ncol(tl)-1
  tauLin <- theta[seq_len(k+1)]
  kappaLin <- theta[k+1 + seq_len(k+1)]
  lambdaLin <- theta[2*(k+1) + seq_len(k+1)]
  print(head(cbind(tau=tauLin, kappa=kappaLin, lambda=lambdaLin)))
  
  tau <- cbind(1, as.matrix(tl[, -1])) %*% as.matrix(tauLin)
  kappa <- cbind(1, as.matrix(tl[, -1])) %*% as.matrix(kappaLin)
  lambda <- cbind(1, as.matrix(tl[, -1])) %*% as.matrix(lambdaLin)
  
  -sum(log(HVDensity(t, tau, kappa, lambda)))
}

FitHV <- function(tl) {
  ### t from 0 to <2*pi
  k <- ncol(tl)-1
  const <- c(1, rep(0, k))
  fit <- mymle(HVLogl, start = list(theta=rep(const, times = 3)), fixed = list(tl = tl), control = list())
  list(coef = fit@coef, vcov = vcov(fit), logl = logLik(fit), aic = AIC(fit))
}


parCount <- 2
d <- read.csv('~/data-repository/crime/crime_incidents_2013_CSV.csv')
datetime <- mdy_hms(d$START_DATE)
timeOfDay <- (hour(datetime) + minute(datetime)/60) /12*pi
timeLoc <- data.frame(t = timeOfDay, l = d$CENSUS_TRACT, c = d$OFFENSE)
timeLoc <- subset(timeLoc, !is.na(t))
timeLoc <- merge(timeLoc, demog_variates[, 1:(parCount+2)], by.x = 'l', by.y = 'TRACT')
timeLoc$l <- NULL
timeLoc$STATE <- NULL
offenses <- levels(d$OFFENSE)
tl <- subset(timeLoc, c == offenses[3]) # works for 3 but not for 8
tl$c <- NULL
fit <- FitHV(tl)

lower=fit$coef-1.96*sqrt(diag(fit$vcov))
higher=fit$coef+1.96*sqrt(diag(fit$vcov))
any(sign(lower)*sign(higher) > 0)
