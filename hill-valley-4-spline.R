HVDensity <- function(x, hill, valley, eps, delta) {
  if (abs(delta) >= 1/pi)
    return(NA)
  
  Mod <- function(x, m) x - floor(x/m)*m 
  alpha <- Mod(valley, 2*pi)
  beta <- Mod(hill-valley, 2*pi)
  print(c(alpha=alpha, beta=beta, delta=delta))

  x <- ifelse(x>alpha, x-alpha, x-alpha+2*pi)

  left <- function(x) 
    eps + 
    (2 *(beta^3* delta - beta^2* delta *pi - 18* beta* delta* pi^2 + 24 *delta *pi^3)* x^2)/
    (beta^2 *(beta - 2 *pi)^2 *(beta + 2 *pi)) -
    (4 *(beta^2 *delta* pi - 14 *beta *delta *pi^2 + 16 *delta *pi^3)* x^3)/
    (beta^3 *(beta - 2 *pi)^2 *(beta + 2* pi)) -
    ((beta^3* delta - 4 *beta^2 *delta *pi + 24 *beta *delta *pi^2 - 24 *delta *pi^3) *x^4)/
    (beta^4 *(beta - 2* pi)^2 *(beta + 2 *pi))
  
  right <- function(x)
    -((-beta^3 *eps + 2 *beta^2* eps *pi + 8* beta *delta *pi^2 + 4* beta* eps *pi^2 - 8 *delta* pi^3 - 8* eps *pi^3)/
        ((beta - 2* pi)^2 *(beta + 2 *pi))) + 
    (2 *(beta^2* delta - beta *delta* pi + 6 *delta *pi^2)* x^2)/
    (beta *(beta - 2 *pi)^2 *(beta + 2* pi)) - 
    (4* delta *pi *x^3)/
    (beta^2 *(beta - 2 *pi)^2) - 
    ((beta *delta - 4 *delta* pi) *x^4)/
    (beta^2 *(beta - 2* pi)^2 *(beta + 2 *pi))
  
  sft_x <- ifelse(x>alpha, x-alpha, x-alpha+2*pi)
  ifelse(sft_x<beta, left(sft_x), right(sft_x))
}

xx <- seq(0, 2*pi, by=0.01)
yy <- HVDensity(xx, 0.2*2*pi, 0, 0, 0.5/pi)
plot(xx, yy, type='l')
