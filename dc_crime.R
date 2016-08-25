library(ggplot2)
library(RANN)
library(lubridate)
# library(stringr)
library(sp)
library(corrplot)
# library(ggmap)

d <- read.csv('~/Downloads/crime_incidents_2013_CSV.csv')

# blockMatch <- str_match(d$BLOCKSITEADDRESS, '([0-9]+) +((- +)?[0-9]+ +)?BLOCK OF (.+) (SW|NW|SE|NE)?')
# blockAddress <- str_c(blockMatch[, 2], blockMatch[, 5], 'DC', sep = ' ')
# crossMatch <- str_match(d$BLOCKSITEADDRESS, '(.+) +(SW|NW|SE|NE) +(SOUTHBOUND +)?AND +(.+) +(SW|NW|SE|NE)')
# crossAddress <- str_c(crossMatch[, 2], crossMatch[, 5], ',DC', sep = ' ')
# address <- ifelse(!is.na(blockAddress), blockAddress, crossAddress)

HeatMap <- function(coord, grid, kn = 10, subtitle = '') {
  query <- grid
  nearestCrimes <- nn2(coord, query, k = kn)
  radiusInd <- max.col(nearestCrimes$nn.dists)
  radius <- nearestCrimes$nn.dists[cbind(seq_along(radiusInd), radiusInd)]
  area <- pi * radius^2
  crimeMap <- cbind(query, density = kn/nrow(coord)/area)
  plt <- ggplot(crimeMap, aes(x = x, y = y, fill = rank(density)))+ 
    geom_tile(color="white", size=0) + 
    coord_equal() + 
    labs(title = sprintf('Crime density inside DC (%d data points)\n%s', nrow(coord), subtitle))
  return(list(plot = plt, data = crimeMap$density))
}

gridSize <- 50
grid <- expand.grid(x = seq(min(d$BLOCKXCOORD), max(d$BLOCKXCOORD), length.out = gridSize), y = seq(min(d$BLOCKYCOORD), max(d$BLOCKYCOORD), length.out = gridSize))
chl <- chull(d$BLOCKXCOORD, d$BLOCKYCOORD)
isInDC <- as.logical(point.in.polygon(grid$x, grid$y, d$BLOCKXCOORD[c(chl, chl[1])], d$BLOCKYCOORD[c(chl, chl[1])]))
grid <- grid[isInDC, ]

coord <- data.frame(x = d$BLOCKXCOORD, y =d$BLOCKYCOORD)
plt <- HeatMap(coord, grid, 10)
plt$plot

pltLst <- list()
offenses <- levels(d$OFFENSE)
crimeRank <- matrix(0, nrow = nrow(grid), ncol = length(offenses))
for (i in seq_along(offenses)) {
  cat(offenses[i], '\n')
  coord <- subset(data.frame(x = d$BLOCKXCOORD, y =d$BLOCKYCOORD), d$OFFENSE == offenses[i])
  pltLst[[i]] <- HeatMap(coord, grid, 10, offenses[i])
  crimeRank[, i] <- pltLst[[i]]$data
}

pltLst[[1]]$plot
pltLst[[2]]$plot
pltLst[[3]]$plot
pltLst[[4]]$plot
pltLst[[5]]$plot
pltLst[[6]]$plot
pltLst[[7]]$plot
pltLst[[8]]$plot
pltLst[[9]]$plot

colnames(crimeRank) <- str_sub(offenses, 1, 10)
M <- cor(crimeRank)
corrplot(M, method = 'circle', type = 'lower')

Density <- function(tm, grid, subtitle = '') {
  query <- grid
  tm <- tm[!is.na(tm)]
  tmCyc <- sort(c(tm, tm[tm<1]+24, tm[tm>23]-24))
  radius <- 1
  nearestCrimesCount <- sapply(query, function(x) sum(abs(x - tmCyc) < radius))
  area <- 2*radius
  crimeTiming <- data.frame(t = query, density = nearestCrimesCount/length(tm)/area)
  plt <- ggplot(crimeTiming, aes(x = t, y = density))+ 
    geom_point() +
    labs(title = sprintf('Crime density inside DC (%d data points)\n%s', length(tm), subtitle))
  return(list(grid = query, plot = plt, data = crimeTiming$density))
}

grid <- seq(0, 24, by = 0.1)

datetime <- mdy_hms(d$START_DATE)
timeOfDay <- hour(datetime) + minute(datetime)/60
plt <- Density(timeOfDay, grid)
plt$plot

pltLst <- list()
offenses <- levels(d$OFFENSE)
crimeDensity <- matrix(0, nrow = length(grid), ncol = length(offenses))
for (i in seq_along(offenses)) {
  cat(offenses[i], '\n')
  tm <- timeOfDay[d$OFFENSE == offenses[i]]
  pltLst[[i]] <- Density(tm, grid, offenses[i])
  crimeDensity[, i] <- pltLst[[i]]$data
}

pltLst[[1]]$plot
pltLst[[2]]$plot
pltLst[[3]]$plot
pltLst[[4]]$plot
pltLst[[5]]$plot
pltLst[[6]]$plot
pltLst[[7]]$plot
pltLst[[8]]$plot
pltLst[[9]]$plot

rob <- data.frame(t = grid, density = pltLst[[6]]$data, offense = 'robbery')
burg <- data.frame(t = grid, density = pltLst[[3]]$data, offense = 'burglary')
ggplot(rbind(rob, burg), aes(x = t, y = density)) + geom_point(aes(color = factor(offense))) +
  labs(title = sprintf('Crime timing inside DC (%d data points)', sum(!is.na(timeOfDay))))

#### Hill Valley density fitting
library(lubridate)

d <- read.csv('~/Downloads/crime_incidents_2013_CSV.csv')
datetime <- mdy_hms(d$START_DATE)
timeOfDay <- (hour(datetime) + minute(datetime)/60) /12*pi
offenses <- levels(d$OFFENSE)

plt <- list()
for (i in seq_along(offenses)) {
  tm <- timeOfDay[d$OFFENSE == offenses[i]]
  tm <- tm[!is.na(tm)]
  xx <- seq(0, 2*pi, by=0.01)
  fit <- NULL
  param <- NULL
  try({
    fit <- FitHV(tm)
    sd <- sqrt(diag(fit$vcov))
    yy <- HVDensity(xx, fit$coef[1], fit$coef[2], fit$coef[3])
    param <- ggplot(data.frame(x=xx, y=yy), aes(x=x, y=y)) +
      geom_path() + ylim(0, max(yy))
    }, silent = TRUE)
  plt[[i]] <- list(size = length(tm),
                   coef = fit$coef,
                   sd = sd,
                   param = param,
                   nonparam = ggplot(data.frame(t=c(-2*pi+tm, tm, 2*pi+tm)), aes(t)) +
                     geom_density() + geom_vline(xintercept = c(0, 2*pi), colour="green", linetype = "longdash"))
}

for (i in seq_along(offenses)) {
  upper <- plt[[i]]$coef + 1.96 * plt[[i]]$sd
  lower <- plt[[i]]$coef - 1.96 * plt[[i]]$sd
  cat(sprintf('Size of %s: %d\n', offenses[i], plt[[i]]$size))
  cat(sprintf('[%f, %f, %f]', lower, plt[[i]]$coef, upper))
  cat('\n')
}

multiplot(plt[[1]]$nonparam, plt[[1]]$param,
          plt[[2]]$nonparam, plt[[2]]$param,
          plt[[3]]$nonparam, plt[[3]]$param,
          plt[[4]]$nonparam, plt[[4]]$param,
          plt[[5]]$nonparam, plt[[5]]$param
          , cols = 5)
multiplot(plt[[6]]$nonparam, plt[[6]]$param,
          plt[[7]]$nonparam, plt[[7]]$param,
          plt[[8]]$nonparam, plt[[8]]$param,
          plt[[9]]$nonparam, plt[[9]]$param
          , cols = 4)

### time density of a particular crime in different locations
coord <- data.frame(x = d$BLOCKXCOORD, y =d$BLOCKYCOORD)
coordSub <- subset(coord, d$OFFENSE == offenses[5])
tm <- timeOfDay[d$OFFENSE == offenses[5]]
coordSub <- subset(coordSub, !is.na(tm))
tm <- tm[!is.na(tm)]
gridSize <- 3
coordSubCoarse <- data.frame(x = cut(coordSub$x, breaks = seq(min(coord$x), max(coord$x), length.out = gridSize+1), include.lowest = TRUE),
                             y = cut(coordSub$y, breaks = seq(min(coord$y), max(coord$y), length.out = gridSize+1), include.lowest = TRUE))
PlotDensity <- function(t) ggplot(data.frame(t=c(-2*pi+t, t, 2*pi+t)), aes(t)) +
  geom_density() + geom_vline(xintercept = c(0, 2*pi), colour="green", linetype = "longdash")

plt <- tapply(seq_len(length(tm)), list(coordSubCoarse$x, coordSubCoarse$y), function(x) PlotDensity(tm[x]))
multiplot(plt[[1]], plt[[2]], plt[[3]], plt[[4]], plt[[5]], plt[[6]], plt[[7]], plt[[8]], plt[[9]], cols = 3)
table(coordSubCoarse$x, coordSubCoarse$y)

plt <- tapply(seq_len(length(tm)), list(coordSubCoarse$x, coordSubCoarse$y), function(x) {
  size <- length(x)
  coef <- NA
  lower <- NA
  upper <- NA
  try({
    fit <- FitHV(tm[x])
    coef <- fit$coef
    sd <- sqrt(diag(fit$vcov))
    upper <- coef + 1.96 * sd
    lower <- coef - 1.96 * sd
  }, silent = TRUE)
  list(size=size, coef=coef, lower=lower, upper=upper)
})

fitregion <- data.frame()
for (i in seq_len(gridSize)) {
  for (j in seq_len(gridSize)) {
    fitregion <- rbind(fitregion, data.frame(size=plt[[i, j]]$size,
               alpha=plt[[i, j]]$coef[1],
               beta=plt[[i, j]]$coef[2],
               delta=plt[[i, j]]$coef[3]))
    cat(sprintf('Size of (%d,%d): %d\n', i, j, plt[[i, j]]$size))
    cat(sprintf('[%f, %f, %f]', plt[[i, j]]$lower, plt[[i, j]]$coef, plt[[i, j]]$upper))
    cat('\n\n')
  }
}
plot(fitregion$size, fitregion$alpha)
plot(fitregion$size, fitregion$alpha+fitregion$beta)
