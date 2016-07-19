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

