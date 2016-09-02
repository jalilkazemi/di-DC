### Select demographics from ACS
source('~/r-workspace/ACS/acs-extract.R')
library(stringr)

locations <- data.frame(state = 'DC', county = 1)

varLst <- read.csv('~/data-repository/ACS/lookup/Sequence_Number_and_Table_Number_Lookup.txt', as.is = TRUE)
names(varLst) <- c('file', 'table', 'seqNo', 'line', 'startCol', 'totCell1', 'totCell2', 'title', 'subject')
varLst <- subset(varLst, str_trim(varLst$totCell1) != '')
varLst$totCell1 <- as.numeric(str_replace(varLst$totCell1, ' *CELL(S)?', ''))
varLst <- subset(varLst, select = c('table', 'seqNo', 'totCell1'))
names(varLst)[3] <- 'colCount'

VariableGroup <- function(table, seqNo, lines) {
  n <- length(lines)
  data.frame(table = rep(table, n),
             seqNo = rep(seqNo, n),
             line = lines,
             stringsAsFactors = FALSE)
}

rowCount <- sum(varLst$colCount)
variables <- data.frame(table = rep('', rowCount), seqNo = rep(0, rowCount), line = rep(0, rowCount), stringsAsFactors = FALSE)
r <- 0
for (i in seq_len(nrow(varLst))) {
  v <- VariableGroup(varLst$table[i],
                varLst$seqNo[i],
                seq_len(varLst$colCount[i]))
  index <- r+seq_len(varLst$colCount[i])
  variables[index, ] <- v
  
  r <- r+varLst$colCount[i]
}
variables <- as.data.frame(variables)

selVars <- head(variables, n = 10000)
demographics <- ExtractCensusData(selVars, data.frame(state = 'DC', county = 1, stringsAsFactors = FALSE))

### Principal components of demographics
demographics <- subset(demographics, is.na(demographics$BLKGRP))
X <- demographics[, -c(1:5)]
for (i in seq_len(ncol(X))) {
  if (any(is.na(X[[i]]))) {
    cat(sprintf('%d: %f\n', i, mean(is.na(X[[i]])))) # elaborate on the message (i.e. this variable has this much missing values)
    X[[i]] <- 0    
  }
}
X <- log(1+X)

demog_pca <- prcomp(X, center = TRUE, retx = TRUE)
demog_variates <- cbind(subset(demographics, select = c('STATE', 'TRACT')), as.data.frame(demog_pca$x))
