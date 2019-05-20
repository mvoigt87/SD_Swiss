
# Life tables based on the SNC
# 2000 - 2014
# Adrien Remund, Feb 2018

#### load & format data ####

# set your own working drive
setwd("D:/switchdrive/SNC 2.0/Data/R")

# load SNC data from 2000 census onwards
load("SNC_00_14.RData")

# 7280246 individuals and 178 variables
dim(snc)

# status (dead/alive) at the end of observation
snc$fail  <-  as.numeric(!is.na(snc$dod))

# reduce to necessary variables
# dstop_imputed = end of observation (including imputed deaths = bad match)
snc <- snc[,c("snc_ge_id","fail","sex","dob","dstop_imputed")]

# start of observation = 2000 census
snc$entry <- as.Date("2000-12-05")

# rename variables
names(snc) <- c("id","fail","sex","birth","exit","entry")

#### compute exposures and death counts by sex, age and period ####

library(Epi)

# create a list to receive the results
out <- list(males = NULL, females = NULL)

# loop over sex
# this will take a LONG time due to the size of the dataset!
for(i in unique(snc$sex)){
  
  # select the subset by sex
  data <- snc[snc$sex == i,]
  
  # format the variables to be used by the Epi package
  data$bt <- cal.yr(data$birth, format="%Y-%m-%d")
  data$en <- cal.yr(data$entry, format="%Y-%m-%d")
  data$ex <- cal.yr(data$exit , format="%Y-%m-%d")
  omega <- ceiling(max(data$ex - data$bt))
  
  # generate a lexis object
  lex <- Lexis( entry = list(per = en),
                exit = list(per = ex, age = ex-bt),
                exit.status = fail,
                id = id,
                data = data)
  
  # split along time axis
  # we want to compare 2000-2004 vs 2010-2014, so we break in 2005 and 2010
  # to get yearly estimates, we'd replace by "breaks = 2000:2014"
  x2 <- splitLexis(lex, breaks = c(2000,2005,2010), time.scale = "per")
  
  # split along age axis
  x2 <- splitLexis(x2, breaks = seq(0,omega,5), time.scale = "age")
  
  # tabulate death counts
  dx <- tapply(status(x2, "exit") == 1, list(timeBand(x2, "age", "left"), timeBand(x2, "per", "left")), sum)
  dx[is.na(dx)] <- 0
  
  # tabulate exposures
  nx <- tapply(dur(x2), list(timeBand(x2,"age","left"), timeBand(x2,"per","left")), sum)
  nx[is.na(nx)] <- 0
  
  # compute rates
  mx <- dx / nx
  
  # store results in list
  out[[which(unique(snc$sex) == i)]] <- list(dx = dx, nx = nx, mx = mx)
  
}

#### compute life tables ####

# look at the age-specific death rates
matplot(log(out$males$mx), type = "s", lty = 1, las = 1, col = heat.colors(ncol(out$males$mx)))

# extract one year
mx <- out$males$mx[,1]

# build table
radix <- 1
x <- as.numeric(rownames(out$males$mx))
N <- length(x)
Widths <- rep(5,N)
ax <- rep(Widths / 2) # with 10th of years, we can try this, but this is debatable
qx <- (Widths * mx) / (1 + (Widths - ax) * mx)
qx[N] <- 1 # last qx must be 1 (everyone dies eventually)
qx[qx > 1] <- 1 # can't have qx > 1 at any age
px <- 1 - qx
lx <- c(radix, radix * cumprod(px))[1:N]
dx <- c(-diff(lx),lx[N])
Lx <- c(Widths[1:(N - 1)] * lx[2:N] + ax[1:(N - 1)] * dx[1:(N - 1)], lx[N] * ax[N])
Lx[is.infinite(Lx)] <- 1 # in case we have a hickup
Lx[is.na(Lx)] <- 0 # same
Tx <- rev(cumsum(rev(Lx)))
ex <- Tx / lx
ex[N] <- ifelse(mx[N] == 0, ax[N], {1 / mx[N]})
ex

#### save all outputs ####

save(out, file = "D:/switchdrive/SNC 2.0/TOMCOD/Output/2000vs2010.RData")