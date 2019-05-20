#####
##### 2.0 - Functions for Estimating life span variability indicators
#####
##### Exploring the relationship between the modal age at deaths and the standard deviation above it #####
##### ---------------------------------------------------------------------------------------------- #####

## Following Cheung et al (2009), the goal of this project is to investigate the characteristics of the 
## longevity development in Switzerland. We assume that the increase in life expectancy which went along with
## the compression around the modal age of death (decreasing standard deviation) is increasingly driven by the
## shift of the modal age at death to higher ages. Looking at the variability measures could proof that hypothesis


# 0.1 load packages
library(HMDHFDplus)
library(MortalitySmooth)
library(MortHump)
library(tidyverse)
# library(plyr)

#### set working directory ####
if(Sys.info()["nodename"] == "CD-5VV9QK2"){
  
  setwd("D:/switchdrive/SNC 2.0/TOMCOD")
  
}


###################################################################
#### function to obtain the standard deviation, mean ages etc. ####

sdfun <- function(x, smooth = FALSE, plot = FALSE, trun = 0, inter = FALSE){
  
  if(plot == TRUE){
    plot(x$Age, x$dx, las = 1, xlab = "age", ylab = "dx")
  }
  
  if(is.numeric(trun)){
    
    x <- x[x$Age >= trun,]
    
    if(plot == TRUE){
      
      abline(v = trun)
      
    }
    
  }
  
  if(smooth == TRUE){
    
    if(is.numeric(inter)){
      
      dx <- predict(smooth.spline(x = x$Age, y = x$dx, control.spar = list(low = 0.3, high = 1)), x = inter)$y
      
      x <- data.frame(Year = rep(x$Year,length(inter)), Age = inter)
      
      x$dx <- dx #* c(diff(inter),1)
      
    }else{
      
      x$dx <- predict(smooth.spline(x = x$Age, y = x$dx, control.spar = list(low = 0.3, high = 1)))$y
      
    }
    
    if(plot == TRUE){
      lines(x$Age, x$dx)
    }
    
  }
  
  if(trun == "mode"){
    
    ranks <- order(x$dx, decreasing = TRUE)
    
    i <- 1
    idx <- ranks[i]
    mode <- x$Age[idx]
    while(mode < 5){
      i <- i + 1
      idx <- ranks[i]
      mode <- x$Age[idx]
    }
    
    x <- x[x$Age >= mode,]
    
    if(plot == TRUE){
      
      abline(v = mode)
      
    }
    
  }
  
  mean <- sum(x$Age * x$dx) / sum(x$dx)
  
  sdev <- sqrt(sum((x$Age - mean)^2 * x$dx) / sum(x$dx))
  
  return(sdev)
  
}

### save function for standard deviation around the mode

dump("sdfun", file = "sdfun.R")


#### apply - test ####

load("FLT_HMD.Rdata")

data <- fem.smooth[,c("Year","Age","dx","ax","lx","ex")]

## Standard deviation above the mode
sdfun(data[data$Year == 1900,], smooth = TRUE, plot = TRUE, trun = "mode", inter = seq(0,110,0.1))

### time series

sd <- by(data = data, INDICES = data$Year, FUN = sdfun, trun = "mode", smooth = TRUE, inter = seq(0,110,0.1))

plot(unique(data$Year), sd, type = "l", xlab = "", ylab = "sd", las = 1)

##################################################################################
##################################################################################
##################################################################################

#####################################
### CV - coefficient of variation ###
#####################################

# source Sholnikov/Andreev 2010: http://www.demogr.mpg.de/papers/technicalreports/tr-2010-001.pdf

### Absolute measure for dispersion

## Assuming the dx are smoothed in a previous step

CV.FUN <- function(x, smooth=FALSE){
  
  ### Smoothing the values / or not
  if(smooth == TRUE){
    
    if(is.numeric(inter)){
      
      dx <- predict(smooth.spline(x = x$Age, y = x$dx, control.spar = list(low = 0.3, high = 1)), x = inter)$y
      
      x <- data.frame(Year = rep(x$Year,length(inter)), Age = inter)
      
      x$dx <- dx #* c(diff(inter),1)
      
    }else{
      
      x$dx <- predict(smooth.spline(x = x$Age, y = x$dx, control.spar = list(low = 0.3, high = 1)))$y
      
    }}
  
  mean <- sum(x$Age * x$dx) / sum(x$dx)
  
  sdev <- sqrt(sum((x$Age - mean)^2 * x$dx) / sum(x$dx))
  
  # coefficient of variance - CV
  cv <- ((sdev)/mean) * 100
  
  return(cv)
  
}


#### save as a global function ####

dump("CV.FUN", file="CV_FUN.R")


  # ## test
  # data <- mal.smooth %>% select(dx,Age,Year,ex)
  # data$Year <- as.numeric(data$Year)
  # 
  # CV.1990 <- CV.FUN(data[data$Year == 1990,])


##################################################################################
##################################################################################
##################################################################################

###########
### IQR ###
###########

### Attempt to write a similar efficient function compared to the sd.fun
## Input: life table-data frame (more specific, dx values, x, and years)

IQR.FUN <- function(x, smooth=FALSE){
  
  ### Smoothing the values / or not
  if(smooth == TRUE){
    
    if(is.numeric(inter)){
      
      dx <- predict(smooth.spline(x = x$Age, y = x$dx, control.spar = list(low = 0.3, high = 1)), x = inter)$y
      
      x <- data.frame(Year = rep(x$Year,length(inter)), Age = inter)
      
      x$dx <- dx #* c(diff(inter),1)
      
    }else{
      
      x$dx <- predict(smooth.spline(x = x$Age, y = x$dx, control.spar = list(low = 0.3, high = 1)))$y
      
    }}
  
    ### IQR function:
  x$rel.DX.CUM <- cumsum(x$dx)/max(cumsum(x$dx))
  Q1 <- x$Age[min(which(x$rel.DX.CUM>0.25))]
  Q3 <- x$Age[min(which(x$rel.DX.CUM>0.75))]
  IQR <- Q3 - Q1
  return(IQR)
}

### Save the function in the R Profile

dump("IQR.FUN", file="IQR_FUN.R")

  ### test
  # data <- mal.smooth %>% select(dx,Age,Year)
  # data$Year <- as.numeric(data$Year)
  # 
  # IQR1990 <- IQR.FUN(subset(data,Year==1990))



#######################################################################################################
#######################################################################################################
#######################################################################################################

################
### e-Dagger ###
################

### ---------------------------------------------------------
# formula based on: Shkolnikov et al. 2011 Losses in life expectancy in the USA
# https://link.springer.com/article/10.1007%2Fs13524-011-0015-6
### alternative continuous version: (http://pages.stern.nyu.edu/~dbackus/BCH/demography/ZhangVaupel_ageseparating_DR_09.pdf)

EDAG.FUN <- function(x, smooth = FALSE, inter = FALSE){
  
  ### If smooth dx' are provided
  if(smooth == TRUE){
    
    if(is.numeric(inter)){
      
      dx <- predict(smooth.spline(x = x$Age, y = x$dx, control.spar = list(low = 0.3, high = 1)), x = inter)$y
      
      x <- data.frame(Year = rep(x$Year,length(inter)), Age = inter)
      
      x$dx <- dx #* c(diff(inter),1)
    }else{
      
      x$dx <- predict(smooth.spline(x = x$Age, y = x$dx, control.spar = list(low = 0.3, high = 1)))$y
      
    }
  }
  # computing e-dagger 
  
  y <- last(x$Age)
  part.one <- sum(x$dx[-y]*x$ex[-1])
  part.two <- 1-(sum(x$dx[-y]*x$ax[-y]))
  edagger <- part.one + part.two
  #H <- edagger/x$ex[1]
  return(edagger)
  
}

### save e+ - function
dump("EDAG.FUN", file="EDAG_FUN.R")


### apply test
data1990 <- data %>% filter(Year==2005)
edag1990 <- EDAG.FUN(data1990, smooth = TRUE)


## e+(10) is much higher than it would be due to missing infant and child mortality
## might have to re-run it with life tables which start at age 0


##################################################################################
##################################################################################
##################################################################################


#######################################
### Modal Age at death above age 10 ###
#######################################

## Calculate the mode over the the dx after age 5 - here age 10 by default
## based on Canudas-Romo 2010: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3000019/) - Appendix A (formula A2)

## and Kannisto (2001) - Mode and Dispersion of the Length of Life
## http://www.jstor.org/stable/3030264

## Function to calculate the mode (after age 5) 
#  x is the age with the highest number of deaths in the life table at the time
MDA.fun <- function(lt) {
  x <- lt$Age[which.max(lt$dx)]
  M <- x + ((lt$dx[x]-lt$dx[x-1])/(lt$dx[x]-lt$dx[x-1])+(lt$dx[x]-lt$dx[x+1]))
  return(M)
}

dump("MDA.fun","MA5_FUN.R")

# test
test <- MDA.fun(subset(fem.smooth,Year==1960))



#### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### !!!!!!!!!!!!!!!!!!!!!!! Every below this is a Test range



## attempt to copy Cheung et al. 2009
SD.plus.FUN <- function(x, smooth = FALSE, inter = FALSE){
  
  ### If smooth dx' are provided
  if(smooth == TRUE){
    
    if(is.numeric(inter)){
      
      dx <- predict(smooth.spline(x = x$Age, y = x$dx, control.spar = list(low = 0.3, high = 1)), x = inter)$y
      
      x <- data.frame(Year = rep(x$Year,length(inter)), Age = inter)
      
      x$dx <- dx #* c(diff(inter),1)
    }else{
      
      x$dx <- predict(smooth.spline(x = x$Age, y = x$dx, control.spar = list(low = 0.3, high = 1)))$y
      
    }
  }
# function by Cheung et al.
  a <- x$Age[which.max(x$dx)]
  M <- a + ((x$dx[a]-x$dx[a-1])/(x$dx[a]-x$dx[a-1])+(x$dx[a]-x$dx[a+1]))
  b <- seq(ceiling(M),last(x$Age),0.1)
  part.one <- c(rep(0,length(b)))
  
  ## for loop (for now)
  for (k in b) {
    part.one <- (sum(b[k]-M))^2/length(b) 
    }
  SD.plus <- sqrt(part.one)
  ### !!! For comparison with other estimates - IF decimal ages are used
  SD.plus <- SD.plus*10
  return(SD.plus)
}

dump("SD.plus.FUN","Cheung_FUN.R")

# test
test.dos <- SD.plus.FUN(subset(fem.smooth,Year==1995), smooth = TRUE)


td <- subset(mal.smooth,Year==1950)

a <- td$Age[which.max(td$dx)]
M <- a + ((td$dx[a]-td$dx[a-1])/(td$dx[a]-td$dx[a-1])+(td$dx[a]-td$dx[a+1]))
b <- seq(ceiling(M),last(td$Age),0.1)
part.uno <- c(rep(0,length(b)))

for (k in b) {
  part.uno <- ((sum(b[k]-M))^2)/(length(b))            # length of the age interval above (not very elegant)
}


SD.plus <- sqrt(part.uno)*10

  
#########################################################################################################
#########################################################################################################
#########################################################################################################
#########################################################################################################
#########################################################################################################  

########################
### GINI Coefficient ###
########################

## based on MPIDR working paper by Shkolnikov/Andreev (2010)
## http://www.demogr.mpg.de/papers/technicalreports/tr-2010-001.pdf
## or alternatively: Shkolnikov et al 2003: https://www.demographic-research.org/volumes/vol8/11/8-11.pdf


GINI.FUN <- function(x){
  # F(x) - share of the life table population
  # Omega(x) - cumulative share
  y <- last(x$Age)
  
  x$Fx <- c(rep(0,length(x$Age[-1])),1)
  for (i in 2:length(x$Age[-1])){
    x$Fx[i] <- 1-(x$lx[i+1]/x$lx[1])
  }
  
  x$Omega <- c(rep(0,length(x$Age[-1])),1)
  for (m in 2:length(x$Age[-1])){
    x$Omega[m] <- sum(head(x$dx,m)*x$Age[m])/sum(x$dx*y)
  }
  x$GINI = c(rep(0,length(x$Age)))
  for (q in 1:length(x$Age[-1])){
    x$GINI[q] <- 1 - (sum((x$Fx[q+1]-x$Fx[q])*(x$Omega[q+1]+x$Omega[q])))
  }
  return(GINI)
}



### apply test
data1990 <- data %>% filter(Year==1880)
gini1990 <- GINI.FUN(data1990)

data1990$Fx <- c(rep(0,length(data1990$Age[-1])),1)
for (i in 2:length(data1990$Age[-1])){
  data1990$Fx[i] <- 1-(data1990$lx[i+1]/data1990$lx[1])
}

data1990$Omega <- c(rep(0,length(data1990$Age[-1])),1)
for (m in 2:length(data1990$Age[-1])){
  data1990$Omega[m] <- sum(head(data1990$dx,m)*data1990$Age[m])/sum(data1990$dx*y)
}

### Lorenz Curve
plot(x=data1990$Fx,y=data1990$Omega, type = "l")
abline(a = 0, b = 1, col = 2)


# Discrete GINI

part.uno <- 1-(1/(data1990$ex * data1990$lx))
part.dos <- c(rep(1,length(data1990$Age)))
for (h in 2:length(data1990$Age[-1])) {
  part.dos[h] <- sum(((data1990$lx[h-1])^2 + data1990$ax[h] * ((data1990$lx[h])^2 - data1990$lx[h-1])^2))
}

GINI1990 <- part.uno * part.dos

# continuous GINI
data1990$GINI = c(rep(1,length(data1990$Age)))
for (q in 2:length(data1990$Age[-1])){
  data1990$GINI[q] <- 1 - (sum((data1990$Fx[q+1]-data1990$Fx[q]))*(data1990$Omega[q+1]+data1990$Omega[q]))
}


