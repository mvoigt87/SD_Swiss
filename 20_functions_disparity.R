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


#### apply ####

load("FLT_HMD.Rdata")

data <- fem.smooth[,c("Year","Age","dx","lx","ex")]

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

CV.FUN <- function(x){
  
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

IQR.FUN <- function(x){
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
# formula based on: www.demogr.mpg.de/papers/technicalreports/tr-2012-002.pdf
### alternative: (http://pages.stern.nyu.edu/~dbackus/BCH/demography/ZhangVaupel_ageseparating_DR_09.pdf)


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
  # computing e-dagger and Keyfitz' entropy based on Sholnikov/Andreev
  for(i in min(x$Age):max(x$Age)) {
    e.dagger <- sum(x$dx[x$Age==i]*1/2*(x$ex[x$Age==i]+x$ex[x$Age==i+1]))
  }
  return(e.dagger)
  
}

### apply
data1990 <- data %>% filter(Year==1990)

for (g in seq(10,110,0.1)) {
data1990$part.one[data1990$Age==g] <- data1990$dx[data1990$Age==g]*(data1990$ex[data1990$Age==g]+
                                                                      data1990$ex[data1990$Age==g+0.1])
}

for (g in min(data1990$Age):max(data1990$Age)) {
  e.dagr[g] <- (1/(2*lx[Age==g])) * cumsum()
}





data1990 %>% data1990 %>% mutate(part.one[Age] = dx[Age]*ex[Age+0.1]) %>% 
  mutate(part.two)

## Check diff - command
x <- cumsum(cumsum(1:10))
diff(x, lag = 2)



cumsum(data1990$dx[data1990$Age==10]:data1990$dx[max(data1990$Age)])


# single year
EDAG.FUN(data[data$Year == 1990,], smooth = TRUE, inter = seq(0,110,0.1))

for (i in 10:max(data.ch.f$Age)) {
  ED.test <- sum(data.ch.f$dx[data.ch.f$Year==1990 & data.ch.f$Age==i])*1/2*(data.ch.f$ex)
}



### time series for males and females

EDAG.CH.M <- by(data = data.ch.m, INDICES = data.ch.m$Year, FUN = IQR.FUN, smooth = TRUE, 
                inter = seq(0,110,0.1))
EDAG.CH.F <- by(data = data.ch.f, INDICES = data.ch.f$Year, FUN = IQR.FUN, smooth = TRUE, 
                inter = seq(0,110,0.1))


# creating a new data frame for easier plotting and handling the summarized values 
e.d.fem <- as.data.frame(unique(females$Year))
colnames(e.d.fem)[1] <- "Year"
e.d.fem <- e.d.fem %>% mutate(edagger=NA)
# for now with a for loop
for (i in 1876:2014) {
  e.d.fem$edagger[e.d.fem$Year==i] <- e.dagger.fun(subset(females,Year==i & Age<95))
}

## E-dagger plot (reasonable)
e.d.fem %>% ggplot(aes(x=Year, y=edagger))+ 
  geom_point() +
  scale_y_continuous(name = "e+ in years") +
  theme_bw()


# dplyr version (to be finished)

# EDAG <- females %>% group_by(Year) %>% summarise(edagger = sum(dxs[Age]*1/2*(ex[Age]+ex[Age+1])))


#######################
### Keyfitz entropy ###
#######################

e.d.fem$ex <- females$exs[females$Age==0]
e.d.fem$H <- e.d.fem$edagger/e.d.fem$ex

e.d.fem %>% ggplot(aes(x=Year, y=H))+ 
  geom_point() +
  scale_y_continuous(name = "Keyfitz' LT Entropy") +
  theme_bw()



##########################################
### From Modal Age at death to the SD+ ###
##########################################

## Calculate the mode over the the dx after age 5 (based on Canudas-Romo 2010)
## https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3000019/) - Appendix A (formula A2)

## Function to calculate the mode (after age 5) 
#  x is the age with the highest number of deaths in the life table at the time
MDA.fun <- function(lt) {
  x <- lt$Age[which.max(lt$dxs)]
  M <- x + ((lt$dxs[x]-lt$dxs[x-1])/(lt$dxs[x]-lt$dxs[x-1])+(lt$dxs[x]-lt$dxs[x+1]))
  return(M)
}


# test <- MDA.fun(subset(females,Year==1960 & Age>5))
# test 
# Seems to work fine with different years

# for each year - give me the age after 5 where the most deaths occur
# the not so elegant way for now:
mdsd.fem <- as.data.frame(unique(females$Year))
colnames(mdsd.fem)[1] <- "Year"
mdsd.fem <- mdsd.fem %>% mutate(M=NA)

for (j in 1876:2014) {
  mdsd.fem$M[mdsd.fem$Year==j] <- MDA.fun(subset(females, Year==j & Age>5))
}

summary(mdsd.fem)

### Plot

mdsd.fem %>% ggplot(aes(x=Year, y=M)) +
  geom_point() +
  scale_x_continuous(name="year") + 
  scale_y_continuous(name="modal age at death (after age 5)") +
  theme_bw()



## dplyr - version Calculate the modal age at death after age 5 by year
#  fem.test <- females %>% filter(Age>5) %>% group_by(Year) %>% 
#   dplyr::summarise(MDA = MDA.fun(lt))

############################
#### "source" functions ####
############################

