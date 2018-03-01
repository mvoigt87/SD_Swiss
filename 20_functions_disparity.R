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

### load smoothed life tables

load("FLT_HMD.Rdata")

load("MLT_HMD.RData")

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


#### apply ####

data <- readHMDweb(CNTRY = "CHE", item = "mltper_1x1", username = name.m, password = pw.m)

data <- data[,c("Year","Age","dx")]

sdfun(data[data$Year == 1988,], smooth = TRUE, plot = TRUE, trun = "mode", inter = seq(0,110,0.1))

### time series

sd <- by(data = data, INDICES = data$Year, FUN = sdfun, trun = "mode", smooth = TRUE, inter = seq(0,110,0.1))

plot(unique(data$Year), sd, type = "l", xlab = "", ylab = "sd", las = 1)


#### comparison ####

par(mfrow = c(1,2))

for(trun in c(0,"mode")){
  
  plot(1, 1, type = "n", ylim = c(ifelse(trun == 0, 10, 3),ifelse(trun == 0, 35, 6)), xlim = c(1960, 2020), xlab = "", ylab = "sd", las = 1, main = paste("truncated at",trun))
  
  for(ctry in c("CHE","NLD")){
    
    cat("\n",ctry)
    
    for(sex in c("m","f")){
      
      cat("\n","-",sex)
      
      data <- readHMDweb(CNTRY = ctry, item = paste(sex,"ltper_1x1",sep = ""), username = name.m, password = pw.m)
      
      data <- data[,c("Year","Age","dx")]
      
      sd <- by(data = data, INDICES = data$Year, FUN = sdfun, trun = trun, smooth = TRUE, inter = c(seq(0,110,0.1)))
      
      col <- ifelse(sex == "m", "blue", "red")
      
      lty <- ifelse(ctry == "CHE", 1, 2)
      
      lines(unique(data$Year), sd, col = col, lty = lty)
      
    }
    
  }
  
  legend("topright", legend = c("CH females", "CH males", "NL females", "NL males"), lty = c(1,1,2,2), col = c("red","blue","red","blue"))
  
}

par(mfrow = c(1,1))

dev.off()

#####

###########
### IQR ###
###########

### Attempt to write a similar efficient function compared to the sd.fun
## Input: life table-data frame (more specific, dx values, x, and years)

IQR.FUN <- function(x, smooth = FALSE, inter = FALSE){
  
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
  ### Now the IQR function:
  x$rel.DX.CUM <- cumsum(x$dx)/max(cumsum(x$dx))
  Q1 <- x$Age[min(which(x$rel.DX.CUM>0.25))]
  Q3 <- x$Age[min(which(x$rel.DX.CUM>0.75))]
  IQR <- Q1 - Q3
  return(IQR)
}


#### apply IQR ####

# male
data.ch.m <- readHMDweb(CNTRY = "CHE", item = "mltper_1x1", username = name.m, password = pw.m)
data.ch.m <- data.ch.m[,c("Year","Age","dx")]

# female
data.ch.f <- readHMDweb(CNTRY = "CHE", item = "fltper_1x1", username = name.m, password = pw.m)
data.ch.f <- data.ch.f[,c("Year","Age","dx")]

# single year
IQR.FUN(data.ch.m[data.ch.m$Year == 1988,], smooth = TRUE, inter = seq(0,110,0.1))

### time series for males and females

IQR.CH.M <- by(data = data.ch.m, INDICES = data.ch.m$Year, FUN = IQR.FUN, smooth = TRUE, 
               inter = seq(0,110,0.1))
IQR.CH.F <- by(data = data.ch.f, INDICES = data.ch.f$Year, FUN = IQR.FUN, smooth = TRUE, 
               inter = seq(0,110,0.1))

# plot(unique(data.ch.m$Year), IQR.CH.M, type = "l", xlab = "", ylab = "IQR", las = 1)
# lines(unique(data.ch.m$Year),IQR.CH.F, col="red")


### Plot IQR (ggplot)

IQR.M <- as.data.frame(cbind(IQR.CH.M,unique(data.ch.m$Year))) %>% mutate(sex="male")
colnames(IQR.M) <- c("IQR","Year","sex")
IQR.F <- as.data.frame(cbind(IQR.CH.F,unique(data.ch.f$Year))) %>% mutate(sex="female")
colnames(IQR.F) <- c("IQR","Year","sex")

IQR <- bind_rows(IQR.F,IQR.M)

IQR %>% ggplot(aes(x=Year,y=IQR,color=sex)) +
  geom_line() +
  scale_x_continuous(name="") +
  scale_y_continuous(name="IQR in years") +
  scale_colour_manual(values = c("orange", "darkgrey"), name="") +
  theme_bw()


#######################################################################################################
#######################################################################################################
#######################################################################################################


#####################################
### CV - coefficient of variation ###
#####################################

# source: www.demogr.mpg.de/papers/technicalreports/tr-2012-002.pdf

### 1. First step

## assign new dx values (smoothed)

dx.FUN <- function(x, smooth = FALSE, inter = FALSE){
  
  ### If smooth dx' are provided
  if(smooth == TRUE){
    
    if(is.numeric(inter)){
      
      dxs <- predict(smooth.spline(x = x$Age, y = x$dx, control.spar = list(low = 0.3, high = 1)), x = inter)$y
      
      x <- data.frame(Year = rep(x$Year,length(inter)), Age = inter)
      
      x$dx <- dx #* c(diff(inter),1)
      
    }else{
      
      dxs <- predict(smooth.spline(x = x$Age, y = x$dx, control.spar = list(low = 0.3, high = 1)))$y
    }
  }
  
  return(dxs)
  
}

### 2. Second step - build the life tables on the "new" dx values
## !!! For the indicators which require information on life expectancy

CV.FUN <- function(x){
  
  # where x is the life table with the smoothed dx
  
  # Obtain the lx from the dxs
  x$lx <- 1:length(x$Age)
  x$lx[1] <- 1
  for (i in 2:length(x$Age)) {
    x$lx[i] <- x$lx[i - 1] - x$dxs[i - 1]
  }
  # Lx values
  x$Lx <- x$lx - (0.5*x$dx)
  # Tx values
  x$Tx <- rev(cumsum(rev(x$Lx)))
  # ex values
  x$ex[x$Age] <- x$Tx[x$Age] / x$Lx[x$Age]
  
  # build from the SD and the life expectancy
  
  mean <- sum(x$Age * x$dx) / sum(x$dxs)
  
  sdev <- sqrt(sum((x$Age - mean)^2 * x$dxs) / sum(x$dxs))
  
  # coefficient of variance - CV
  cv <- (sdev)/(x$Age + x$ex)
  
  return(CV)
}



#### apply ####

# male
data.ch.m <- readHMDweb(CNTRY = "CHE", item = "mltper_1x1", username = name.m, password = pw.m)
data.ch.m <- data.ch.m[,c("Year","Age","dx")]

# female
data.ch.f <- readHMDweb(CNTRY = "CHE", item = "fltper_1x1", username = name.m, password = pw.m)
data.ch.f <- data.ch.f[,c("Year","Age","dx")]

### single year

dxs.m.1990 <- dx.FUN(data.ch.m[data.ch.m$Year == 1990,], smooth = TRUE, inter = seq(0,110,0.1))

### time series for males and females

dxs.m <- by(data = data.ch.m, INDICES = data.ch.m$Year, FUN = dx.FUN, smooth = TRUE, 
            inter = seq(0,110,0.1))
data.ch.f$dxs <- by(data = data.ch.f, INDICES = data.ch.f$Year, FUN = dx.FUN, smooth = TRUE, 
                    inter = seq(0,110,0.1))

### Apply the CV function

CV.CH.mal <- by(data = data.ch.m, INDICES = data.ch.m$Year, FUN = CV.FUN, smooth = TRUE, 
                inter = seq(0,110,0.1))


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
    e.dagger <- sum(x$dxs[x$Age==i]*1/2*(x$ex[x$Age==i]+x$ex[x$Age==i+1]))
    # H <- e.dagger/x$ex[x$Age==0]
  }
  return(e.dagger)
  
}

### apply

# single year
EDAG.FUN(data.ch.m[data.ch.m$Year == 1990,], smooth = TRUE, inter = seq(0,110,0.1))

for (i in 0:max(data.ch.f$Age)) {
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

