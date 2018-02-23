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

set.seed(17952)

# 0.3 set passworts for HMDHFDplus
name.m <- "mathias.voigt4@uni-rostock.de"
name.f <- "m.voigt87@gmx.de"
pw.m <- "1320270854"
pw.f <- "43700"


deaths <- readHMDweb(CNTRY = "CHE", item = "Deaths_1x1", username = name.m, password = pw.m)

deaths <- deaths[deaths$Age >= 10,]

expos  <- readHMDweb(CNTRY = "CHE", item = "Exposures_1x1", username = name.m, password = pw.m)

expos <- expos[expos$Age >= 10,]

## putting the deaths and exposures in a matrix format (age x year)
D <- do.call(cbind,tapply(X = deaths$Female, INDEX = deaths$Year, FUN = identity))

E <- do.call(cbind,tapply(X = expos$Female, INDEX = expos$Year, FUN = identity))

#### interpolate ####

y <- unique(deaths$Year)
x <- unique(deaths$Age)
m <- length(x)
n <- length(y)

# some E are equal to zer0, -> weights are necessary  (Wheight matrix)
W <- matrix(1, m, n)
W[E==0] <- 0

# subjective lambdas (previously optmize for swiss males by GC) ---> same for females??

fit <- Mort2Dsmooth(x = x, y = y,
                    Z = D, offset=log(E),
                    W=W,
                    method=3, lambdas=c(3.2, 32))
fit$lambdas
# plots
plot(fit)
plot(fit$logmortality, log(D/E))

## compute the density from the fitted hazard at finer grid
## finer-grid age - decimal ages (0.1)
delta <- 0.1
xs <- seq(min(x), max(x), delta)
ms <- length(xs)

## new basis over new ages
xl <- min(x)
xr <- max(x)
xmax <- xr + 0.01 * (xr - xl)
xmin <- xl - 0.01 * (xr - xl)
# generate an equally-spaced B-Splines basis over the abscissa (applying Mort1Dsmooth)
Bxs <- MortSmooth_bbase(xs, xmin, xmax, fit$ndx[1], fit$deg[1])

## over years are the same
By <- fit$By
## fitted coefficients
betas <- fit$coef
## log-mortality (linear predictor) over new ages and years (output: log hazard)
ln.h <- MortSmooth_BcoefB(Bxs, By, betas)
## hazard
h <- exp(ln.h)
## cumulative hazard using cumsum
H <- matrix(0, ms, n)
for(i in 1:n){
  H[,i] <- cumsum(h[,i]*delta)
}

# ## fitted hazard functions and its cumulative
# h0 <- exp(fit$logmortality)
# h <- matrix(NA, ms, n)
# H <- matrix(0, ms, n)
# for(i in 1:n){
#   fun <- splinefun(x, h0[,i])
#   h[,i] <- fun(xs)
#   for(j in 2:ms){
#     H[j,i] <- integrate(fun, xs[1], xs[j])$value
#   }
# }
## fitted survival functions
S <- apply(H, 2, function(x){exp(-x)})
## fitted density functions
f <- h * S

image(t(f))


#### building life tables from the fx ####

dim(f)
# radix of the table by year
colSums(f)

# survival function by year (beautiful ;) )
matplot(S, type = "l", lty = 1, col = heat.colors(n = ncol(S)))

# rough estimate of the life expectancy at birth
colSums(S) / 10
plot(y, colSums(S) / 10, type = "b", las = 1, main = "e0")


# life table
# reminder of the main functions
# library(MortHump)
# LT

# build life table

##  help variables
# 1001 age groups
N <- length(xs)
# the 0.1 differences
Widths <- rep(delta, length(xs))



## FEMALES
## large data frame which contains values in HMD format
fem.smooth <- as.data.frame(rep(xs,n))
colnames(fem.smooth) <- "xs"  
fem.smooth <- fem.smooth %>% 
  mutate(Year=rep(min(y):max(y), times=1, each=1001)) %>% 
  # ax values (may be to be changed for the highest age groups)
  mutate(ax = rep(Widths / 2, times=n))
  # ... step in between

# ------------------------------------------------------------
# obtain the mx values from the smoothed hazard function
dim(h)
      ## get the hx in the right format
      h.new <- as.data.frame(h)
      h.new <- data.frame(mx=unlist(h.new, use.names = FALSE))
# and add them as column to the femsmooth data frame
# ------------------------------------------------------------  

# continuing building the life table      
fem.smooth <- fem.smooth %>% bind_cols(h.new) %>% 
    # qx
    mutate(qx = (Widths * mx) / (1 + (Widths - ax) * mx)) 
# ------------------------------------------------------------
 ## a little trick which would not work with a data frame to make the last qx=1
 qx <- matrix(fem.smooth$qx)
 qx[1:(0+1001)==(0+1001)] <- 1
# ------------------------------------------------------------      
 fem.smooth <- fem.smooth %>% select(-qx) %>%  bind_cols(as.data.frame(qx))
 colnames(fem.smooth)[5] <- "qx"
 fem.smooth <- fem.smooth %>% mutate(qx = ifelse(qx>1,1,qx)) %>% 
 ## add the px
 mutate(px = 1 - qx)
 ## lx beginning at the estimated radix (sum f[,1])
# ------------------------------------------------------------   
    ## matrix operations: sum over the columns of the estimated f-values to obtain
    ## the base/radix for the life table
    radix.mat <- as.data.frame(matrix(data=colSums (f, na.rm = FALSE, dims = 1), nrow = 1)) %>% 
    ## now making filling dummie values in between to make it the same length as the data frame
    bind_rows(as.data.frame(matrix(data = 0,nrow = 1000, ncol = 139))) 
    ## stack them in order and delete the extra variable
    radix.mat <- stack(radix.mat) %>% select(-ind)
# ------------------------------------------------------------   
    
## now add them to the large life table
 fem.smooth <- fem.smooth %>% bind_cols(as.data.frame(radix.mat))
 colnames(fem.smooth)[7] <- "lx"
## use the dplyr group_by command to calculate the rest of the lx from the px
 fem.smooth <- fem.smooth %>% group_by(Year) %>% mutate(lx = c(lx[1],lx[1] * cumprod(px))[1:N]) %>% 
## dx values from the lx (alternatively from the smoothing algorithm)
 group_by(Year) %>%  mutate(dx = c(-diff(lx),lx[N])) %>% 
## Create the Lx from the lx and the dx
 group_by(Year) %>% mutate(Lx = c(Widths[1:(N - 1)] * lx[2:N] + ax[1:(N - 1)] * dx[1:(N - 1)], lx[N] * ax[N])) %>% 
 ## account for infinite Lx and NA
 mutate(Lx = ifelse(is.infinite(Lx),1,Lx)) %>% mutate(Lx = ifelse(is.na(Lx),0,Lx)) %>% 
## Calculate the Tx from the Lx
 group_by(Year) %>% mutate(Tx = rev(cumsum(rev(Lx)))) %>% 
## Finally obtain the life expectancy from the Tx and lx
 group_by(Year) %>% mutate(ex = Tx / lx)
      

 ## test plot
 
 fem.smooth %>% filter(xs==10) %>% ggplot(aes(x=Year,y=ex)) +
                    geom_line() +
                    theme_bw()
 

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


###################################################################################################
###################################################################################################



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
