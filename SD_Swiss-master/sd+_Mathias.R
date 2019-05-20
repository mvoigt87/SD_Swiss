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
library(plyr)
# library(LifeTables)

# 0.2 working directory
setwd("C:/Users/Mathias/Documents/LongPop_Madrid/Secondment/Geneva/Compression of Mortality/Data and Code")

# 0.3 set passworts for HMDHFDplus
name.m <- "mathias.voigt4@uni-rostock.de"
name.f <- "m.voigt87@gmx.de"
pw.m <- "1320270854"
pw.f <- "43700"


females <- readHMDweb(CNTRY = "CHE", item = "fltper_1x1", username = name.m, password = pw.m)

# Smooth the mortality surface with P-Splines
fit <- Mort2Dsmooth(x = unique(females$Age), y = unique(females$Year), Z = matrix(females$mx, nrow = 111), 
                    method = 3, lambdas = c(.01,.01))
# plots
plot(fit)

# use the predict function to interpolate
females$mxs <- as.vector(predict(fit, type = "response"))

lt <- by(data = females$mxs, INDICES = females$Year, FUN = function(x){LT(Mx = x, mxsmooth = FALSE, 
                                                                          axsmooth = FALSE)})
## Compare the estimated values with observed
females$dxs <- unlist(lapply(lt, function(x){x$dx}))
females$exs <- unlist(lapply(lt, function(x){x$ex}))

plot(females$ex, females$exs)
head(females)

#femalesplus <- data.frame(Year = rep(unique(females$Year), each = 1100), x = rep(seq(0,109.9,0.1),length(unique(females$Year))))
#femalesplus$dxs <- predict(fit, newdata = data.frame(y = femalesplus$Year, x = femalesplus$x), type = "response")



                      ######## Interpolation by GC/Adrien ##########
###################################################################################################
###################################################################################################
###################################################################################################
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

## fitted hazard functions and its cumulative (H)

# baseline hazard
h0 <- exp(fit$logmortality)
# prepare matrizes
h <- matrix(NA, ms, n)
H <- matrix(0, ms, n)

# Apply a cubic spline function to the h0 [at time i] and then over the decimal ages
for(i in 1:n){
  fun <- splinefun(x, h0[,i])
  h[,i] <- fun(xs)
  for(j in 2:ms){
    H[j,i] <- integrate(fun, xs[1], xs[j])$value
  }
}
## fitted survival functions
S <- apply(H, 2, function(x){exp(-x)})
## fitted density functions
f <- h * S

image(t(f))

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

sdfun(data[data$Year == 1942,], smooth = TRUE, plot = TRUE, trun = "mode", inter = seq(0,110,0.1))

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
###################################################################################################





### ---------
### For Males
### ---------

males <- readHMDweb(CNTRY = "CHE", item = "mltper_1x1", username = name.m, password = pw.m)
# Smooth the mortality surface with P-Splines
fit <- Mort2Dsmooth(x = unique(males$Age), y = unique(males$Year), Z = matrix(males$mx, nrow = 111), 
                    method = 3, lambdas = c(.01,.01))
plot(fit) # Spanish flu visible
# extract smoothed life table mx and calculate the rest of the indicators
males$mxs <- as.vector(predict(fit, type = "response"))
lt <- by(data = males$mxs, INDICES = males$Year, FUN = function(x){LT(Mx = x, mxsmooth = FALSE, axsmooth = FALSE)})
males$dxs <- unlist(lapply(lt, function(x){x$dx}))
males$exs <- unlist(lapply(lt, function(x){x$ex}))
plot(males$ex, males$exs)
head(males)

## test if outliers were generated by the Spanish flu (year 1918)

males.without <- males %>% filter(Year!=1918)
par(mfrow=c(1,2))
plot(males$ex, males$exs)
plot(males.without$ex, males.without$exs)
par(mfrow=c(1,1))

# looks like

# functions to compute e-dagger and sd+ (check this one though!)

# e-dagger - the average life expectancy losses due to death. 

# edagger <- function(lt){
#   
#   ed <- sum(lt$dxs / sum(lt$dxs) * lt$exs)
#   
#   return(ed)
#   
# }

 # sdplus <- function(lt){
 #  
 #  dx <- lt$dxs / sum(lt$dxs)
 #  x <- lt$Age
 #  m <- which.max(dx[x > 50]) + 51
 #  dx <- dx[m:length(dx)]
 #  x  <-  x[m:length(x)]
 #  mean <- sum(x * dx) / sum(dx)
 #  s <- sqrt(sum((x - m)^2 * dx) / sum(dx))
 # 
 #  return(s)
 # 
 # }
 # 
 # sdplus(females)
 
# 
# # apply both functions
# png("EPC/edagger.png", width = 30, height = 15, units = "cm", res = 600)
# par(mfrow = c(1,2), mar = c(3,4,3,1) + 0.1)
# plot(unique(females$Year), as.vector(by(data = females, INDICES = females$Year, FUN = edagger)), type = "l", xlab = "", ylab = expression(e^dagger), las = 1, main = "Swiss females", ylim = c(5,25))
# plot(unique(males$Year), as.vector(by(data = males, INDICES = males$Year, FUN = edagger)), type = "l", xlab = "", ylab = expression(e^dagger), las = 1, main = "Swiss males", ylim = c(5, 25))
# par(mfrow = c(1,1), mar = c(5,4,4,1) + 0.1)
# dev.off()
# 
# plot(unique(females$Year), as.vector(by(data = females, INDICES = females$Year, FUN = sdplus)), type = "l")




### ----------------- Experimenting !!!

#######################################################
### Standard deviation (Tuljarpurka & Edwards 2005) ###
#######################################################

# dx matrix - females
dx.fem <- tapply(X=females$dxs,
                 INDEX=list(Age=females$Age,
                            Year=females$Year),
                 FUN=sum)

# ex matrix - females
ex.fem <- tapply(X=females$exs,
                 INDEX=list(Age=females$Age,
                            Year=females$Year),
                 FUN=sum)

# To avoid the noise generated by infant mortality = suggested to start at age 10
 # Note: Margin=2 in matrix stands for columns

sd10.fem <- apply(X=dx.fem[-(1:11),], MARGIN = 2, FUN = function(x)
                  sd(c(x))) 
summary(sd10.fem) 
plot(x=unique(females$Year),sd10.fem)

# from age 0
sd.fem <- apply(X=dx.fem, MARGIN = 2, FUN = function(x)
                sd(c(x))) 
plot(x=unique(females$Year),sd.fem)

###########
### IQR ###
###########
IQR.fem <- apply(X=dx.fem, MARGIN = 2, FUN = function(x)
  IQR(x))
summary(IQR.fem)
plot(x=unique(females$Year),log(IQR.fem))

IQR.fem <- as.data.frame(IQR.fem) %>% mutate(Year=seq(1876,2014,1))
  
IQR.fem %>%  ggplot(aes(x=Year,y=IQR.fem)) +
  geom_point() +
  scale_y_continuous(name="IQR") +
  theme_bw()

################
### e-Dagger ###
################

### ---------------------------------------------------------
# formula based on: www.demogr.mpg.de/papers/technicalreports/tr-2012-002.pdf
### alternative: (http://pages.stern.nyu.edu/~dbackus/BCH/demography/ZhangVaupel_ageseparating_DR_09.pdf)

e.dagger.fun <- function(lt) {
  e.dagger <- sum(lt$dxs[lt$Age]*1/2*(lt$ex[lt$Age]+lt$ex[lt$Age+1]))
  return(e.dagger)
}

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
# fem.test <- females %>% filter(Age>5) %>% group_by(Year) %>% 
#   dplyr::summarise(MDA = MDA.fun(lt))







## --------------------------------------------------------------------------- ##
## general mode
# Mode <- function(x) {
#   ux <- unique(x)
#   ux[which.max(tabulate(match(x, ux)))]
# }

# Kannisto (2001) - M(t) = x + (d(x,t)-d(x-1,t))/(2d(x,t)-d(x-1,t)-d(x+1,t))
## --------------------------------------------------------------------------- ##








###### function to fit a Gompertz model for adult mortality

# DENSITY FUNCTION Gompertz
# NOTE: x=x, shape=a, scale=b
x <- seq(5,35,1)
dwei <- function(x, a, b){
  part.1   <- (a/b)
  part.2   <- (x/b)^(a-1)
  part.3   <- exp(-(x/b)^a)
  f.x.a.b  <- part.1*part.2*part.3
  return(f.x.a.b)
}





