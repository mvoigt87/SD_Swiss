### Applying the functions for different lifte table variability measures to to life tables of choice

## 1. Standard deviation above the mode (Cheung et al. 2009)
## 2. Coefficient of variation
## 3. Interquartile Range
## 4. E-dagger (Vaupel, Canudas-Romo)
## 5. Keyfitz' entropy

## load necessary packages

library(tidyverse)
library(HMDHFDplus)
library(MortalitySmooth)
library(MortHump)

### load smoothed life tables (decimal ages, HMD)

load("FLT_HMD.Rdata")

load("MLT_HMD.RData")

data.f <- fem.smooth[,c("Year","Age","dx","ax","lx","ex")]
data.m <- mal.smooth[,c("Year","Age","dx","ax","lx","ex")]


### load function
source("MA5_FUN.R")
source("sdfun.R")
source("Cheung_FUN.R")
source("IQR_FUN.R")
source("CV_FUN.R")
source("EDAG_FUN.R")

####################################################
#### Mode above age 5/10 based on Kannisto 2001 ####
####################################################


## apply modal age above age 5/10
# male
M.male <- by(data = data.m, INDICES = data.m$Year, FUN = MDA.fun)
# female
M.female <- by(data = data.f, INDICES = data.f$Year, FUN = MDA.fun)

## Extract information and generate dataframe for furhter analysis and plotting
# male
M.M <- as.data.frame(cbind(M.male,unique(data.m$Year))) %>% mutate(sex="male")
colnames(M.M) <- c("M10","Year","sex")
# female
M.F <- as.data.frame(cbind(M.female,unique(data.f$Year))) %>% mutate(sex="female")
colnames(M.F) <- c("M10","Year","sex")

### Plot

M10 <- bind_rows(M.M,M.F)

m10.plot <- M10 %>% ggplot(aes(x=Year, y=M10, color=sex)) +
  geom_line() +
  scale_x_continuous(name="Year") + 
  scale_y_continuous(name="Modal age above age 10") +
  scale_colour_manual(values = c("orange", "darkgrey"), name="") +
  theme_bw()


#############################################
#### Standard deviation above the mode ######
#############################################

## apply the the standard deviation above the mode to time series of HMD life tables
# male
SDM.male <- by(data = data.m, INDICES = data.m$Year, FUN = sdfun, trun = "mode", smooth = TRUE, inter = seq(0,110,0.1))
# female
SDM.female <- by(data = data.f, INDICES = data.f$Year, FUN = sdfun, trun = "mode", smooth = TRUE, inter = seq(0,110,0.1))

## Extract information and generate dataframe for furhter analysis and plotting
# male
SDM.M <- as.data.frame(cbind(SDM.male,unique(data.m$Year))) %>% mutate(sex="male")
colnames(SDM.M) <- c("SDM+","Year","sex")
# female
SDM.F <- as.data.frame(cbind(SDM.female,unique(data.f$Year))) %>% mutate(sex="female")
colnames(SDM.F) <- c("SDM+","Year","sex")

## combine the values for the two sexes

SDM <- bind_rows(SDM.M,SDM.F)

sdm.plot <- SDM %>% ggplot(aes(x=Year,y=`SDM+`,color=sex)) +
                    geom_line() +
                    scale_x_continuous(name="") +
                    scale_y_continuous(name="SD above the mode") +
                    scale_colour_manual(values = c("orange", "darkgrey"), name="") +
                    theme_bw()

#####################################################
##### and sd+ from the Cheung paper et al. 2009 #####
#####################################################


# ## apply own calculation of sd+
# # male
# SDP.male <- by(data = data.m, INDICES = data.m$Year, FUN = SD.plus.FUN)
# # female
# SDP.female <- by(data = data.f, INDICES = data.f$Year, FUN = SD.plus.FUN)
# 
# ## Extract information and generate dataframe for furhter analysis and plotting
# # male
# SDP.M <- as.data.frame(cbind(SDP.male,unique(data.m$Year))) %>% mutate(sex="male")
# colnames(SDP.M) <- c("SDM+","Year","sex")
# # female
# SDP.F <- as.data.frame(cbind(SDP.female,unique(data.f$Year))) %>% mutate(sex="female")
# colnames(SDP.F) <- c("SDM+","Year","sex")
# 
# ## combine the values for the two sexes
# 
# SDP <- bind_rows(SDP.M,SDP.F)
# 
# sdp.plot <- SDP %>% ggplot(aes(x=Year,y=`SDM+`,color=sex)) +
#   geom_line() +
#   scale_x_continuous(name="") +
#   scale_y_continuous(name="SD(M+)") +
#   scale_colour_manual(values = c("orange", "darkgrey"), name="") +
#   theme_bw()

#######################################################################################################
#######################################################################################################
#######################################################################################################

#################################################
### CV - coefficient of variation (at age 10) ###
#################################################

### single year test

data <- mal.smooth %>% select(dx,Age,Year,ex) %>% filter(Age==10)
data$Year <- as.numeric(data$Year)

CV.1990 <- CV.FUN(data[data$Year == 1990,])

# male
data.ch.m <- mal.smooth %>% select(dx,Age,Year,ex)
data.ch.m$Year <- as.numeric(data.ch.m$Year)

# female
data.ch.f <- fem.smooth %>% select(dx,Age,Year,ex)
data.ch.f$Year <- as.numeric(data.ch.f$Year)


### time series for males and females

CV.CH.M <- by(data = data.ch.m, INDICES = data.ch.m$Year, FUN = CV.FUN)
CV.CH.F <- by(data = data.ch.f, INDICES = data.ch.f$Year, FUN = CV.FUN)

### as dataframe for ggplot

CV.M <- as.data.frame(cbind(CV.CH.M,unique(data.ch.m$Year))) %>% mutate(sex="male")
colnames(CV.M) <- c("CV","Year","sex")
CV.F <- as.data.frame(cbind(CV.CH.F,unique(data.ch.f$Year))) %>% mutate(sex="female")
colnames(CV.F) <- c("CV","Year","sex")

CV.CH <- bind_rows(CV.M,CV.F)

cv.plot <- CV.CH %>% ggplot(aes(x=Year,y=CV, color=sex)) +
          geom_line() +
          scale_x_continuous(name="") +
          scale_y_continuous(name="Coefficient of Variation") +
          scale_colour_manual(values = c("orange", "darkgrey"), name="") +
          theme_bw()

#######################################################################################################
#######################################################################################################
#######################################################################################################

###################
#### apply IQR ####
###################

# male
data.ch.m <- mal.smooth %>% select(dx,Age,Year)
data.ch.m$Year <- as.numeric(data.ch.m$Year)

IQR1990 <- IQR.FUN(subset(data.ch.m,Year==1990))

# female
data.ch.f <- fem.smooth %>% select(dx,Age,Year)
data.ch.f$Year <- as.numeric(data.ch.f$Year)



### time series for males and females

IQR.CH.M <- by(data = data.ch.m, INDICES = data.ch.m$Year, FUN = IQR.FUN)
IQR.CH.F <- by(data = data.ch.f, INDICES = data.ch.f$Year, FUN = IQR.FUN)

# plot(unique(data.ch.m$Year), IQR.CH.M, type = "l", xlab = "", ylab = "IQR", las = 1)
# lines(unique(data.ch.m$Year),IQR.CH.F, col="red")


### Plot IQR (ggplot)

IQR.M <- as.data.frame(cbind(IQR.CH.M,unique(data.ch.m$Year))) %>% mutate(sex="male")
colnames(IQR.M) <- c("IQR","Year","sex")
IQR.F <- as.data.frame(cbind(IQR.CH.F,unique(data.ch.f$Year))) %>% mutate(sex="female")
colnames(IQR.F) <- c("IQR","Year","sex")

IQR <- bind_rows(IQR.F,IQR.M)

iqr.plot <- IQR %>% ggplot(aes(x=Year,y=IQR,color=sex)) +
  geom_line() +
  scale_x_continuous(name="") +
  scale_y_continuous(name="IQR in years") +
  scale_colour_manual(values = c("orange", "darkgrey"), name="") +
  theme_bw()


#######################################################################################################
#######################################################################################################
#######################################################################################################

#######################
### e+ ("e-dagger") ###
#######################

## e+(10) is much higher than it would be due to missing infant and child mortality
## might have to re-run it with life tables which start at age 0

# male
data.ch.m <- mal.smooth %>% select(dx,Age,Year,ax,ex)
data.ch.m$Year <- as.numeric(data.ch.m$Year)

# female
data.ch.f <- fem.smooth %>% select(dx,Age,Year,ax,ex)
data.ch.f$Year <- as.numeric(data.ch.f$Year)

### time series for males and females

EDAG.CH.M <- by(data = data.ch.m, INDICES = data.ch.m$Year, FUN = EDAG.FUN, smooth = TRUE)
EDAG.CH.F <- by(data = data.ch.f, INDICES = data.ch.f$Year, FUN = EDAG.FUN, smooth = TRUE)


# creating a new data frame for easier plotting and handling the summarized values 
EDAG.M <- as.data.frame(cbind(EDAG.CH.M,unique(data.ch.m$Year))) %>% mutate(sex="male")
colnames(EDAG.M) <- c("edagger","Year","sex")
EDAG.F <- as.data.frame(cbind(EDAG.CH.F,unique(data.ch.f$Year))) %>% mutate(sex="female")
colnames(EDAG.F) <- c("edagger","Year","sex")

# since the unit for edagger are decimal ages the values appeared relatively big (10 times the size)
# For now the values are presented in years (simple conversion by dividing through 10)

EDAG <- bind_rows(EDAG.M,EDAG.F) %>% mutate(edagger=edagger/10)
## E-dagger plot (reasonable)
edag.plot <- EDAG %>%  ggplot(aes(x=Year, y=edagger,color=sex))+ 
  geom_line() +
  scale_y_continuous(name = "e+(10) in years") +
  scale_colour_manual(values = c("orange", "darkgrey"), name="") +
  theme_bw()

###### Closely related = Keyfitz' entropy #######
#################################################

# e+ divided by the life expectancy

## male
H.CH.M <- bind_cols(EDAG.M, as.data.frame(data.ch.m$ex[data.ch.m$Age==10])) %>% 
  # quick calculation
  mutate(H = edagger/data.ch.m$ex[data.ch.m$Age==10])
colnames(H.CH.M)[4] <- "ex10"

## female
H.CH.F <- bind_cols(EDAG.F, as.data.frame(data.ch.f$ex[data.ch.f$Age==10])) %>% 
  # quick calculation
  mutate(H = edagger/data.ch.f$ex[data.ch.f$Age==10])
colnames(H.CH.F)[4] <- "ex10"


H.CH <- bind_rows(H.CH.F,H.CH.M)
entropy.plot <- H.CH  %>% ggplot(aes(x=Year, y=H, color=sex))+ 
  geom_line() +
  scale_y_continuous(name = "Keyfitz entropy") +
  scale_colour_manual(values = c("orange", "darkgrey"), name="") +
  theme_bw()


#######################################################################################################
#######################################################################################################
#######################################################################################################

############
### GINI ###
############



multiplot(cv.plot,iqr.plot,edag.plot,entropy.plot, cols = 2)
