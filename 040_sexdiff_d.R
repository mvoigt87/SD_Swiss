#########################################################
## Descriptive Stats for explaining life table variation 
#########################################################

# 1. Mortality surfaces for men and women

# 2. Lifespan variation measures

library(tidyverse)
library(HMDHFDplus)
library(MortalitySmooth)
library(MortHump)


set.seed(17952)

### load smoothed (interpolated) life tables (decimal ages, source: HMD)

load("FLT_HMD.Rdata")

load("MLT_HMD.RData")

data.f <- fem.smooth[,c("Year","Age","dx","ax","lx","ex")]
data.m <- mal.smooth[,c("Year","Age","dx","ax","lx","ex")]


### load function for life span disparity
# source("MA5_FUN.R")
# source("sdfun.R")
# source("Cheung_FUN.R")
source("IQR_FUN.R")
source("CV_FUN.R")
source("EDAG_FUN.R")

###############################################################################
#### Sex Gap in Life Expectancy and Mortality Surfaces using the mx values ####
###############################################################################


## Calculating the gender gap in LE over time (x-female LE, y-male LE)
    GAP <- function (x,y) {
      z <- x-y
      return (z)
    }

## Switzerland - Age 0
## -------------------

  # SUI_0 <- as.data.frame(cbind(rep(NA,139), seq(1876,2014,1)))
  # names(SUI_0) <- c("gap","year")
  # 
  # for (i in min(data.f$Year):max(data.f$Year)) {
  #   SUI_0$gap[SUI_0$year==i] <- GAP(data.f$ex[data.f$Year==i & data.f$Age==0],data.m$ex[data.m$Year==i & data.m$Age==0])
  # }

## Switzerland - Age 65
## ---------------------

SUI_65 <- as.data.frame(cbind(rep(NA,139), seq(1876,2014,1)))
names(SUI_65) <- c("gap","year")

for (i in min(data.f$Year):max(data.f$Year)) {
  SUI_65$gap[SUI_65$year==i] <- GAP(data.f$ex[data.f$Year==i & data.f$Age==65],data.m$ex[data.m$Year==i & data.m$Age==65])
}

## Plotting the gap in LE (SUI) 
## ----------------------------

GAP_SUI <- ggplot(SUI_65, aes(year)) +  geom_line(aes(y = gap, colour = " ")) 
GAP_SUI  + ggtitle("Sex Gap in LE at age 65") + 
  scale_y_continuous(name="Gender Gap in LE in years") + 
  scale_x_continuous(name="Year", breaks=c(1875, 1900, 1925, 1950, 1975, 2000)) +
  scale_color_discrete(guide=FALSE) +
#  labs(colour=" " , guide=FALSE) + 
  theme_bw()


### Mortality Surfaces
### ------------------

## generate the basic matrices (rows = age groups, columns = years)

SUI_FEM_qx <- matrix(NA, nrow=24, ncol=56)
SUI_MAL_qx<- matrix(NA, nrow=24, ncol=56)

### fill the matrices
## females
for (i in 1:56){
  EGER.F.ma.5[,i] <- EGER.fem.5$qx[EGER.fem.5$Year==1955+i]
}

### males
for (i in 1:56){
  EGER.M.ma.5[,i] <- EGER.mal.5$qx[EGER.mal.5$Year==1955+i]
}

### gap = male - female mortality rates (in order to have positive numbers)
EGER.GAP.ma.5 <- EGER.M.ma.5 - EGER.F.ma.5

## Create a dataEGERme for ggplot - only for the >5 population
EGER.GAP.DF.5 <- data.frame(
  mx.gap = as.numeric(EGER.GAP.ma.5),
  year = as.numeric(rep(1956:2011, each=24)),
  age = as.numeric(rep(seq(0,23,1),56))
)



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
