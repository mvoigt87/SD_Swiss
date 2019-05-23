### ------------------------------------------------------------------ ###
###  Life span disparity differences - data and plots for NDS 2019     ###
### ------------------------------------------------------------------ ###

# set working directory
dir()
setwd("C:/Users/y4956294S/Documents/Workshops and Conferences/Demographic Symposium Iceland 2019/SD_Swiss/SD_Swiss-master")

set.seed(17952)

## use hmd/hfd package for load the data (Spain)

# LIBRARIES #
library(ggplot2)
library(gcookbook)
library(HMDHFDplus)
library(plyr)
library(reshape2)
library(grid)
library(gridExtra)
library(tidyr)
library(ggplot2)
library(readxl)
library(dplyr)
library(scales)
library(RColorBrewer)
library(MortalitySmooth)
# library(MortHump)
library(latex2exp)

### load function
source("MA5_FUN.R")
source("sdfun.R")
source("IQR_FUN.R")
source("EDAG_FUN.R")


### ------------------------------------------------------------------------------------------------- ###
### ------------------------------------------------------------------------------------------------- ###
### ------------------------------------------------------------------------------------------------- ###

### 2. Using the same life tables for calculating differences in e-dagger

# Switzerland

LT.SUI.fem <- readHMD("C:/Users/y4956294S/Documents/Workshops and Conferences/Demographic Symposium Iceland 2019/data/LT_SUI_fem.txt", fixup = T)
LT.SUI.mal <- readHMD("C:/Users/y4956294S/Documents/Workshops and Conferences/Demographic Symposium Iceland 2019/data/LT_SUI_mal.txt", fixup = T)

### 3. Calculate E-Dagger

  # females (up to age 100)
  LT.SUI.fem <- LT.SUI.fem %>% filter(Age<=100)
          # Test for 1 year
          LT.SUI.fem2000 <- LT.SUI.fem %>% filter(Year==2000)
          ED_FEM2000 <- EDAG.FUN(LT.SUI.fem2000)/100000
  # female
  ED_FEM_SUI <- by(data = LT.SUI.fem, INDICES = LT.SUI.fem$Year, FUN = EDAG.FUN)
  # This step is necessary because our formula is based on l_0 = 1 and the HMD standard l_0 = 100000
  ED_FEM_SUI <- ED_FEM_SUI/100000
  #ED_FEM_SUI <- as.matrix(ED_FEM_SUI)
  
  # males
  LT.SUI.mal <- LT.SUI.mal %>% filter(Age<=100)
  ED_MAL_SUI <- by(data = LT.SUI.mal, INDICES = LT.SUI.mal$Year, FUN = EDAG.FUN)
  # This step is necessary because our formula is based on l_0 = 1 and the HMD standard l_0 = 100000
  ED_MAL_SUI <- ED_MAL_SUI/100000
  #ED_MAL_SUI <- as.matrix(ED_MAL_SUI)
  
  # data frame
  
  edag_both_SUI <- as.data.frame(cbind(unique(LT.SUI.fem$Year),ED_FEM_SUI, ED_MAL_SUI))
  
  colnames(edag_both_SUI)[1] <- "Year"
  
  # plot edagger Switzerland
  edag_both_plot <- edag_both_SUI %>% ggplot(aes(x=Year)) +
    geom_line(aes(y = ED_MAL_SUI, color="Male")) + 
    geom_line(aes(y = ED_FEM_SUI, color="Female")) +
    scale_y_continuous(name = TeX('$e^\\dagger$')) +
    scale_color_manual(name=" ", values = c("#FF6600","#0D3BB2"), guide=F) +
    theme_bw()
  
  edag_both_plot <- edag_both_plot + theme(axis.text=element_text(size=12),
                                           axis.title=element_text(size=12,face="bold"))
  



## Calculating the sex gap in disparity over time
  
GAP <- function (x,y) {
  z <- x-y
  return (z)
}

# Suisse
Gap_SUI_LSD <- as.data.frame(cbind(rep(NA,140), seq(1876,2016,1)))
names(Gap_SUI_LSD) <- c("gap","year")

for (i in min(LT.SUI.fem$Year):max(LT.SUI.fem$Year)) {
  Gap_SUI_LSD$gap[Gap_SUI_LSD$year==i] <- GAP(edag_both_SUI$ED_MAL_SUI[edag_both_SUI$Year==i],edag_both_SUI$ED_FEM_SUI[edag_both_SUI$Year==i])
}

# Plot GAP
plotgap_LSD <- Gap_SUI_LSD %>% ggplot() +
  # line plot
  geom_line(aes(x = year, y = gap, color="#FF6600"))  +
  geom_point(aes(x = year, y = gap, color="#FF6600")) +
  scale_y_continuous(name = "Male-Female Gap in LSD") +
  scale_x_continuous(name = " ") +
  scale_colour_manual(values = c("#FF6600"), name="", guide=F) +
  #scale_alpha_discrete(range = c(0.25, 0.85), name="", guide=F) +
  theme_bw()

plotgap_LSD <- plotgap_LSD + theme(legend.position = c(0.85, 0.25)) + 
  scale_shape_discrete(guide=FALSE) + theme(axis.text=element_text(size=12),
                                            axis.title=element_text(size=12,face="bold"))
