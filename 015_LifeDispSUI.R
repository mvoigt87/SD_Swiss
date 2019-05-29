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

# France

LT.FRA.fem <- readHMD("C:/Users/y4956294S/Documents/lt_female/fltper_1x1/FRACNP.fltper_1x1.txt", fixup = T)
LT.FRA.mal <- readHMD("C:/Users/y4956294S/Documents/lt_male/mltper_1x1/FRACNP.mltper_1x1.txt", fixup = T)


# Spain

LT.ESP.fem <- readHMD("C:/Users/y4956294S/Documents/lt_female/fltper_1x1/ESP.fltper_1x1.txt", fixup = T)
LT.ESP.mal <- readHMD("C:/Users/y4956294S/Documents/lt_male/mltper_1x1/ESP.mltper_1x1.txt", fixup = T)


# Denmark

LT.DEN.fem <- readHMD("C:/Users/y4956294S/Documents/lt_female/fltper_1x1/DNK.fltper_1x1.txt", fixup = T)
LT.DEN.mal <- readHMD("C:/Users/y4956294S/Documents/lt_male/mltper_1x1/DNK.mltper_1x1.txt", fixup = T)


# Poland

LT.POL.fem <- readHMD("C:/Users/y4956294S/Documents/lt_female/fltper_1x1/POL.fltper_1x1.txt", fixup = T)
LT.POL.mal <- readHMD("C:/Users/y4956294S/Documents/lt_male/mltper_1x1/POL.mltper_1x1.txt", fixup = T)

# RUS

LT.RUS.fem <- readHMD("C:/Users/y4956294S/Documents/lt_female/fltper_1x1/RUS.fltper_1x1.txt", fixup = T)
LT.RUS.mal <- readHMD("C:/Users/y4956294S/Documents/lt_male/mltper_1x1/RUS.mltper_1x1.txt", fixup = T)

### 3. Calculate E-Dagger

### Switzerland
### -----------

  # females (up to age 100)
  LT.SUI.fem <- LT.SUI.fem %>% filter(Age<=100) %>% filter(Year >= 1960) 
          # Test for 1 year
          LT.SUI.fem2000 <- LT.SUI.fem %>% filter(Year==2000)
          ED_FEM2000 <- EDAG.FUN(LT.SUI.fem2000)/100000
          
  # female
  ED_FEM_SUI <- by(data = LT.SUI.fem, INDICES = LT.SUI.fem$Year, FUN = EDAG.FUN)
  # This step is necessary because our formula is based on l_0 = 1 and the HMD standard l_0 = 100000
  ED_FEM_SUI <- ED_FEM_SUI/100000
  #ED_FEM_SUI <- as.matrix(ED_FEM_SUI)
  
  # males
  LT.SUI.mal <- LT.SUI.mal %>% filter(Age<=100) %>% filter(Year >= 1960) 
  ED_MAL_SUI <- by(data = LT.SUI.mal, INDICES = LT.SUI.mal$Year, FUN = EDAG.FUN)
  # This step is necessary because our formula is based on l_0 = 1 and the HMD standard l_0 = 100000
  ED_MAL_SUI <- ED_MAL_SUI/100000
  #ED_MAL_SUI <- as.matrix(ED_MAL_SUI)
  
  # bind data frame and give a country stamp
  
  edag_both_SUI <- as.data.frame(cbind(unique(LT.SUI.fem$Year),ED_FEM_SUI, ED_MAL_SUI))
  
  colnames(edag_both_SUI) <- c("Year","ed_fem","ed_mal")
  
  edag_both_SUI <- edag_both_SUI %>% mutate(cntry = "SUI")
  
  
  ### France
  ### -------
  
  # females (up to age 100)
  LT.FRA.fem <- LT.FRA.fem %>% filter(Age<=100) %>% filter(Year >= 1960) 
  # Test for 1 year
  LT.FRA.fem2000 <- LT.FRA.fem %>% filter(Year==2000)
  ED_FEM2000 <- EDAG.FUN(LT.FRA.fem2000)/100000
  
  # female
  ED_FEM_FRA <- by(data = LT.FRA.fem, INDICES = LT.FRA.fem$Year, FUN = EDAG.FUN)
  # This step is necessary because our formula is based on l_0 = 1 and the HMD standard l_0 = 100000
  ED_FEM_FRA <- ED_FEM_FRA/100000
  #ED_FEM_FRA <- as.matrix(ED_FEM_FRA)
  
  # males
  LT.FRA.mal <- LT.FRA.mal %>% filter(Age<=100) %>% filter(Year >= 1960) 
  ED_MAL_FRA <- by(data = LT.FRA.mal, INDICES = LT.FRA.mal$Year, FUN = EDAG.FUN)
  # This step is necessary because our formula is based on l_0 = 1 and the HMD standard l_0 = 100000
  ED_MAL_FRA <- ED_MAL_FRA/100000
  #ED_MAL_FRA <- as.matrix(ED_MAL_FRA)
  
  # bind data frame and give a country stamp
  
  edag_both_FRA <- as.data.frame(cbind(unique(LT.FRA.fem$Year),ED_FEM_FRA, ED_MAL_FRA))
  
  colnames(edag_both_FRA) <- c("Year","ed_fem","ed_mal")
  
  edag_both_FRA <- edag_both_FRA %>% mutate(cntry = "FRA")
  
  
  ### Spain
  ### -----
  
  # females (up to age 100)
  LT.ESP.fem <- LT.ESP.fem %>% filter(Age<=100) %>% filter(Year >= 1960) 
  # Test for 1 year
  LT.ESP.fem2000 <- LT.ESP.fem %>% filter(Year==2000)
  ED_FEM2000 <- EDAG.FUN(LT.ESP.fem2000)/100000
  
  # female
  ED_FEM_ESP <- by(data = LT.ESP.fem, INDICES = LT.ESP.fem$Year, FUN = EDAG.FUN)
  # This step is necessary because our formula is based on l_0 = 1 and the HMD standard l_0 = 100000
  ED_FEM_ESP <- ED_FEM_ESP/100000
  #ED_FEM_ESP <- as.matrix(ED_FEM_ESP)
  
  # males
  LT.ESP.mal <- LT.ESP.mal %>% filter(Age<=100) %>% filter(Year >= 1960) 
  ED_MAL_ESP <- by(data = LT.ESP.mal, INDICES = LT.ESP.mal$Year, FUN = EDAG.FUN)
  # This step is necessary because our formula is based on l_0 = 1 and the HMD standard l_0 = 100000
  ED_MAL_ESP <- ED_MAL_ESP/100000
  #ED_MAL_ESP <- as.matrix(ED_MAL_ESP)
  
  # bind data ESPme and give a country stamp
  
  edag_both_ESP <- as.data.frame(cbind(unique(LT.ESP.fem$Year),ED_FEM_ESP, ED_MAL_ESP))
  
  colnames(edag_both_ESP) <- c("Year","ed_fem","ed_mal")
  
  edag_both_ESP <- edag_both_ESP %>% mutate(cntry = "ESP") 
  
  
  ### Denmark
  ### -------
  
  # females (up to age 100)
  LT.DEN.fem <- LT.DEN.fem %>% filter(Age<=100) %>% filter(Year >= 1960) 
  # Test for 1 year
  LT.DEN.fem2000 <- LT.DEN.fem %>% filter(Year==2000)
  ED_FEM2000 <- EDAG.FUN(LT.DEN.fem2000)/100000
  
  # female
  ED_FEM_DEN <- by(data = LT.DEN.fem, INDICES = LT.DEN.fem$Year, FUN = EDAG.FUN)
  # This step is necessary because our formula is based on l_0 = 1 and the HMD standard l_0 = 100000
  ED_FEM_DEN <- ED_FEM_DEN/100000
  #ED_FEM_DEN <- as.matrix(ED_FEM_DEN)
  
  # males
  LT.DEN.mal <- LT.DEN.mal %>% filter(Age<=100) %>% filter(Year >= 1960) 
  ED_MAL_DEN <- by(data = LT.DEN.mal, INDICES = LT.DEN.mal$Year, FUN = EDAG.FUN)
  # This step is necessary because our formula is based on l_0 = 1 and the HMD standard l_0 = 100000
  ED_MAL_DEN <- ED_MAL_DEN/100000
  #ED_MAL_DEN <- as.matrix(ED_MAL_DEN)
  
  # bind data DENme and give a country stamp
  
  edag_both_DEN <- as.data.frame(cbind(unique(LT.DEN.fem$Year),ED_FEM_DEN, ED_MAL_DEN))
  
  colnames(edag_both_DEN) <- c("Year","ed_fem","ed_mal")
  
  edag_both_DEN <- edag_both_DEN %>% mutate(cntry = "DEN")  
  
  
  
  ### Poland
  ### ------
  
  # females (up to age 100)
  LT.POL.fem <- LT.POL.fem %>% filter(Age<=100) %>% filter(Year >= 1960) 
  # Test for 1 year
  LT.POL.fem2000 <- LT.POL.fem %>% filter(Year==2000)
  ED_FEM2000 <- EDAG.FUN(LT.POL.fem2000)/100000
  
  # female
  ED_FEM_POL <- by(data = LT.POL.fem, INDICES = LT.POL.fem$Year, FUN = EDAG.FUN)
  # This step is necessary because our formula is based on l_0 = 1 and the HMD standard l_0 = 100000
  ED_FEM_POL <- ED_FEM_POL/100000
  #ED_FEM_POL <- as.matrix(ED_FEM_POL)
  
  # males
  LT.POL.mal <- LT.POL.mal %>% filter(Age<=100) %>% filter(Year >= 1960) 
  ED_MAL_POL <- by(data = LT.POL.mal, INDICES = LT.POL.mal$Year, FUN = EDAG.FUN)
  # This step is necessary because our formula is based on l_0 = 1 and the HMD standard l_0 = 100000
  ED_MAL_POL <- ED_MAL_POL/100000
  #ED_MAL_POL <- as.matrix(ED_MAL_POL)
  
  # bind data POLme and give a country stamp
  
  edag_both_POL <- as.data.frame(cbind(unique(LT.POL.fem$Year),ED_FEM_POL, ED_MAL_POL))
  
  colnames(edag_both_POL) <- c("Year","ed_fem","ed_mal")
  
  edag_both_POL <- edag_both_POL %>% mutate(cntry = "POL") 
  
  
  ### Russia
  ### ------
  
  # females (up to age 100)
  LT.RUS.fem <- LT.RUS.fem %>% filter(Age<=100) %>% filter(Year >= 1960) 
  # Test for 1 year
  LT.RUS.fem2000 <- LT.RUS.fem %>% filter(Year==2000)
  ED_FEM2000 <- EDAG.FUN(LT.RUS.fem2000)/100000
  
  # female
  ED_FEM_RUS <- by(data = LT.RUS.fem, INDICES = LT.RUS.fem$Year, FUN = EDAG.FUN)
  # This step is necessary because our formula is based on l_0 = 1 and the HMD standard l_0 = 100000
  ED_FEM_RUS <- ED_FEM_RUS/100000
  #ED_FEM_RUS <- as.matrix(ED_FEM_RUS)
  
  # males
  LT.RUS.mal <- LT.RUS.mal %>% filter(Age<=100) %>% filter(Year >= 1960) 
  ED_MAL_RUS <- by(data = LT.RUS.mal, INDICES = LT.RUS.mal$Year, FUN = EDAG.FUN)
  # This step is necessary because our formula is based on l_0 = 1 and the HMD standard l_0 = 100000
  ED_MAL_RUS <- ED_MAL_RUS/100000
  #ED_MAL_RUS <- as.matrix(ED_MAL_RUS)
  
  # bind data RUSme and give a country stamp
  
  edag_both_RUS <- as.data.frame(cbind(unique(LT.RUS.fem$Year),ED_FEM_RUS, ED_MAL_RUS))
  
  colnames(edag_both_RUS) <- c("Year","ed_fem","ed_mal")
  
  edag_both_RUS <- edag_both_RUS %>% mutate(cntry = "RUS") 
  
  
  
  
  
  
  ## Create a long dataset with all the countries
  
  LSD <- bind_rows(edag_both_SUI,edag_both_FRA, edag_both_ESP, edag_both_DEN, edag_both_POL, edag_both_RUS) %>% 
    ## Cut to a uniform age range
    #filter(Year>1950) #%>% 
    # highlight Spanish values
    mutate(highlight_flag = ifelse(cntry=="SUI",T,F))
  
  # For plotting
  #LSD$sex <- as.factor(LSD$sex)
  LSD$cntry <- as.factor(LSD$cntry)
  
  summary(LSD)

# ------------ #
# plot edagger #
# ------------ #
  
edag_plot <- LSD %>% ggplot() +
    geom_line(aes(x = Year, y = ed_fem, color = cntry, alpha = highlight_flag), linetype=1) + 
    geom_line(aes(x = Year, y = ed_mal, color = cntry, alpha = highlight_flag), linetype=2) +
    scale_y_continuous(name = TeX('$e^\\dagger$')) +
    scale_colour_manual(values = c("#26B7FF","#26FF57", "#FFD846","#FF6600","#0D3BB2","#FF6934"), name="") +
    scale_alpha_discrete(range = c(0.35, 0.85), name="", guide=F) +
    theme_bw()
  
edag_plot <- edag_plot + theme(axis.text=element_text(size=12), axis.title=element_text(size=12,face="bold"))


## ----------------------------------------------- ##
## Calculating the sex gap in disparity over time  ##
## ----------------------------------------------- ##

# calculate the gap in LE
LSD <- LSD %>% mutate(gap = ed_mal - ed_fem)

# Life Expectancy at age 65
plotgap_LSD <- LSD %>% ggplot() +
  # line plot
  geom_line(aes(x = Year, y = gap, color=cntry, alpha = highlight_flag))  +
  geom_point(aes(x = Year, y = gap, color=cntry, alpha = highlight_flag)) +
  scale_y_continuous(name = "Male-Female Gap in Lifespan Disparity") +
  scale_x_continuous(name = " ") +
  scale_colour_manual(values = c("#26B7FF","#26FF57", "#FFD846","#FF6600","#0D3BB2","#FF6934"), name="") +
  scale_alpha_discrete(range = c(0.25, 0.85), name="", guide=F) +
  theme_bw()

plotgap_LSD <- plotgap_LSD +  scale_shape_discrete(guide=FALSE) + theme(axis.text=element_text(size=12),
                                                                axis.title=element_text(size=12,face="bold"))











###########################################################################################################################





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
