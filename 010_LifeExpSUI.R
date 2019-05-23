
### ------------------------------------------------------------------ ###
###  Life expectancy differences - data and plots for NDS 2019         ###
### ------------------------------------------------------------------ ###

# set working directory
dir()
setwd("C:/Users/y4956294S/Documents/Workshops and Conferences/Demographic Symposium Iceland 2019/data")

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

### -------------------------------------------------- ###
### South European Countries LE at Birth and at age 50 ###
### -        -        -         -        -         -   ###
### and later for the life span disparity              ###

# Switzerland

LT.SUI.fem <- readHMD("LT_SUI_fem.txt", fixup = T)
LT.SUI.mal <- readHMD("LT_SUI_mal.txt", fixup = T)

### ------------------------------------------------------------------------------------------------- ###
### ------------------------------------------------------------------------------------------------- ###
### ------------------------------------------------------------------------------------------------- ###

### 2. Binding the data together and show development for LE 0 and LE 65

#### Age 0


## Females
LT_SUI_0_FEM <- LT.SUI.fem %>% select(ex,Age,Year) %>% mutate(sex = "female") %>% 
  # subset only for the ex at birth (changeable)
  filter(Age == 0) %>% 
  # to assure the same end year (2012)
  #filter(Year < 2013) %>% 
  # age won't be needed
  select(-Age)

## Males
LT_SUI_0_MAL <- LT.SUI.mal %>% select(ex,Age,Year) %>% mutate(sex = "male") %>% 
  # subset only for the ex at birth (changeable)
  filter(Age == 0) %>% 
  # to assure the same end year (2012)
  # filter(Year < 2013) %>% 
  # age won't be needed
  select(-Age)


## Create a long dataset with all the countries

LE_0 <- bind_rows(LT_SUI_0_FEM,LT_SUI_0_MAL)
  ## Cut to a uniform age range
  #filter(Year>1950) #%>% 
  # highlight Spanish values
  # mutate(highlight_flag = ifelse(country=="Spain",T,F))
summary(LE_0)

### ------------------------------------------------------------------------------------------------- ###

#### Age 65

## Females
LT_SUI_65_FEM <- LT.SUI.fem %>% select(ex,Age,Year) %>% mutate(sex = "female") %>% 
  # subset only for the ex at birth (changeable)
  filter(Age == 65) %>% 
  # to assure the same end year (2012)
  filter(Year > 1950) # %>% 
  # age won't be needed
  select(-Age)

## Males
LT_SUI_65_MAL <- LT.SUI.mal %>% select(ex,Age,Year) %>% mutate(sex = "male") %>% 
  # subset only for the ex at birth (changeable)
  filter(Age == 65) %>% 
  # to assure the same end year (2012)
  filter(Year > 1950) # %>% 
  # age won't be needed
  select(-Age)


## Create a long dataset with all the countries

LE_65 <- bind_rows(LT_SUI_65_FEM,LT_SUI_65_MAL)



 rm(LT_SUI_65_FEM, LT_SUI_65_MAL)

### 3. Plot (using multiplot)

# Life Expectancy at age zero
plotLE_zero <- LE_0 %>% ggplot() +
  # line plot
  geom_line(aes(x = Year, y = ex, color = sex))  +
  geom_point(aes(x = Year, y = ex, color = sex )) +
  scale_y_continuous(name = "Life expectancy at birth") +
  scale_x_continuous(name = " ") +
  scale_colour_manual(values = c("#FF6600","#0D3BB2"), name="", guide=F) +
#  scale_alpha_discrete(range = c(0.25, 0.85), name="", guide=F) +
  theme_bw()
plotLE_zero <- plotLE_zero + theme(axis.text=element_text(size=12),
                   axis.title=element_text(size=12,face="bold"))

# ----------------------------------------------------------------------------------- #
# changing colors and legends
# scale_colour_manual(values = c("#FF6934", "#26B7FF","#26FF57", "#FFD846"), name="")
# theme(legend.position = c(0.85, 0.25))
#### -----------------------------------
# colors
#        blue  = #0D3BB2
#         red  = #D52513
# yellow/gold  = #D5AD21
#      purple  = #8D11B2
#       green  = #20D51A
#### ----------------------------------
# c("#0D3BB2","#D5AD21", "#8D11B2","#20D51A", "#D52513")
# ----------------------------------------------------------------------------------- #

# Life Expectancy at age 65
plotLE_65 <- LE_65 %>% ggplot() +
  # line plot
  geom_line(aes(x = Year, y = ex, color = sex))  +
  geom_point(aes(x = Year, y = ex, color = sex)) +
  scale_y_continuous(name = "Life expectancy at age 65") +
  scale_x_continuous(name = " ") +
  scale_colour_manual(values = c("#FF6600","#0D3BB2"), name="", guide=F) +
  #scale_alpha_discrete(range = c(0.25, 0.85), name="", guide=F) +
  theme_bw()

plotLE_65 <- plotLE_65 + theme(legend.position = c(0.85, 0.25)) + 
  scale_shape_discrete(guide=FALSE) + theme(axis.text=element_text(size=12),
                                           axis.title=element_text(size=12,face="bold"))



### ------------------------------------------------ ###
###    For the gender gap in LE the LT by sex        ###

## Calculating the gender gap in LE over time (x-female LE, y-male LE)
GAP <- function (x,y) {
  z <- x-y
  return (z)
}

# Suisse
Gap_SUI <- as.data.frame(cbind(rep(NA,140), seq(1876,2016,1)))
names(Gap_SUI) <- c("gap","year")

for (i in min(LT.SUI.fem$Year):max(LT.SUI.fem$Year)) {
  Gap_SUI$gap[Gap_SUI$year==i] <- GAP(LT.SUI.fem$ex[LT.SUI.fem$Year==i],LT.SUI.mal$ex[LT.SUI.mal$Year==i])
}

### -------------------------------------------------- ###

# Plot the gap at the 2 time points

# Life Expectancy at age 65
plotgap<- Gap_SUI %>% ggplot() +
  # line plot
  geom_line(aes(x = year, y = gap, color="#FF6600"))  +
  geom_point(aes(x = year, y = gap, color="#FF6600")) +
  scale_y_continuous(name = "Female-Male Gap in LE") +
  scale_x_continuous(name = " ") +
  scale_colour_manual(values = c("#FF6600"), name="", guide=F) +
  #scale_alpha_discrete(range = c(0.25, 0.85), name="", guide=F) +
  theme_bw()

plotgap <- plotgap + theme(legend.position = c(0.85, 0.25)) + 
  scale_shape_discrete(guide=FALSE) + theme(axis.text=element_text(size=12),
                                            axis.title=element_text(size=12,face="bold"))
