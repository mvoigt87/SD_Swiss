
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

# France

LT.FRA.fem <- readHMD("C:/Users/y4956294S/Documents/lt_female/fltper_1x1/FRACNP.fltper_1x1.txt", fixup = T)
LT.FRA.mal <- readHMD("C:/Users/y4956294S/Documents/lt_male/mltper_1x1/FRACNP.mltper_1x1.txt", fixup = T)


# # Spain
# 
# LT.ESP.fem <- readHMD("C:/Users/y4956294S/Documents/lt_female/fltper_1x1/ESP.fltper_1x1.txt", fixup = T)
# LT.ESP.mal <- readHMD("C:/Users/y4956294S/Documents/lt_male/mltper_1x1/ESP.mltper_1x1.txt", fixup = T)

# Italy

LT.ITA.fem <- readHMD("C:/Users/y4956294S/Documents/lt_female/fltper_1x1/ITA.fltper_1x1.txt", fixup = T)
LT.ITA.mal <- readHMD("C:/Users/y4956294S/Documents/lt_male/mltper_1x1/ITA.mltper_1x1.txt", fixup = T)

# Denmark

LT.DEN.fem <- readHMD("C:/Users/y4956294S/Documents/lt_female/fltper_1x1/DNK.fltper_1x1.txt", fixup = T)
LT.DEN.mal <- readHMD("C:/Users/y4956294S/Documents/lt_male/mltper_1x1/DNK.mltper_1x1.txt", fixup = T)


# Poland

LT.POL.fem <- readHMD("C:/Users/y4956294S/Documents/lt_female/fltper_1x1/POL.fltper_1x1.txt", fixup = T)
LT.POL.mal <- readHMD("C:/Users/y4956294S/Documents/lt_male/mltper_1x1/POL.mltper_1x1.txt", fixup = T)

# RUS

LT.RUS.fem <- readHMD("C:/Users/y4956294S/Documents/lt_female/fltper_1x1/RUS.fltper_1x1.txt", fixup = T)
LT.RUS.mal <- readHMD("C:/Users/y4956294S/Documents/lt_male/mltper_1x1/RUS.mltper_1x1.txt", fixup = T)

### ------------------------------------------------------------------------------------------------- ###
### ------------------------------------------------------------------------------------------------- ###
### ------------------------------------------------------------------------------------------------- ###

### 2. Binding the data together and show development for LE 0 and LE 65

###########
#### Age 0
###########

## Females - SUI
LT_SUI_FEM <- LT.SUI.fem %>% select(ex,Age,Year) %>% mutate(sex = "female") %>% mutate(cntry = "SUI") %>% 
  # subset only for the ex at birth (changeable)
  filter(Age == c(0,65)) %>% 
  # to assure the same end year (2012)
  filter(Year >= 1960) 
 

## Females - FRA
LT_FRA_FEM <- LT.FRA.fem %>% select(ex,Age,Year) %>% mutate(sex = "female") %>% mutate(cntry = "FRA") %>% 
  # subset only for the ex at birth (changeable)
  filter(Age == c(0,65)) %>% 
  # to assure the same end year (2012)
  filter(Year >= 1960)

# ## Females - ESP
# LT_ESP_FEM <- LT.ESP.fem %>% select(ex,Age,Year) %>% mutate(sex = "female") %>% mutate(cntry = "ESP") %>% 
#   # subset only for the ex at birth (changeable)
#   filter(Age == c(0,65)) %>% 
#   # to assure the same end year (2012)
#   filter(Year >= 1960)

## Females - ITA
LT_ITA_FEM <- LT.ITA.fem %>% select(ex,Age,Year) %>% mutate(sex = "female") %>% mutate(cntry = "ITA") %>% 
  # subset only for the ex at birth (changeable)
  filter(Age == c(0,65)) %>% 
  # to assure the same end year (2012)
  filter(Year >= 1960)

## Females - DEN
LT_DEN_FEM <- LT.DEN.fem %>% select(ex,Age,Year) %>% mutate(sex = "female") %>% mutate(cntry = "DEN") %>% 
  # subset only for the ex at birth (changeable)
  filter(Age == c(0,65)) %>% 
  # to assure the same end year (2012)
  filter(Year >= 1960)

## Females - POL
LT_POL_FEM <- LT.POL.fem %>% select(ex,Age,Year) %>% mutate(sex = "female") %>% mutate(cntry = "POL") %>% 
  # subset only for the ex at birth (changeable)
  filter(Age == c(0,65)) %>% 
  # to assure the same end year (2012)
  filter(Year >= 1960)


## Females - RUS
LT_RUS_FEM <- LT.RUS.fem %>% select(ex,Age,Year) %>% mutate(sex = "female") %>% mutate(cntry = "RUS") %>% 
  # subset only for the ex at birth (changeable)
  filter(Age == c(0,65)) %>% 
  # to assure the same end year (2012)
  filter(Year >= 1960)



## Males
########

## SUI
LT_SUI_MAL <- LT.SUI.mal %>% select(ex,Age,Year) %>% mutate(sex = "male") %>% mutate(cntry = "SUI") %>% 
  # subset only for the ex at birth (changeable)
  filter(Age == c(0,65)) %>% 
  # to assure the same end year (2012)
  filter(Year >= 1960) 


## FRA
LT_FRA_MAL <- LT.FRA.mal %>% select(ex,Age,Year) %>% mutate(sex = "male") %>% mutate(cntry = "FRA") %>% 
  # subset only for the ex at birth (changeable)
  filter(Age == c(0,65)) %>% 
  # to assure the same end year (2012)
  filter(Year >= 1960)

# ## ESP
# LT_ESP_MAL <- LT.ESP.mal %>% select(ex,Age,Year) %>% mutate(sex = "male") %>% mutate(cntry = "ESP") %>% 
#   # subset only for the ex at birth (changeable)
#   filter(Age == c(0,65)) %>% 
#   # to assure the same end year (2012)
#   filter(Year >= 1960)

## ITA
LT_ITA_MAL <- LT.ITA.mal %>% select(ex,Age,Year) %>% mutate(sex = "male") %>% mutate(cntry = "ITA") %>% 
  # subset only for the ex at birth (changeable)
  filter(Age == c(0,65)) %>% 
  # to assure the same end year (2012)
  filter(Year >= 1960)

## Females - DEN
LT_DEN_MAL <- LT.DEN.mal %>% select(ex,Age,Year) %>% mutate(sex = "male") %>% mutate(cntry = "DEN") %>% 
  # subset only for the ex at birth (changeable)
  filter(Age == c(0,65)) %>% 
  # to assure the same end year (2012)
  filter(Year >= 1960)

## POL
LT_POL_MAL <- LT.POL.mal %>% select(ex,Age,Year) %>% mutate(sex = "male") %>% mutate(cntry = "POL") %>% 
  # subset only for the ex at birth (changeable)
  filter(Age == c(0,65)) %>% 
  # to assure the same end year (2012)
  filter(Year >= 1960)


## RUS
LT_RUS_MAL <- LT.RUS.mal %>% select(ex,Age,Year) %>% mutate(sex = "male") %>% mutate(cntry = "RUS") %>% 
  # subset only for the ex at birth (changeable)
  filter(Age == c(0,65)) %>% 
  # to assure the same end year (2012)
  filter(Year >= 1960)


## Create a long dataset with all the countries

LE <- bind_rows(LT_SUI_FEM,LT_SUI_MAL,LT_DEN_FEM,LT_DEN_MAL,LT_ITA_FEM,LT_ITA_MAL,LT_FRA_FEM,LT_FRA_MAL,
                LT_POL_FEM, LT_POL_MAL, LT_RUS_FEM, LT_RUS_MAL) %>% 
  # LT_ESP_FEM,LT_ESP_MAL,
  ## Cut to a uniform age range
  #filter(Year>1950) #%>% 
  # highlight Spanish values
  mutate(highlight_flag = ifelse(cntry=="SUI",T,F))

# For plotting
LE$sex <- as.factor(LE$sex)
LE$cntry <- as.factor(LE$cntry)

summary(LE)

### ------------------------------------------------------------------------------------------------- ###
### 3. Plot (using multiplot)
### ------------------------------------------------------------------------------------------------- ###

# Life Expectancy at age zero
plotLE_zero <- LE %>% filter(Age==0) %>% ggplot() +
  # line plot
  geom_line(aes(x = Year, y = ex, color = cntry, alpha = highlight_flag))  +
  geom_point(aes(x = Year, y = ex, color = cntry, alpha = highlight_flag )) +
  scale_y_continuous(name = "Life expectancy at birth") +
  scale_x_continuous(name = " ") +
  scale_colour_manual(values = c("#26B7FF","#26FF57", "#FFD846","#FF6600","#0D3BB2","#FF6934"), name="") +
  scale_alpha_discrete(range = c(0.25, 0.85), name="", guide=F) +
  facet_wrap(~ sex) +
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

plotLE_65 <- LE %>% filter(Age==65) %>% ggplot() +
  # line plot
  geom_line(aes(x = Year, y = ex, color = cntry, alpha = highlight_flag))  +
  geom_point(aes(x = Year, y = ex, color = cntry, alpha = highlight_flag )) +
  scale_y_continuous(name = "Life expectancy at Age 65") +
  scale_x_continuous(name = " ") +
  scale_colour_manual(values = c("#26B7FF","#26FF57", "#FFD846","#FF6600","#0D3BB2","#FF6934"), name="") +
  scale_alpha_discrete(range = c(0.25, 0.85), name="", guide=F) +
  facet_wrap(~ sex) +
  theme_bw()
plotLE_65 <- plotLE_65 + theme(axis.text=element_text(size=12),
                                   axis.title=element_text(size=12,face="bold"))


### ------------------------------------------------ ###
###    For the gender gap in LE the LT by sex        ###
### ------------------------------------------------ ###


LE_M <- LE %>% filter(Age==0) %>% filter(sex=="male") %>% mutate(ex_M=ex) %>% select(- sex, - ex, - highlight_flag)
LE_F <- LE %>% filter(Age==0) %>% filter(sex=="female") %>% select(- highlight_flag)

LE_GAP <- left_join(LE_F, LE_M, by=c("Age", "Year", "cntry"))

  # calculate the gap in LE
LE_GAP <- LE_GAP %>% mutate(gap = ex - ex_M) %>% mutate(highlight_flag = ifelse(cntry=="SUI",T,F))

# For plotting
LE_GAP$cntry <- as.factor(LE_GAP$cntry)

### -------------------------------------------------- ###

# Plot the gap at the 2 time points

# Life Expectancy at age 65
plotgap <- LE_GAP %>% ggplot() +
  # line plot
  geom_line(aes(x = Year, y = gap, color=cntry, alpha = highlight_flag))  +
  geom_point(aes(x = Year, y = gap, color=cntry, alpha = highlight_flag)) +
  scale_y_continuous(name = "Female-Male Gap in LE") +
  scale_x_continuous(name = " ") +
  scale_colour_manual(values = c("#26B7FF","#26FF57", "#FFD846","#FF6600","#0D3BB2","#FF6934"), name="") +
  scale_alpha_discrete(range = c(0.25, 0.85), name="", guide=F) +
  theme_bw()

plotgap <- plotgap +  scale_shape_discrete(guide=FALSE) + theme(axis.text=element_text(size=12),
                                            axis.title=element_text(size=12,face="bold"))
# theme(legend.position = c(0.85, 0.25)) +


# Maximum Gap after 1960

maxGaP <- max(LE_GAP$gap[LE_GAP$cntry=="SUI"])
LE_GAP$Year[LE_GAP$gap==maxGaP] # 1992


###########################################################################################################################




## Calculating the gender gap in LE over time (x-female LE, y-male LE)
GAP <- function (x,y) {
  z <- x-y
  return (z)
}


# Suisse
Gap_SUI <- as.data.frame(cbind(rep(NA,29), seq(1960,2016,2),rep("SUI",29)))
names(Gap_SUI) <- c("gap","year","cntry")

for (i in min(LT_SUI_FEM$Year):max(LT_SUI_FEM$Year)) {
  Gap_SUI$gap[Gap_SUI$year==i] <- GAP(LT_SUI_FEM$ex[LT_SUI_FEM$Year==i & LT_SUI_FEM$Age==0],LT_SUI_MAL$ex[LT_SUI_MAL$Year==i & LT_SUI_MAL$Age==0])
}

# FRANCE
Gap_FRA <- as.data.frame(cbind(rep(NA,57), seq(1960,2016,1),rep("FRA",57)))
names(Gap_FRA) <- c("gap","year","cntry")

for (i in min(LT.FRA.fem$Year):max(LT.FRA.fem$Year)) {
  Gap_FRA$gap[Gap_FRA$year==i] <- GAP(LT.FRA.fem$ex[LT.FRA.fem$Year==i & LT.FRA.fem$Age==0],LT.FRA.mal$ex[LT.FRA.mal$Year==i & LT.FRA.mal$Age==0])
}

# Spain
Gap_ESP <- as.data.frame(cbind(rep(NA,57), seq(1960,2016,1),rep("ESP",57)))
names(Gap_ESP) <- c("gap","year","cntry")

for (i in min(LT.ESP.fem$Year):max(LT.ESP.fem$Year)) {
  Gap_ESP$gap[Gap_ESP$year==i] <- GAP(LT.ESP.fem$ex[LT.ESP.fem$Year==i & LT.ESP.fem$Age==0],LT.ESP.mal$ex[LT.ESP.mal$Year==i & LT.ESP.mal$Age==0])
}

# Denmark
Gap_DEN <- as.data.frame(cbind(rep(NA,57), seq(1960,2016,1), rep("DEN",57)))
names(Gap_DEN) <- c("gap","year","cntry")

for (i in min(LT.DEN.fem$Year):max(LT.DEN.fem$Year)) {
  Gap_DEN$gap[Gap_DEN$year==i] <- GAP(LT.DEN.fem$ex[LT.DEN.fem$Year==i & LT.DEN.fem$Age==0],LT.DEN.mal$ex[LT.DEN.mal$Year==i  & LT.DEN.mal$Age==0])
}

# Poland
Gap_POL <- as.data.frame(cbind(rep(NA,57), seq(1960,2016,1), rep("POL",57)))
names(Gap_POL) <- c("gap","year","cntry")

for (i in min(LT.POL.fem$Year):max(LT.POL.fem$Year)) {
  Gap_POL$gap[Gap_POL$year==i] <- GAP(LT.POL.fem$ex[LT.POL.fem$Year==i & LT.POL.fem$Age==0],LT.POL.mal$ex[LT.POL.mal$Year==i & LT.POL.mal$Age==0])
}

# Russia
Gap_RUS <- as.data.frame(cbind(rep(NA,54), seq(1960,2014,1),rep("RUS",54)))
names(Gap_RUS) <- c("gap","year","cntry")

for (i in min(LT.RUS.fem$Year):max(LT.RUS.fem$Year)) {
  Gap_RUS$gap[Gap_RUS$year==i] <- GAP(LT.RUS.fem$ex[LT.RUS.fem$Year==i & LT.RUS.fem$Age==0],LT.RUS.mal$ex[LT.RUS.mal$Year==i & LT.RUS.mal$Age==0])
}

### Combine the gaps ###
### ---------------- ###

Gap <- bind_rows(Gap_SUI, Gap_FRA, Gap_ESP, Gap_DEN, Gap_POL, Gap_RUS)

