### -------------------------------------------------------------------- ###
###  Cause of Death Decomposition - e0 - data and plots for NDS 2019     ###
### -------------------------------------------------------------------- ###

# set working directory
dir()
setwd("C:/Users/y4956294S/Documents/Workshops and Conferences/Demographic Symposium Iceland 2019/SD_Swiss/SD_Swiss-master")

set.seed(17952)

## use hmd/hfd package for load the data (Spain)

# LIBRARIES #
library(HMDHFDplus)
library(plyr)
library(reshape2)
library(grid)
library(gridExtra)
library(tidyverse)
library(ggplot2)
library(readxl)
library(scales)
library(RColorBrewer)
library(MortalitySmooth)
# library(MortHump)
library(latex2exp)
library(DemoDecomp)


# Functions:
source("C:/Users/y4956294S/Documents/R/EDAG_FUN.R")
source("C:/Users/y4956294S/Documents/R/LT_WHO.R")


##########################################################################################################################

# CoD data
deaths <- read.csv("C:/Users/y4956294S/Documents/Workshops and Conferences/Demographic Symposium Iceland 2019/data/SUI_COD.csv", header = T, skip = 4)

# Respective Population data
pop <- read.csv("C:/Users/y4956294S/Documents/Workshops and Conferences/Demographic Symposium Iceland 2019/data/pop.csv", header = T)

##########################################################################################################################

# generate a cause-of-death typology
causes <- unique(deaths$Cause) # all causes in the dataset
smoking <- as.character(causes[substr(causes,1,3) %in% c("C33","C34","J40", "J41", "J42", "J43", "J44", "J45", "J46", "J47", "J12", "J13", "J14", "J15", "J16", "J17", "J18")])
obesity <-  as.character(causes[substr(causes,1,3) %in% c("E10","E11","E12","E13","E14")])
alcohol <-  as.character(causes[substr(causes,1,3) %in% c("K70","K73","K74","F10") | substr(causes,1,4) %in% c("G312", "G621", "G721", "I426", "K292","K860", "Q869")])
transport <- as.character(causes[substr(causes,1,1) == "V"])
cardiovasc <- as.character(causes[substr(causes,1,3) %in% c("I10","I11","I12","I13","I15", "I20","I21","I22","I23","I24","I25", "I47", "I48", "I49", "I50", "I60","I61","I62","I63","I64","I65","I66","I67","I68","I69", "I70")])
suicide <- as.character(causes[substr(causes,1,2) %in% c("X6","X7") | substr(causes,1,3) %in% c("X80","X81","X82","X83","X84", "X45","X65")])
typ <- data.frame(name = rep(c("smoking","obesity","alcohol","transport","cardiovasc","suicide"), times = c(length(smoking), length(obesity), length(alcohol), length(transport), length(cardiovasc), length(suicide))),
                  no = rep(1:6, times = c(length(smoking), length(obesity),length(alcohol) ,length(transport), length(cardiovasc), length(suicide))), 
                  code = c(smoking, obesity, alcohol, transport, cardiovasc, suicide))


# add typology to the original dataset
deaths$typ <- typ$no[match(deaths$Cause,typ$code)]
table(deaths$typ, useNA = "always")
deaths$typ[is.na(deaths$typ)] <- length(unique(typ$no)) + 1 # deaths that do not fall into one of the types are labeled as others

# age groups
x <- as.numeric(substr(names(deaths)[9:30], 2,3))



# ------------------------------------------------------------------------------------------------------------- # 

# Objective function to compute e0 from stacked rates

e0 <- function(rates){
  
  mxc <- matrix(rates, ncol = 7, byrow = FALSE)
  
  mx <- rowSums(mxc, na.rm = TRUE)
  
  lt <- LT(Mx = mx, Age = x)
  
  return(lt$ex[1])
  
}
# ------------------------------------------------------------------------------------------------------------- #   

# Prepare color (make it colorblind friendly)

# The palette with black:
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# To use for fills, add
# scale_fill_manual(values=cbbPalette)

# To use for line and point colors, add
# scale_colour_manual(values=cbPalette)



# 3. Compute cause- and age-specific death rates for each sex from WHO data #
# ------------------------------------------------------------------------- #

# 1995 #  first year with ICD 10
# ---- #

########### FOR NOW 1995 (first year when ICD 10 was implemented)
pop_2 <- pop[pop$Country == 4300 & pop$Year == 1995,]
deaths_2 <- deaths[deaths$Year == 1995 & deaths$Cause != "A000",]

mxc.males95 <- do.call(cbind,by(data = deaths_2[deaths_2$Sex == 1,9:30], INDICES = deaths_2$typ[deaths_2$Sex == 1], FUN = colSums)) / unlist(pop_2[pop_2$Sex == 1,8:29])
mxc.females95 <- do.call(cbind,by(data = deaths_2[deaths_2$Sex == 2,9:30], INDICES = deaths_2$typ[deaths_2$Sex == 2], FUN = colSums)) / unlist(pop_2[pop_2$Sex == 2,8:29])


cod95_M <- matplot(x, log(mxc.males95), type = "l", lty= 1)
cod95_F <- matplot(x, log(mxc.females95), type = "l", lty= 1)

# stack into a vector
r.males95 <- as.vector(mxc.males95)
r.females95 <- as.vector(mxc.females95)


# Decomposition for e0 #
# -------------------- #

dec_95 <- stepwise_replacement(func = e0, pars1 = r.males95, pars2 = r.females95)

sum(dec_95) ; e0(r.females95) - e0(r.males95) # sum of decomposition is equal to difference
cxc <- matrix(dec_95, ncol = 7, byrow = FALSE) # rearange as matrix
barplot(colSums(cxc), las = 2) # contribution by cause of death
barplot(t(cxc)) # contribution by cause and age

# Make a dataframe out of it to apply ggplot

cause95 <- as.data.frame(cxc)
names(cause95) <- c("smoking", "obesity","alcohol" ,"transport", "cardiovasc" ,"suicide", "others")
# add an age vector
cause95$Age <- c(0,1,2,3,4,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85)
# add a variable for all causes combined
cause95 <- cause95 %>% mutate(allcause = smoking + obesity + alcohol + transport + cardiovasc + suicide + others)

# death rate differences over age groups
line_causes <- cause95 %>% ggplot() + geom_line(aes(x = Age, y = allcause, color="all others")) +
  geom_line(aes(x = Age, y = smoking, color="smoking")) +
  geom_line(aes(x = Age, y = obesity, color="obesity")) +
  geom_line(aes(x = Age, y = obesity, color="alcohol")) +
  geom_line(aes(x = Age, y = transport, color="transport")) +
  geom_line(aes(x = Age, y = transport, color="cardiovasc")) +
  geom_line(aes(x = Age, y = suicide, color="suicide")) +
  scale_x_continuous(name = "Age") +
  scale_y_continuous(name = "Difference in mx by cause")  +
  scale_colour_manual(values=cbbPalette, name=" ") +
  theme_bw()

line_causes + theme(axis.text=element_text(size=12),
                    axis.title=element_text(size=14,face="bold"), strip.text.y = element_text(size=12, face="bold"))


# long data format is easier to use in tidyR language

c95 <- melt(cause95,id.vars = "Age")
names(c95) <- c("Age","Cause","Contribution")

# bar plot - contribution by age and cause of death
bar_causes <- c95 %>% filter(Cause!="allcause") %>% ggplot(aes(as.factor(Age))) + 
  geom_bar(aes(weight = Contribution, fill=Cause))+
  scale_x_discrete(name = "Age") +
  scale_y_continuous(name = "Contribution by Cause")  +
  scale_fill_manual(values=cbbPalette, name=" ") +
  theme_bw()

bar_causes + theme(axis.text=element_text(size=12),
                   axis.title=element_text(size=14,face="bold"), strip.text.y = element_text(size=12, face="bold"))

# contribution by cause (combined over ages)
bar_causes2 <- c95 %>% filter(Cause!="allcause") %>% ggplot(aes(x=Cause, y=Contribution, fill=Cause)) + 
  geom_bar( stat="identity")+
  scale_x_discrete(name = "Causes") +
  scale_fill_manual(values=cbbPalette, name=" ") +                    
  #scale_fill_brewer(palette = "Set2") +
  theme_bw()

bar_causes2 + theme(axis.text=element_text(size=12),
                    axis.title=element_text(size=14,face="bold"), strip.text.y = element_text(size=12, face="bold"))


# ---- #
# 2013 #
# ---- #

# Extract rates
pop_3 <- pop[pop$Country == 4300 & pop$Year == 2013,]
deaths_3 <- deaths[deaths$Year == 2013 & deaths$Cause != "A000",]
### Next step - make the age groups in comparable (by adding up the last two age groups to the third last in 2013)
pop_3 <- pop_3 %>% mutate(Pop23 = Pop23+Pop24+Pop25) %>% mutate(Pop24=NA) %>% mutate(Pop25=NA)

# Calculate age specific and cause-specific rates
mxc.males13 <- do.call(cbind,by(data = deaths_3[deaths_3$Sex == 1,9:30], INDICES = deaths_3$typ[deaths_3$Sex == 1], FUN = colSums)) / unlist(pop_3[pop_3$Sex == 1,8:29])
mxc.females13 <- do.call(cbind,by(data = deaths_3[deaths_3$Sex == 2,9:30], INDICES = deaths_3$typ[deaths_3$Sex == 2], FUN = colSums)) / unlist(pop_3[pop_3$Sex == 2,8:29])

cod13_M <- matplot(x, log(mxc.males13), type = "l", lty= 1)
cod13_F <- matplot(x, log(mxc.females13), type = "l", lty= 1)

# stack into a vector
r.males13 <- as.vector(mxc.males13)
r.females13 <- as.vector(mxc.females13)

# age groups
x <- as.numeric(substr(names(deaths)[9:30], 2,3))

# run decomposition for e0
dec13 <- stepwise_replacement(func = e0, pars1 = r.males13, pars2 = r.females13)

sum(dec13) ; e0(r.females13) - e0(r.males13) # sum of decomposition is equal to difference
cxc <- matrix(dec13, ncol = 7, byrow = FALSE) # rearange as matrix
barplot(colSums(cxc), las = 2) # contribution by cause of death
barplot(t(cxc)) # contribution by cause and age


# 3. Plot and add color

cause13 <- as.data.frame(cxc)
names(cause13) <- c("smoking", "obesity","alcohol" ,"transport", "cardiovasc" , "suicide", "others")
# add an age vector
cause13$Age <- c(0,1,2,3,4,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85)
# add a variable for all causes combined
cause13 <- cause13 %>% mutate(allcause = smoking + obesity + alcohol + transport + cardiovasc + suicide + others)

# death rate differences over age groups
line_causes <- cause13 %>% ggplot() + geom_line(aes(x = Age, y = allcause, color="all others")) +
  geom_line(aes(x = Age, y = smoking, color="smoking")) +
  geom_line(aes(x = Age, y = obesity, color="obesity")) +
  geom_line(aes(x = Age, y = obesity, color="alcohol")) +
  geom_line(aes(x = Age, y = transport, color="transport")) +
  
  geom_line(aes(x = Age, y = suicide, color="suicide")) +
  scale_x_continuous(name = "Age") +
  scale_y_continuous(name = "Difference in mx by cause")  +
  scale_color_manual(values=cbbPalette, name=" ") +
  theme_bw()

line_causes + theme(axis.text=element_text(size=12),
                    axis.title=element_text(size=14,face="bold"), strip.text.y = element_text(size=12, face="bold"))


# long data format is easier to use in tidyR language

c13 <- melt(cause13,id.vars = "Age")
names(c13) <- c("Age","Cause","Contribution")

# bar plot - contribution by age and cause of death
bar_causes13 <- c13 %>% filter(Cause!="allcause") %>% ggplot(aes(as.factor(Age))) + 
  geom_bar(aes(weight = Contribution, fill=Cause))+
  scale_x_discrete(name = "Age") +
  scale_y_continuous(name = "Contribution by Cause")  +
  scale_fill_manual(values=cbbPalette, name=" ") +
  theme_bw()

bar_causes13 + theme(axis.text=element_text(size=12),
                     axis.title=element_text(size=14,face="bold"), strip.text.y = element_text(size=12, face="bold"))

# contribution by cause (combined over ages)
bar_causes2 <- c13 %>% filter(Cause!="allcause") %>% ggplot(aes(x=Cause, y=Contribution, fill=Cause)) + 
  geom_bar( stat="identity")+
  scale_x_discrete(name = "Causes") +
  scale_fill_manual(values=cbbPalette, name=" ") +
  theme_bw()

bar_causes2 + theme(axis.text=element_text(size=12),
                    axis.title=element_text(size=14,face="bold"), strip.text.y = element_text(size=12, face="bold"))

# -------------------------------------------------------------------------------------------------------------------- #
# Plot both years with the same scale

# Add year variable
c95 <- c95 %>% mutate(Year=1995)
c13 <- c13 %>% mutate(Year=2013)

ctotal <- union(c95, c13)


# contribution by cause (combined over ages) + facetwrap
bar_causes3 <- ctotal %>% filter(Cause!="allcause") %>% ggplot(aes(x=Cause, y=Contribution, fill=Cause)) + 
  geom_bar( stat="identity")+
  scale_x_discrete(name = "Causes") +
  scale_fill_manual(values=cbbPalette, name=" ") +
  theme_bw()

bar_causes3 + facet_grid(. ~ Year) + theme(axis.text=element_text(size=12),
                                           axis.title=element_text(size=14,face="bold"), strip.text.y = element_text(size=12, face="bold"))


# bar plot - contribution by age and cause of death
bar_causestotal <- ctotal %>% filter(Cause!="allcause") %>% ggplot(aes(as.factor(Age))) + 
  geom_bar(aes(weight = Contribution, fill=Cause))+
  scale_x_discrete(name = "Age") +
  scale_y_continuous(name = "Contribution by Cause")  +
  scale_fill_manual(values=cbbPalette, name=" ") +
  theme_bw()

bar_causestotal + facet_grid(. ~ Year) + theme(axis.text=element_text(size=12),
                                               axis.title=element_text(size=14,face="bold"), strip.text.y = element_text(size=12, face="bold"))
