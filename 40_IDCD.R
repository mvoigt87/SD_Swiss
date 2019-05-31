# Project:  Ill-defined causes of death
# Author: BWAvdL
# Created: 23-02-2018

# Load data

  if(Sys.info()[['nodename']] == "CD-5VV9QK2"){

    setwd("D:/switchdrive/SNC 2.0/Data/R")

  }else{

  setwd("C:/Users/Rose van der Linden/Documents/Projects/TOMCOD/IDCD/R")

  }

  load("SNC_00_14.RData")

  snc$status  <-  !is.na(snc$dod)

  names(snc)

  table(snc$cause_prim_icd10s)

  sum(snc$cause_prim_icd10 == "" & snc$status == 1 & snc$yod != 2014)

  # recode with hcd typology
  if(Sys.info()[['nodename']] == "CD-5VV9QK2"){

    setwd("D:/switchdrive/SNC 2.0/TOMCOD/R/IDCD")

  }else{

    setwd("C:/Users/Rose van der Linden/Documents/Projects/TOMCOD/IDCD")

  }
  typ <- read.csv("hcd.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)
  allcause <- sort(unique(c(snc$cause_prim_icd10,typ$start,typ$end)))
  head(allcause, 20)
  snc$cod <- NA

  for(i in 1:nrow(typ)){
    st <- typ$start[i]
    en <- typ$end[i]
    idx1 <- which(allcause == st)
    idx2 <- which(allcause == en)
    snc$cod[snc$cause_prim_icd10 %in% allcause[idx1:idx2]] <- typ$label[i]
  }
  snc$cod[is.na(snc$cod) & snc$status == 1] <- "Other"
  snc$cod[snc$status == 0] <- NA

table(snc$cod)

table(snc$cod, snc$sex)

table(snc$cod, snc$v0_canton)

table(snc$cod, snc$cause_prim_icd10s)


library(questionr)
cramer.v(tab = table(snc$cod, snc$cause_prim_icd10s))

snc$hcd <- snc$cod ; rm(typ)

# recode with eurostat typology
#setwd("C:/Users/Rose van der Linden/Documents/Projects/TOMCOD/IDCD")
typ <- read.csv("eurostat.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)
allcause <- sort(unique(c(snc$cause_prim_icd10,typ$start,typ$end)))
head(allcause, 20)
snc$cod <- NA

for(i in 1:nrow(typ)){
  st <- typ$start[i]
  en <- typ$end[i]
  idx1 <- which(allcause == st)
  idx2 <- which(allcause == en)
  snc$cod[snc$cause_prim_icd10 %in% allcause[idx1:idx2]] <- typ$label[i]
}
snc$cod[is.na(snc$cod) & snc$status == 1] <- "Other"
snc$cod[snc$status == 0] <- NA

table(snc$cod)

table(snc$cod, snc$sex)

table(snc$cod, snc$v0_canton)

table(snc$cod, snc$cause_prim_icd10s)


cramer.v(tab = table(snc$cod, snc$cause_prim_icd10s))

snc$eur <- snc$cod ; rm(typ)

# recode with fso typology
#setwd("C:/Users/Rose van der Linden/Documents/Projects/TOMCOD/IDCD")
typ <- read.csv("fso.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)
allcause <- sort(unique(c(snc$cause_prim_icd10,typ$start,typ$end)))
head(allcause, 20)
snc$cod <- NA

for(i in 1:nrow(typ)){
  st <- typ$start[i]
  en <- typ$end[i]
  idx1 <- which(allcause == st)
  idx2 <- which(allcause == en)
  snc$cod[snc$cause_prim_icd10 %in% allcause[idx1:idx2]] <- typ$label[i]
}
snc$cod[is.na(snc$cod) & snc$status == 1] <- "Other"
snc$cod[snc$status == 0] <- NA

table(snc$cod)

table(snc$cod, snc$sex)

table(snc$cod, snc$v0_canton)

table(snc$cod, snc$cause_prim_icd10s)


cramer.v(tab = table(snc$cod, snc$cause_prim_icd10s))

snc$sfo <- snc$cod ; rm(typ)

# recode with icd10 typology
#setwd("C:/Users/Rose van der Linden/Documents/Projects/TOMCOD/IDCD")
typ <- read.csv("icd10.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)
allcause <- sort(unique(c(snc$cause_prim_icd10,typ$start,typ$end)))
head(allcause, 20)
snc$cod <- NA

for(i in 1:nrow(typ)){
  st <- typ$start[i]
  en <- typ$end[i]
  idx1 <- which(allcause == st)
  idx2 <- which(allcause == en)
  snc$cod[snc$cause_prim_icd10 %in% allcause[idx1:idx2]] <- typ$label[i]
}
snc$cod[is.na(snc$cod) & snc$status == 1] <- "Other"
snc$cod[snc$status == 0] <- NA

table(snc$cod)

table(snc$cod, snc$sex)

table(snc$cod, snc$v0_canton)

table(snc$cod, snc$cause_prim_icd10s)


cramer.v(tab = table(snc$cod, snc$cause_prim_icd10s))

snc$icd <- snc$cod ; rm(typ)


#### compare ####

mat <- matrix(1, ncol = 4, nrow = 4)
colnames(mat) <- rownames(mat) <- names(snc)[181:184]
for(i in 1:4){
  for(j in 1:4){
    if(i == j){next}else{
      mat[i,j] <- cramer.v(tab = table(snc[,180+i], snc[,180+j]))
    }
  }
}
mat

library(corrplot)

if(Sys.info()[['nodename']] == "CD-5VV9QK2"){
  
  setwd("D:/switchdrive/SNC 2.0/TOMCOD")
  
}else{
  
  setwd("C:/Users/Rose van der Linden/Documents/Projects/TOMCOD")
  
}

png(filename = "Figures/typologies.png", width = 20, height = 20, units = "cm", res = 600)

corrplot(corr = mat, method = "color", is.corr = FALSE)

dev.off()

table(snc$hcd, snc$icd)
table(snc$icd, snc$hcd)

# conclusions
# 1. Eurostat and ICD are identical (V = 1)
# 2. SFO is the odd one out and has multually non-exclusive items, so we discart it
# 3. ICD/EUR and HCD mainly differ in the fact that circularoty and respiratory diseases are broken down in hcd
# 4. we keep hcd because circulatory diseases make up a smaller share of deaths (25% vs. 34%), so we can work with finer granularity


