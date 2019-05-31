# Decomposition of changes in Life Expectancy (LE) and Lifespan Disparities (LD)
# 2000 - 2014
# Adrien Remund, May 2018

library(RColorBrewer)
cols <- c(brewer.pal(n = 9, name = "Set1"),brewer.pal(n = 8, name = "Set3"))

#### load data ####

if(Sys.info()[['nodename']] == "CD-5VV9QK2"){
  
  setwd("D:/switchdrive/SNC 2.0/Data/R")
  
}else{ if(Sys.info()[['nodename']] == "REMUND"){
  
  setwd("C:/Users/Adrien/switchdrive/SNC 2.0/Data/R")
  
}else{
  
  setwd("C:/Users/Rose van der Linden/switchdrive/SNC 2.0/Data/R")
  
}}

load("SNC_00_14.RData")

#### clean snc data ####

snc$status  <-  !is.na(snc$dod) 

if(Sys.info()[['nodename']] == "CD-5VV9QK2"){
  
  setwd("D:/switchdrive/SNC 2.0/TOMCOD")
  
}else{ if(Sys.info()[['nodename']] == "REMUND"){
  
  setwd("C:/Users/Adrien/switchdrive/SNC 2.0/TOMCOD")
  
}else{
  
  setwd("C:/Users/Rose van der Linden/switchdrive/SNC 2.0/TOMCOD")
  
}}

# recode with HCD typology
typ <- read.csv("R/IDCD/hcd.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)
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
# missing causes of death are treated as ill-defined
snc$cod[is.na(snc$cod) & snc$status == 1] <- "Ill-defined"
snc$cod[snc$status == 0] <- NA


# since cod is missing in 2014, censor on 31.12.2013
# and recode remaining 88 Missings into Other causes
table(snc$cod, snc$yod)
snc$status[snc$yod == 2014] <- 0
snc$dstop_imputed[snc$dstop_imputed > as.Date("2013-12-31")] <- as.Date("2014-01-01")
snc$cod[snc$yod == 2014] <- NA

# reduce to necessary variables
data <- snc[,c("snc_ge_id","status","imputed","link","sex","dob","dod","dstop_imputed","stopcode_imputed","cod","v0_canton")]
data$entry <- as.Date("2000-12-05")
data <- data[,c("snc_ge_id", "sex", "dob", "status", "entry", "dstop_imputed", "cod","v0_canton")]
names(data) <- c("id","sex","birth","fail","entry","exit","cod","canton")
data$fail <- as.numeric(data$fail)
cantons  <- read.csv("R/cantons.csv", sep = ";", header = TRUE)
data$canton <- cantons$code[match(data$canton, cantons$canton)]

rm(list = ls()[-which(ls() %in% c("data","cols"))])

#### compute rates ####

years <- 2000:2014

x <- seq(0,110,5)

dxc_m <- table(cut(as.numeric((data$exit - data$birth))[data$sex == "Male"] / 365.25, breaks = x, right = FALSE), data$cod[data$sex == "Male"], cut(lubridate::year(data$exit[data$sex == "Male"]), breaks = years, right = FALSE))

dxc_f <- table(cut(as.numeric((data$exit - data$birth))[data$sex == "Female"] / 365.25, breaks = x, right = FALSE), data$cod[data$sex == "Female"], cut(lubridate::year(data$exit[data$sex == "Female"]), breaks = years, right = FALSE))

#### redistribute ill-defined causes ####

load(file = "R/IDCD/redistribute.RData")

coefs
colSums(coefs)

apply(dxc_m, 2, sum)

for(i in 1:length(years)){
  dxc_m[,,i] <- dxc_m[,,i] + dxc_m[,rownames(coefs) == "Ill-defined",i] %*% t(coefs[,1])
  dxc_m[,rownames(coefs) == "Ill-defined",i] <- 0
  dxc_f[,,i] <- dxc_f[,,i] + dxc_f[,rownames(coefs) == "Ill-defined",i] %*% t(coefs[,2])
  dxc_f[,rownames(coefs) == "Ill-defined",i] <- 0
}

apply(dxc_m, 2, sum)

dxc_m <- dxc_m[,-which(rownames(coefs) == "Ill-defined"),]
dxc_f <- dxc_f[,-which(rownames(coefs) == "Ill-defined"),]
cols <- cols[-which(rownames(coefs) == "Ill-defined")]

png("Figures/pxc_m.png", width = 30, height = 50, units = "cm", res = 600)
par(mfrow = MortHump:::disp(length(years)))
for(i in 1:(length(years) - 1)){
  barplot(t(dxc_m[,,i] / rowSums(dxc_m[,,i])), col = cols, main = years[i], las = 2, cex.names = 0.8, space = 0, border = NA)
}
plot.new()
legend("center", legend = rownames(coefs), fill = cols)
par(mfrow = c(1,1))
dev.off()

png("Figures/pxc_f.png", width = 30, height = 50, units = "cm", res = 600)
par(mfrow = MortHump:::disp(length(years)))
for(i in 1:(length(years) - 1)){
  barplot(t(dxc_f[,,i] / rowSums(dxc_f[,,i])), col = cols, main = years[i], las = 2, cex.names = 0.8, space = 0, border = NA)
}
plot.new()
legend("center", legend = rownames(coefs), fill = cols)
par(mfrow = c(1,1))
dev.off()

#### proportions ####

pxc_m <- dxc_m # apply(dxc_m, 3, function(x){x / rowSums(x)})
for(i in 1:length(years)){pxc_m[,,i] <- pxc_m[,,i] / rowSums(pxc_m[,,i])}

pxc_f <- dxc_f # apply(dxc_f, 3, function(x){x / rowSums(x)})
for(i in 1:length(years)){pxc_f[,,i] <- pxc_f[,,i] / rowSums(pxc_f[,,i])}

xs <- seq(0,110,0.1)

pxc_f_smooth <- pxc_m_smooth <- array(NA, dim = c(length(xs),dim(pxc_f)[2:3]))


# library(MortalitySmooth)
# fit <- Mort2Dsmooth(x = x[-1], y = years[-1], Z = pxc_f[,2,], lambdas = c(1e-5,1e-5), method = 3)
# fit$lambdas
# matplot(x = x[-length(x)], pxc_f[,2,], type = "p", pch = 16, lty = 1, col = heat.colors(length(years)), ylim = c(0,0.2))
# matlines(x = x[-length(x)], exp(predict(fit)), lty = 1, lwd = 2, col = heat.colors(length(years)))


pdf(file = "Figures/smooth_pxc_f.pdf", width = 10, height = 10)
for(i in c(1:dim(pxc_f)[2])){
 plot(0, 0, type = "n", xlim = range(xs), ylim = range(pxc_f[,i,1], na.rm = TRUE), las = 1, xlab = "", ylab = "%", main = colnames(pxc_f[,,1])[i])
 l <- list()
  for(j in 1:dim(pxc_f)[3]){
    y1 <- na.omit(pxc_f[,i,j])
    if(is.null(attr(y1,"na.action"))){x1 <- seq(0,105,5) ; w1 <- dxc_f[,i,j]}else{x1 <- seq(0,105,5)[-attr(y1,"na.action")] ; w1 <- dxc_f[,i,j][-attr(y1,"na.action")]}
    fit <- smooth.spline(x = x1, y = y1, w = w1, control.spar = list(low = ifelse(i == 7,-0.5,0.5), high = ifelse(i == 7,0.2,1)))
    l[[i]] <- fit$spar
    ys <- predict(fit, x = xs)$y
    ys[ys < 0] <- 0
    points(x1, y1, col = heat.colors(length(years))[j])
    lines(xs, ys, lwd = 2, col = heat.colors(length(years))[j])
    pxc_f_smooth[,i,j] <- ys
  }  
 text(0,0, labels = round(mean(unlist(l)),2))
 }
dev.off()

# rescale to 100%
apply(pxc_f_smooth, 3, rowSums, na.rm = TRUE)
for(i in 1:length(years)){
  
  pxc_f_smooth[,,i] <- pxc_f_smooth[,,i] / rowSums(pxc_f_smooth[,,i], na.rm = TRUE)
  
}

png("Figures/pxc_f_smooth.png", width = 20, height = 70, units = "cm", res = 600)
par(mfrow = c(7,2))
for(i in 1:length(years)){
  matplot(xs, pxc_f_smooth[,,i], type = "l", lty = 1, lwd = 2, las = 1, ylab = "%", xlab = "age", col = cols, main = years[i])
  matlines(x[-length(x)], pxc_f[,,i], type = "l", lty = 1, col = cols)
}
par(mfrow = c(1,1))
dev.off()

dim(pxc_f_smooth)
pxc_f_smooth <- do.call(rbind,lapply(seq(dim(pxc_f_smooth)[3]), function(x){pxc_f_smooth[ , , x]}))
dim(pxc_f_smooth)
colnames(pxc_f_smooth) <- colnames(pxc_f[,,1])


pdf(file = "Figures/smooth_pxc_m.pdf", width = 10, height = 10)
for(i in c(1:dim(pxc_f)[2])){
  plot(0, 0, type = "n", xlim = range(xs), ylim = range(pxc_m[,i,1], na.rm = TRUE), las = 1, xlab = "", ylab = "%", main = colnames(pxc_m[,,1])[i])
  l <- list()
  for(j in 1:dim(pxc_m)[3]){
    y1 <- na.omit(pxc_m[,i,j])
    if(is.null(attr(y1,"na.action"))){x1 <- seq(0,105,5) ; w1 <- dxc_f[,i,j]}else{x1 <- seq(0,105,5)[-attr(y1,"na.action")] ; w1 <- dxc_m[,i,j][-attr(y1,"na.action")]}
    fit <- smooth.spline(x = x1, y = y1, w = w1, control.spar = list(low = ifelse(i == 7,-0.5,0.5), high = ifelse(i == 7,0.2,1)))
    l[[i]] <- fit$spar
    ys <- predict(fit, x = xs)$y
    ys[ys < 0] <- 0
    points(x1, y1, col = heat.colors(length(years))[j])
    lines(xs, ys, lwd = 2, col = heat.colors(length(years))[j])
    pxc_m_smooth[,i,j] <- ys
  } 
  text(0,0, labels = round(mean(unlist(l)),2))
}
dev.off()

# rescale to 100%
colMeans(apply(pxc_m_smooth, 3, rowSums, na.rm = TRUE))
for(i in 1:length(years)){
  
  pxc_m_smooth[,,i] <- pxc_m_smooth[,,i] / rowSums(pxc_m_smooth[,,i], na.rm = TRUE)
  
}

png("Figures/pxc_m_smooth.png", width = 20, height = 70, units = "cm", res = 600)
par(mfrow = c(7,2))
for(i in 1:length(years)){
  matplot(xs, pxc_m_smooth[,,i], type = "l", lty = 1, lwd = 2, las = 1, ylab = "%", xlab = "age", col = cols, main = years[i])
  matlines(x[-length(x)], pxc_m[,,i], type = "l", lty = 1, col = cols)
}
par(mfrow = c(1,1))
dev.off()

dim(pxc_m_smooth)
pxc_m_smooth <- do.call(rbind,lapply(seq(dim(pxc_m_smooth)[3]), function(x){pxc_m_smooth[ , , x]}))
dim(pxc_m_smooth)
colnames(pxc_m_smooth) <- colnames(pxc_m[,,1])

rm(list = ls()[-which(ls() %in% c("cols","pxc_m_smooth","pxc_f_smooth","years","x","xs"))])

#### apply proportions to interpolated data ####

load("R/interpolated/FLT_HMD.Rdata")

load("R/interpolated/MLT_HMD.RData")

fem.smooth <- fem.smooth[fem.smooth$Year %in% years[-length(years)],]
mal.smooth <- mal.smooth[mal.smooth$Year %in% years[-length(years)],]

data.f <- cbind(fem.smooth, matrix(NA, ncol = dim(pxc_f_smooth)[2], nrow = nrow(fem.smooth)))
data.m <- cbind(mal.smooth, matrix(NA, ncol = dim(pxc_m_smooth)[2], nrow = nrow(mal.smooth)))
colnames(data.f)[12:ncol(data.f)] <- colnames(pxc_f_smooth)
colnames(data.m)[12:ncol(data.m)] <- colnames(pxc_m_smooth)

data.f[,(ncol(fem.smooth)+1):ncol(data.f)] <- data.f$mx * pxc_f_smooth
data.m[,(ncol(mal.smooth)+1):ncol(data.m)] <- data.m$mx * pxc_m_smooth

save(data.f, data.m, file = "R/final.RData")

dim(data.f)
dim(data.m)

#### Decompose ####

library(DecompHoriuchi)
library(MortHump)

source("R/functions/sdfun.R")
source("R/functions/MA5_FUN.R")

data <- as.vector(as.matrix(data.f[data.f$Year == 2001,12:ncol(data.f)]))

wrapper <- function(data){
  
  mxc <- matrix(data, ncol = 16, byrow = FALSE)
  
  mx <- rowSums(mxc, na.rm = TRUE)
  
  lt <- suppressWarnings(LT(Mx = mx, ages = xs, mxsmooth = FALSE, axmethod = "midpoint"))
  
  sd <- sdfun(x = data.frame(Age = xs, dx = lt$dx), trun = "mode")
  
  return(sd)

}


wrapper(data)

r1 <- as.vector(as.matrix(data.f[data.f$Year == 2001,12:ncol(data.f)]))
r2 <- as.vector(as.matrix(data.f[data.f$Year == 2004,12:ncol(data.f)]))

# about 30min!
dec <- stepwise_replacement(func = wrapper, rates1 = r1, rates2 = r2)

# sum of decomposition is equal to difference in sd => good!
sum(dec)
sdfun(x = data.f[data.f$Year == 2004,], trun = "mode") - sdfun(x = data.f[data.f$Year == 2001,], trun = "mode")

cxc <- matrix(dec, ncol = 16, byrow = FALSE)
colnames(cxc) <- colnames(pxc_f_smooth)

barplot(colSums(cxc), col = cols, las = 2)

#### all decompositions ####

# females
r1 <- as.vector(as.matrix(data.f[data.f$Year == 2001,12:ncol(data.f)]))
r2 <- as.vector(as.matrix(data.f[data.f$Year == 2004,12:ncol(data.f)]))
f_01_04_sd <- stepwise_replacement(func = wrapper, rates1 = r1, rates2 = r2)

r1 <- as.vector(as.matrix(data.f[data.f$Year == 2004,12:ncol(data.f)]))
r2 <- as.vector(as.matrix(data.f[data.f$Year == 2010,12:ncol(data.f)]))
f_04_10_sd <- stepwise_replacement(func = wrapper, rates1 = r1, rates2 = r2)

# males
r1 <- as.vector(as.matrix(data.m[data.m$Year == 2001,12:ncol(data.m)]))
r2 <- as.vector(as.matrix(data.m[data.m$Year == 2004,12:ncol(data.m)]))
m_01_04_sd <- stepwise_replacement(func = wrapper, rates1 = r1, rates2 = r2)

r1 <- as.vector(as.matrix(data.m[data.m$Year == 2004,12:ncol(data.m)]))
r2 <- as.vector(as.matrix(data.m[data.m$Year == 2010,12:ncol(data.m)]))
m_04_10_sd <- stepwise_replacement(func = wrapper, rates1 = r1, rates2 = r2)

save(f_01_04_sd, f_04_10_sd, m_01_04_sd, m_04_10_sd, file = "R/decomposition_sd.RData")


#### decomposition of e0 (stepwise and Arriaga => same) ####

wrapper_e0 <- function(data){
  
  mxc <- matrix(data, ncol = 16, byrow = FALSE)
  
  mx <- rowSums(mxc, na.rm = TRUE)
  
  lt <- suppressWarnings(LT(Mx = mx, ages = xs, mxsmooth = FALSE, axmethod = "midpoint"))
  
  e0 <- lt$ex[1]
  
  return(e0)
  
}

r1 <- as.vector(as.matrix(data.f[data.f$Year == 2001,12:ncol(data.f)]))
r2 <- as.vector(as.matrix(data.f[data.f$Year == 2004,12:ncol(data.f)]))
f_01_04 <- stepwise_replacement(func = wrapper_e0, rates1 = r1, rates2 = r2)


arriaga <- function(data1, data2, k = 10){
  
  mxc1 <- data1
  mx1 <- rowSums(mxc1, na.rm = TRUE)
  lt1 <- suppressWarnings(LT(Mx = mx1, ages = xs, mxsmooth = FALSE, axmethod = "midpoint"))
  lx1 <- lt1$lx
  Lx1 <- lt1$Lx
  
  mxc2 <- data2
  mx2 <- rowSums(mxc2, na.rm = TRUE)
  lt2 <- suppressWarnings(LT(Mx = mx2, ages = xs, mxsmooth = FALSE, axmethod = "midpoint"))
  lx2 <- lt2$lx
  Lx2 <- lt2$Lx

  dx <- vector(length = length(xs))
  
  for(i in 1:(length(xs) - 1)){
    
    dx[i] <- lx1[i]/lx1[1] * (Lx2[i]/lx2[i] - Lx1[i]/lx1[i]) + sum(Lx2[(i+1):length(Lx2)])/lx1[1] * (lx1[i]/lx2[i] - lx1[i+1]/lx2[i+1])
    
  }

  dc <- matrix(0, ncol = ncol(mxc1), nrow = nrow(mxc1))
  for(i in 1:nrow(dc)){
    
    dc[i,] <- unlist(dx[i] * (mxc2[i,] - mxc1[i,]) / (sum(mxc2[i,]) - sum(mxc1[i,])))
    
    # if the change in all-cause rate is much smaller than the cause-specific ones, it can lead to not meaningful results
    # if any of the latter is k time larger than the former, we drop this age alltogether
    test <- any(abs(mxc2[i,] - mxc1[i,]) > k * abs(sum(mxc2[i,]) - sum(mxc1[i,])))
    if(test == TRUE){dc[i,] <- 0}
  
  }
  
  rownames(dc) <- names(dx) <- xs
  colnames(dc) <- colnames(mxc1)
  
  return(dc)

}

dc <- arriaga(data.f[data.f$Year == 2001,12:ncol(data.f)], data.f[data.f$Year == 2004,12:ncol(data.f)])

# Arriaga vs. Andreev et al.
sum(dc) / sum(f_01_04)
barplot(colSums(dc), col = cols, las = 2)
barplot(colSums(matrix(f_01_04, ncol = 16, byrow = FALSE)), col = cols, las = 2)
plot(x = xs, rowSums(dc), type = "l")
lines(x = xs, rowSums(matrix(f_01_04, ncol = 16, byrow = FALSE)), col = 2)

# females
r1 <- data.f[data.f$Year == 2001,12:ncol(data.f)]
r2 <- data.f[data.f$Year == 2004,12:ncol(data.f)]
f_01_04_e0 <- arriaga(r1, r2)

r1 <- data.f[data.f$Year == 2004,12:ncol(data.f)]
r2 <- data.f[data.f$Year == 2010,12:ncol(data.f)]
f_04_10_e0 <- arriaga(r1, r2)

# males
r1 <- data.m[data.m$Year == 2001,12:ncol(data.m)]
r2 <- data.m[data.m$Year == 2004,12:ncol(data.m)]
m_01_04_e0 <- arriaga(r1, r2)

r1 <- data.m[data.m$Year == 2004,12:ncol(data.m)]
r2 <- data.m[data.m$Year == 2010,12:ncol(data.m)]
m_04_10_e0 <- arriaga(r1, r2)

save(f_01_04_e0, f_04_10_e0, m_01_04_e0, m_04_10_e0, file = "R/decomposition_e0.RData")

#### plots ####

load("R/decomposition_sd.RData")
load("R/decomposition_e0.RData")

f_01_04_sd <- matrix(f_01_04_sd, ncol = 16, byrow = FALSE)
colnames(f_01_04_sd) <- colnames(f_01_04_e0)
m_01_04_sd <- matrix(m_01_04_sd, ncol = 16, byrow = FALSE)
colnames(m_01_04_sd) <- colnames(m_01_04_e0)
f_04_10_sd <- matrix(f_04_10_sd, ncol = 16, byrow = FALSE)
colnames(f_04_10_sd) <- colnames(f_04_10_e0)
m_04_10_sd <- matrix(m_04_10_sd, ncol = 16, byrow = FALSE)
colnames(m_04_10_sd) <- colnames(m_04_10_e0)

png("Figures/decomp_e0.png", width = 20, height = 20, units = "cm", res = 600)

par(mar = c(7,4,4,2) + 0.1, mfrow = c(2,2))
barplot(colSums(f_01_04_e0), col = cols, las = 2, main = "Females 2001-04, e0")
barplot(colSums(m_01_04_e0), col = cols, las = 2, main = "Males 2001-04, e0")
barplot(colSums(f_04_10_e0), col = cols, las = 2, main = "Females 2004-10, e0")
barplot(colSums(m_04_10_e0), col = cols, las = 2, main = "Males 2004-10, e0")

dev.off()


png("Figures/decomp_sd.png", width = 20, height = 20, units = "cm", res = 600)

par(mar = c(7,4,4,2) + 0.1, mfrow = c(2,2))
barplot(colSums(f_01_04_sd), col = cols, las = 2, main = "Females 2001-04, sd+")
barplot(colSums(m_01_04_sd), col = cols, las = 2, main = "Males 2001-04, sd+")
barplot(colSums(f_04_10_sd), col = cols, las = 2, main = "Females 2004-10, sd+")
barplot(colSums(m_04_10_sd), col = cols, las = 2, main = "Males 2004-10, sd+")

dev.off()

