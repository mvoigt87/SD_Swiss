# Redistribution of ill-defined causes of death
# 2000 - 2014
# Adrien Remund, April 2018

#### load data ####

if(Sys.info()[['nodename']] == "CD-5VV9QK2"){
  
  setwd("D:/switchdrive/SNC 2.0/TOMCOD")
  
}else{ if(Sys.info()[['nodename']] == "REMUND"){
  
  setwd("C:/Users/Adrien/switchdrive/SNC 2.0/TOMCOD")
  
}else{
  
  setwd("C:/Users/Rose van der Linden/switchdrive/SNC 2.0/TOMCOD")
  
}}

load(file = "R/IDCD/bycanton.RData")

# reorganize from cantons to sexes
males <- as.data.frame(do.call(rbind,lapply(pc, function(x){x$Male})))
females <- as.data.frame(do.call(rbind,lapply(pc, function(x){x$Female})))

library(RColorBrewer)
cols <- c(brewer.pal(n = 9, name = "Set1"),brewer.pal(n = 8, name = "Set3"))

#### marginal distributions ####

# distribution by cause
png(filename = "Figures/pc.png", width = 40, height = 20, units = "cm", res = 600)
par(mfrow = c(1,2), mar = c(6,4,4,2) + 0.1)

tab <- prop.table(table(data$cod[data$sex == "Male"]))
rank <- order(tab, decreasing = TRUE)
barplot(tab[rank], col = cols[rank], main = "Males", las = 2, ylab = "% of all deaths", cex.names = 0.8)
tab <- prop.table(table(data$cod[data$sex == "Female"]))
rank <- order(tab, decreasing = TRUE)
barplot(tab[rank], col = cols[rank], main = "Females", las = 2, ylab = "% of all deaths", cex.names = 0.8)

par(mfrow = c(1,1), mar = c(5,4,4,2) + 0.1)
dev.off()

# distribution by canton
png(filename = "Figures/pc_canton.png", width = 30, height = 12, units = "cm", res = 600)
par(mar = c(5,4,4,1))
layout(mat = rbind(c(1,1,3,2,2)))

rank <- order(prcomp(males)$x[,1])
barplot(t(males[rank,]), beside = FALSE, las = 2, cex.names = 0.8, main = "Males", col = cols)
rank <- order(prcomp(females)$x[,1])
barplot(t(females[rank,]), beside = FALSE, las = 2, cex.names = 0.8, main = "Females", col = cols)
plot(0,0, type = "n", axes = FALSE, xlab = "", ylab = "")
legend("topleft", legend = sort(unique(data$cod), decreasing = TRUE), fill = rev(cols), ncol = 1, cex = 1.2)

par(mar = c(5,4,4,2) + 0.1)
dev.off()

# PCA on the distribution by canton to look for different epidemiological regimes within Switzerland
# looks like there is a "Romandie effect" 
cols2 <- unique(data$canton) %in% c("GE","VD","VS","NE","JU","FR") + 1
  
png(filename = "Figures/pca.png", width = 30, height = 15, units = "cm", res = 600)
par(mfrow = c(1,2))
pca <- prcomp(males)
x <- pca$x[,1]
x <- (x - mean(x)) / sd(x)
y <- pca$x[,2]
y <- (y - mean(y)) / sd(y)
plot(x, y, las = 1, ylim = range(y) * 1.1, xlim = range(x) * 1.1, main = "Males", 
     xlab = "", ylab = "", pch = 16, col = cols2)
text(x = x, y = y, labels = rownames(males), cex = 0.7, pos = 3)

pca <- prcomp(females)
x <- pca$x[,1]
x <- (x - mean(x)) / sd(x)
y <- pca$x[,2]
y <- (y - mean(y)) / sd(y)
plot(x, y, las = 1, ylim = range(y) * 1.1, xlim = range(x) * 1.1, main = "Females", 
     xlab = "", ylab = "", pch = 16, col = cols2)
text(x = x, y = y, labels = rownames(males), cex = 0.7, pos = 3)

par(mfrow = c(1,1))
dev.off()

#### Ill-defined ####

# share of ill-defined
# double effect of French-speaking and city-state
# Geneva combines both and is clearly an outlier (esp. for males)
png("Figures/illdefined_all.png", width = 40, height = 20, units = "cm", res = 600)
par(mfrow = c(1,2), mar = c(4,4,3,1) + 0.1)

barplot(sort(males$`Ill-defined`), names.arg = rownames(males)[order(males$`Ill-defined`)], las = 2, cex.names = 0.7, ylab = "% ill-defined", main = "Males", ylim = c(0,0.12))
barplot(sort(females$`Ill-defined`), names.arg = rownames(females)[order(females$`Ill-defined`)], las = 2, cex.names = 0.7, ylab = "% ill-defined", main = "Females", ylim = c(0,0.12))

par(mfrow = c(1,1), mar = c(5,4,4,1) + 0.1)
dev.off()

# linear models
# Heart diseases are the only cause that is strongly negatively associated with ill-defined causes
# Circulatory diseases are also playing a role, but only with males
# The "other" category seems strongly positively associated, so there might be garbage codes in there as well => find them and put them in ill-defined
# Geneva is clearly having a strong influence on the results, we should run sensitivity analyses without it

png("Figures/lm_GE.png", width = 40, height = 90, units = "cm", res = 600)
par(mfcol = c(8,4), mar = c(2,2,4,1))

ge <- which(rownames(males) == "GE")

xx <- seq(0,0.14,0.01)

coefs.males <- matrix(NA, ncol = 2, nrow = ncol(males)) ; coefs.males[10,] <- c(0,NA)
for(i in c(1:9,11:ncol(males))){
  
  fit <- lm(males[-ge,i] ~ males[-ge,which(names(males) == "Ill-defined")])
  coefs.males[i,] <- summary(fit)$coef[2,c(1,4)]
  plot(males[,which(names(males) == "Ill-defined")], males[,i], las = 1, main = names(males)[i], xlab = "", ylab = "", ylim = c(0,0.4))
  text(males[,which(names(males) == "Ill-defined")], males[,i], rownames(males), cex = 0.6, pos = 3)
  lines(xx, coef(fit)[1] + xx * coef(fit)[2], lwd = 2, col = c(8,2)[(summary(fit)$coef[2,4] < 0.05) + 1])
  
}

coefs.females <- matrix(NA, ncol = 2, nrow = ncol(males)) ; coefs.females[10,] <- c(0,NA)
for(i in c(1:9,11:ncol(females))){
  
  fit <- lm(females[-ge,i] ~ females[-ge,which(names(females) == "Ill-defined")])
  coefs.females[i,] <- summary(fit)$coef[2,c(1,4)]
  plot(females[,which(names(females) == "Ill-defined")], females[,i], las = 1, main = names(females)[i], xlab = "", ylab = "", ylim = c(0,0.4))
  text(females[,which(names(females) == "Ill-defined")], females[,i], rownames(females), cex = 0.6, pos = 3)
  lines(xx, coef(fit)[1] + xx * coef(fit)[2], lwd = 2, col = c(8,2)[(summary(fit)$coef[2,4] < 0.05) + 1])
  
}

par(mfcol = c(1,1),mar = c(5,4,4,2) + 0.1)
dev.off()

# looking at the coefficients
# Heart diseases are clearly very correlated
# cerebrovascular, Circulatory and Endocrinian disease too, but to a smaller extent
# The other associations are probably spurious, possibly due to a "Geneva-bias"

png("Figures/coefs_GE.png", width = 20, height = 10, units = "cm", res = 600)

par(mfrow = c(1,2))

mp <- barplot(coefs.males[,1], col = cols, names.arg = names(males), las = 2, cex.names = 0.7, main = "Males", ylim = range(coefs.males[,1]) * 1.2)
text(x = mp, y = coefs.males[,1], labels = c("","*")[(coefs.males[,2] < 0.05) + 1], pos = c(1,3)[(coefs.males[,1] > 0) + 1])
barplot(coefs.females[,1], col = cols, names.arg = names(males), las = 2, cex.names = 0.7, main = "Females", ylim = range(coefs.females[,1]) * 1.2)
text(x = mp, y = coefs.females[,1], labels = c("","*")[(coefs.females[,2] < 0.05) + 1], pos = c(1,3)[(coefs.females[,1] > 0) + 1])

par(mfrow = c(1,1))
dev.off()

#### redistribution coefficients ####

coefs <- data.frame(males = -coefs.males[,1], females = -coefs.females[,1])
rownames(coefs) <- names(males)

coefs[-c(3,4,6,9),] <- 0 # keep only the most pertinent effects

coefs <- t( t(coefs)/ rowSums(t(coefs)))

png("Figures/redistr_GE.png", width = 20, height = 10, units = "cm", res = 600)

par(mfrow = c(1,2))

barplot(coefs[,1], names.arg = names(males), col = cols, las = 2, cex.names = 0.7, main = "Males")

barplot(coefs[,2], names.arg = names(males), col = cols, las = 2, cex.names = 0.7, main = "Females")

par(mfrow = c(1,1))
dev.off()

save(coefs, file = "R/IDCD/redistribute.RData")

#### todo ####

# 1.	Rerun without Geneva
# 2.	Rerun without infant deaths
# 3.	Rerun after age-standardization
# 4.	Try splitting further heart diseases


