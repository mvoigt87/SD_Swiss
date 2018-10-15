
####################################################
# computing standard deviation in the age at death #
#       including an attempt at smoothing and      #
#        interpolation of sub-integer ages         #
#              TOMCOD project                      #
#           A.Remund February 2018                 #
####################################################

#### set working directory ####
if(Sys.info()["nodename"] == "CD-5VV9QK2"){
  
  setwd("D:/switchdrive/SNC 2.0/TOMCOD")
  
}

#### load data ####
library(HMDHFDplus)

data <- readHMDweb(CNTRY = "CHE", item = "mltper_1x1", username = "", password = "")

data <- data[,c("Year","Age","dx")]

#### interpolate ####
# needs to be fixed

mat <- matrix(data$dx, nrow = length(unique(data$Age)))

library(MortalitySmooth)

fit <- Mort2Dsmooth(x = unique(data$Age), y = unique(data$Year), Z = mat, method = 3, lambdas = c(1,1))

plot(fit)

xx <- seq(0,max(data$Age),0.1)

mat <- predict(fit, type = "response", newdata = list(y = unique(data$Year), x = xx)) / 10 # something wrong with interpolation => missing values created

data_interpol <- data.frame(Year = rep(unique(data$Year), each = nrow(mat)), Age = rep(xx, times = ncol(mat)))
data_interpol$dx <- as.vector(mat)
matplot(mat, type = "l")

#### function ####

sdfun <- function(x, smooth = FALSE, plot = FALSE, trun = 0, inter = FALSE){

  if(plot == TRUE){
    plot(x$Age, x$dx, las = 1, xlab = "age", ylab = "dx")
  }
  
  if(is.numeric(trun)){
    
    x <- x[x$Age >= trun,]
    
    if(plot == TRUE){
      
      abline(v = trun)
      
      }
  
    }
  
  if(smooth == TRUE){
    
    if(is.numeric(inter)){
      
      dx <- predict(smooth.spline(x = x$Age, y = x$dx, control.spar = list(low = 0.3, high = 1)), x = inter)$y
      
      x <- data.frame(Year = rep(x$Year,length(inter)), Age = inter)
      
      x$dx <- dx #* c(diff(inter),1)
      
    }else{
    
    x$dx <- predict(smooth.spline(x = x$Age, y = x$dx, control.spar = list(low = 0.3, high = 1)))$y
    
    }
    
    if(plot == TRUE){
      lines(x$Age, x$dx)
    }
    
  }
  
  if(trun == "mode"){
    
    ranks <- order(x$dx, decreasing = TRUE)
    
    i <- 1
    idx <- ranks[i]
    mode <- x$Age[idx]
    while(mode < 5){
      i <- i + 1
      idx <- ranks[i]
      mode <- x$Age[idx]
      }
    
    x <- x[x$Age >= mode,]
    
    if(plot == TRUE){
      
      abline(v = mode)
      
    }
    
    }
  
  mean <- sum(x$Age * x$dx) / sum(x$dx)
  
  sdev <- sqrt(sum((x$Age - mean)^2 * x$dx) / sum(x$dx))
  
  return(sdev)
  
}

#### apply ####

sdfun(data[data$Year == 1942,], smooth = TRUE, plot = TRUE, trun = "mode", inter = seq(0,110,0.1))

sd <- by(data = data, INDICES = data$Year, FUN = sdfun, trun = trun, smooth = TRUE, inter = seq(0,110,0.1))

plot(unique(data$Year), sd, type = "l", xlab = "", ylab = "sd", las = 1)


#### comparison ####

pdf("Output/comparison.pdf", width = 10, height = 5)

par(mfrow = c(1,2))

for(trun in c(0,"mode")){

plot(1, 1, type = "n", ylim = c(ifelse(trun == 0, 10, 3),ifelse(trun == 0, 35, 6)), xlim = c(1900, 2020), xlab = "", ylab = "sd", las = 1, main = paste("truncated at",trun))

for(ctry in c("CHE","NLD")){
  
  cat("\n",ctry)
  
  for(sex in c("m","f")){
    
    cat("\n","-",sex)
    
    data <- readHMDweb(CNTRY = ctry, item = paste(sex,"ltper_1x1",sep = ""), username = "", password = "")
    
    data <- data[,c("Year","Age","dx")]
    
    sd <- by(data = data, INDICES = data$Year, FUN = sdfun, trun = trun, smooth = TRUE, inter = c(seq(0,110,0.1)))
      
    col <- ifelse(sex == "m", "blue", "red")
    
    lty <- ifelse(ctry == "CHE", 1, 2)
      
    lines(unique(data$Year), sd, col = col, lty = lty)

      }
    
  }
  
  legend("topright", legend = c("CH females", "CH males", "NL females", "NL males"), lty = c(1,1,2,2), col = c("red","blue","red","blue"))
  
  }

par(mfrow = c(1,1))

dev.off()


