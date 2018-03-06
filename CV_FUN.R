CV.FUN <-
function(x){
  
  mean <- sum(x$Age * x$dx) / sum(x$dx)
  
  sdev <- sqrt(sum((x$Age - mean)^2 * x$dx) / sum(x$dx))
  
  # coefficient of variance - CV
  cv <- ((sdev)/mean) * 100
  
  return(cv)
  
}
