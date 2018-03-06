IQR.FUN <-
function(x){
    ### IQR function:
  x$rel.DX.CUM <- cumsum(x$dx)/max(cumsum(x$dx))
  Q1 <- x$Age[min(which(x$rel.DX.CUM>0.25))]
  Q3 <- x$Age[min(which(x$rel.DX.CUM>0.75))]
  IQR <- Q3 - Q1
  return(IQR)
}
