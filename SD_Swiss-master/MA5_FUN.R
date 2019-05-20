MDA.fun <-
function(lt) {
  x <- lt$Age[which.max(lt$dx)]
  M <- x + ((lt$dx[x]-lt$dx[x-1])/(lt$dx[x]-lt$dx[x-1])+(lt$dx[x]-lt$dx[x+1]))
  return(M)
}
