
LT <- function (Mx, Age = c(0, 1, cumsum(rep(5, length(Mx) - 2))), radix = 1e+05){

  N <- length(Mx)
  w <- c(diff(Age), 5)
  ax <- c(0.07 + 1.7 * Mx[1], 0.4, rep(0.5, N - 3), 1/Mx[N])
  Nax <- w * ax
  qx <- (w * Mx)/(1 + w * (1 - ax) * Mx)
  qx[N] <- 1
  px <- 1 - qx
  lx <- c(radix, radix * (cumprod(px[-N])))
  dx <- -diff(c(lx, 0))
  Lx <- lx[-1] * w[-N] + dx[-N] * Nax[-N]
  Lx[N] <- ax[N] * lx[N]
  Tx <- rev(cumsum(rev(Lx)))
  ex <- Tx/lx
  return(lt = list(Age = Age, ax = ax, mx = Mx, qx = qx, px = px, lx = lx, dx = dx, Lx = Lx, Tx = Tx, ex = ex))
}