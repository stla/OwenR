T.int <- function(h, a, jmax, cut.point) {
  # fui <- function(h, i)
  #   #(h ^ (2 * i)) / ((2 ^ i) * gamma(i + 1))
  #   exp(2*i*log(h) - i*log(2) - lgamma(i+1))
  seriesL <- seriesH <-  NULL
  i <- 0L:jmax
  low <- (h <= cut.point)
  hL <- h[low]
  hH <- h[!low]
  L <- length(hL)
  if (L > 0L) {
    b <- cbind(1, outer(hL, 1L:jmax, function(h,i) exp(2*i*log(h) - i*log(2) - lgamma(i+1L))))
    cumb <- apply(b, 1L, cumsum)
    b1 <- exp(-0.5 * hL ^ 2) * t(cumb)
    matr <- matrix(1, jmax + 1L, L) - t(b1)
    jk <- rep(c(1L,-1L), jmax)[1L:(jmax + 1L)] / (2 * i + 1)
    matr <- t(matr * jk) %*% a ^ (2 * i + 1)
    seriesL <- (atan(a) - as.vector(matr)) / (2 * pi)
  }
  if (length(hH))
    seriesH <- atan(a) * exp(-0.5 * (hH ^ 2) * a / atan(a)) *
      (1 + 0.00868 * (hH * a) ^ 4) / (2 * pi)
  series <- c(seriesL, seriesH)
  id <- c((1:length(h))[low], (1:length(h))[!low])
  series[id] <- series
  series
}


#' @title Owen T-function
#' @description Evaluates the Owen T-function.
#' @param h numeric vector
#' @param a numeric scalar
#' @param jmax integer scalar which regulates the accuracy of the result
#' @param cut.point scalar number which regulates the behaviour of the algorithm
#' @return A vector of numbers between 0 and 1.
#' @export
#' @importFrom stats pnorm
#' @examples
# # OwenT(h,a) = OwenT(-h,a)
#' OwenT(2,1) == OwenT(-2,1)
# # OwenT(0,a) = atan(a)/2pi
#' a <- runif(1, -1000, 1000)
#' OwenT(0,a) - atan(a)/(2*pi)
# # OwenT(h,1) = Phi(h)(1-Phi(h))/2
#' h <- runif(1, -3, 3)
#' OwenT(h,1) - pnorm(h)*(1-pnorm(h))/2
# # OwenT(h,Inf) = (1-Phi(|h|))/2 :
#' OwenT(1,10000) - (1-pnorm(abs(1)))/2
#' OwenT(1,Inf) == (1-pnorm(abs(1)))/2
#' @export
OwenT <- function (h,
                   a,
                   jmax = 50L,
                   cut.point = 8)
{
  if (!is.vector(a) | length(a) > 1)
    stop("'a' must be a vector of length 1")
  if (!is.vector(h))
    stop("'h' must be a vector")
  aa <- abs(a)
  ah <- abs(h)
  if (is.na(aa))
    stop("parameter 'a' is NA")
  if (aa == Inf)
    return(sign(a) * 0.5 * pnorm(-ah))
  if (aa == 0)
    return(rep(0, length(h)))
  na <- is.na(h)
  inf <- (ah == Inf)
  ah <- replace(ah, (na | inf), 0)
  if (aa <= 1)
    owen <- T.int(ah, aa, jmax, cut.point)
  else
    owen <- (0.5 * pnorm(ah) + pnorm(aa * ah) * (0.5 - pnorm(ah)) -
               T.int(aa * ah, (1 / aa), jmax, cut.point))
  owen <- replace(owen, na, NA)
  owen <- replace(owen, inf, 0)
  return(owen * sign(a))
}
