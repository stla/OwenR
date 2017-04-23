#' @title First Owen Q-function
#' @description Evaluates the first Owen Q-function (integral from \eqn{0} to \eqn{R})
#' for an integer value of the degrees of freedom.
#' @param nu integer greater than \eqn{1}, the number of degrees of freedom
#' @param t finite number, positive or negative
#' @param delta vector of finite numbers, with the same length as \code{R}
#' @param R (upper bound of the integral) vector of finite positive numbers, with the same length as \code{delta}
#' @param jmax,cut.point passed to \code{\link{OwenT}} (when \code{nu} is odd)
#' @return A vector of numbers between \eqn{0} and \eqn{1}, the values of the integral from \eqn{0} to \eqn{R}.
#' @export
#' @importFrom stats pnorm dnorm
#' @note The results are theoretically exact when the number of degrees of freedom is even.
#' When odd, the procedure resorts to the Owen T-function.
#' @examples
#' # OwenQ1(nu, t, delta, Inf) = pt(t, nu, delta)
#' OwenQ1(nu=4, t=3, delta=2, R=100)
#' ptOwen(q=3, nu=4, delta=2)
OwenQ1 <- function(nu, t, delta, R, jmax=50L, cut.point=8){
  J <- length(delta)
  if(J != length(R)){
    stop("`delta` and `R` must have the same length.")
  }
  if(any(R<0)){
    stop("`R` must be positive.")
  }
  if(isNotPositiveInteger(nu)){
    stop("`nu` must be an integer >=1.")
  }
  if(any(is.infinite(R))){
    stop("`R` must be finite.")
  }
  if(any(is.infinite(delta))){
    stop("`delta` must be finite.")
  }
  a <- sign(t)*sqrt(t*t/nu)
  b <- nu/(nu+t*t)
  sB <- sqrt(b)
  ab <- ifelse(is.infinite(t), 0, a*b)
  asB <- ifelse(is.infinite(t), sign(t), sign(t)*sqrt(t*t/(nu+t*t)))
  if(nu==1){
    C <- pnorm(R) - (delta>=0) + 2*.OwenT(delta*sB, a, jmax=jmax, cut.point=cut.point) -
      vapply(seq_len(J), function(i){
        2*.OwenT(R[i], (a*R[i]-delta[i])/R[i], jmax=jmax, cut.point=cut.point) +
          2*.OwenT(delta[i]*sB, (delta[i]*ab-R[i])/b/delta[i], jmax=jmax, cut.point=cut.point)
      }, FUN.VALUE=numeric(1L))
    return(C)
  }
  nu <- as.integer(nu)
  n <- nu-1L
  H <- M <- matrix(numeric(n*J), nrow=n)
  H[1L,] <- -dnorm(R) * pnorm (a*R-delta)
  M[1L,] <- asB*dnorm(delta*sB)*(pnorm(delta*asB)-pnorm((delta*ab-R)/sB))
  if(nu>=3L){
    H[2L,] <- R * H[1L,]
    M[2L,] <- ab*(delta*M[1L,] +
                   dnorm(delta*sB)*(dnorm(delta*asB)-dnorm((delta*ab-R)/sB)))
    if(nu>=4L){
      A <- numeric(n)
      L <- matrix(numeric((n-2L)*J), ncol=J)
      A[1L:2L] <- 1
      L[1L,] <- ab * R * dnorm(R) * dnorm(a*R-delta) / 2
      for(k in 3L:n){
        A[k] <- 1/(k-1L)/A[k-1L]
      }
      if(nu>=5L){
        for(k in 2L:(n-2L)){
          L[k,] <- A[k+2L] * R * L[k-1L,]
        }
      }
      for(k in 3L:n){
        H[k,] <- A[k] * R * H[k-1L,]
        M[k,] <- (k-2L)/(k-1L) * (ab * A[k-2L] * delta * M[k-1L,] + b*M[k-2L,]) -
          L[k-2L,]
      }
    }
  }
  if(nu%%2L==1L){
    C <- pnorm(R) - (delta>=0) + 2*.OwenT(delta*sB, a, jmax=jmax, cut.point=cut.point) -
      vapply(seq_len(J), function(i){
        2*.OwenT(R[i], (a*R[i]-delta[i])/R[i], jmax=jmax, cut.point=cut.point) +
          2*.OwenT(delta[i]*sB, (delta[i]*ab-R[i])/b/delta[i], jmax=jmax, cut.point=cut.point)
      }, FUN.VALUE=numeric(1L))
    indices <- seq(2L, n, by=2L)
    return(C + 2*.colSums(H[indices,]+M[indices,], m=length(indices), n=J))
  }else{
    indices <- seq(1L, n, by=2L)
    return(pnorm(-delta) + sqrt(2*pi) *
             .colSums(H[indices,]+M[indices,], m=length(indices), n=J))
  }
}

#' @useDynLib OwenR
OwenQ1_Rcpp <- function(nu, t, delta, R, jmax=50L, cut.point=8){
  J <- length(delta)
  if(J != length(R)){
    stop("`delta` and `R` must have the same length.")
  }
  if(any(R<0)){
    stop("`R` must be positive.")
  }
  if(isNotPositiveInteger(nu)){
    stop("`nu` must be an integer >=1.")
  }
  if(any(is.infinite(R))){
    stop("`R` must be finite.")
  }
  if(any(is.infinite(delta))){
    stop("`delta` must be finite.")
  }
  RcppOwenQ1(nu, t, delta, R,
             function(h,a) .OwenT(h, a, jmax=jmax, cut.point=cut.point))
}
