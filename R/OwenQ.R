#' @title First Owen Q-function
#' @description Evaluates the first Owen Q-function (integral from \eqn{0} to \eqn{R})
#' for an integer value of the degrees of freedom.
#' @param nu integer greater than \eqn{1}, the number of degrees of freedom
#' @param t finite number, positive or negative
#' @param delta vector of finite numbers, with the same length as \code{R}
#' @param R (upper bound of the integral) vector of finite positive numbers, with the same length as \code{delta}
#' @return A vector of numbers between \eqn{0} and \eqn{1}, the values of the integral from \eqn{0} to \eqn{R}.
#' @export
#' @importFrom stats pnorm dnorm
#' @note The results are theoretically exact when the number of degrees of freedom is even.
#' When odd, the procedure resorts to the Owen T-function.
#' @examples
#' # OwenQ1(nu, t, delta, Inf) = pt(t, nu, delta)
#' OwenQ1(4, 3, 2, 100)
#' ptOwen(3, 4, 2)
OwenQ1 <- function(nu, t, delta, R){
  J <- length(delta)
  if(J != length(R)){
    stop("`delta` and `R` must have the same length.")
  }
  if(any(R<0)){
    stop("R must be positive.")
  }
  if(isNotPositiveInteger(nu)){
    stop("`nu` must be an integer >=1.")
  }
  if(is.infinite(t) || any(is.infinite(delta)) || any(is.infinite(R))){
    stop("Parameters must be finite.")
  }
  a <- sign(t)*sqrt(t*t/nu)
  b <- nu/(nu+t*t)
  sB <- sqrt(b)
  if(nu==1){
    C <- pnorm(R) - (delta>=0) + 2*OwenT(delta*sB, a) -
      vapply(seq_len(J), function(i){
        2*OwenT(R[i], (a*R[i]-delta[i])/R[i]) +
          2*OwenT(delta[i]*sB, (delta[i]*a*b-R[i])/b/delta[i])
      }, FUN.VALUE=numeric(1L))
    return(C)
  }
  nu <- as.integer(nu)
  n <- nu-1L
  H <- M <- matrix(numeric(n*J), nrow=n)
  H[1L,] <- -dnorm(R) * pnorm (a*R-delta)
  M[1L,] <- a*sB*dnorm(delta*sB)*(pnorm(delta*a*sB)-pnorm((delta*a*b-R)/sB))
  if(nu>=3L){
    H[2L,] <- R * H[1L,]
    M[2L,] <- b*(delta*a*M[1L,] +
                   a*dnorm(delta*sB)*(dnorm(delta*a*sB)-dnorm((delta*a*b-R)/sB)))
    if(nu>=4L){
      A <- numeric(n)
      L <- matrix(numeric((n-2L)*J), ncol=J)
      A[1L:2L] <- 1
      L[1L,] <- a * b * R * dnorm(R) * dnorm(a*R-delta) / 2
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
        M[k,] <- (k-2L)/(k-1L) * b * (A[k-2L] * delta * a * M[k-1L,] + M[k-2L,]) -
          L[k-2L,]
      }
    }
  }
  if(nu%%2L==1L){
    C <- pnorm(R) - (delta>=0) + 2*OwenT(delta*sB, a) -
      vapply(seq_len(J), function(i){
        2*OwenT(R[i], (a*R[i]-delta[i])/R[i]) +
          2*OwenT(delta[i]*sB, (delta[i]*a*b-R[i])/b/delta[i])
      }, FUN.VALUE=numeric(1L))
    indices <- seq(2L, n, by=2L)
    return(C + 2*.colSums(H[indices,]+M[indices,], m=length(indices), n=J))
  }else{
    indices <- seq(1L, n, by=2L)
    return(pnorm(-delta) + sqrt(2*pi) *
             .colSums(H[indices,]+M[indices,], m=length(indices), n=J))
  }
}

