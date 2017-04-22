# OwenSequences <- function(n, a=1, b=1, d=1, r=1){
#   H <- M <- numeric(n)
#   H[1L] <- -dnorm(r) * pnorm (a*r-d)
#   sB <- sqrt(b)
#   M[1L] <- a*sB*dnorm(d*sB)*(pnorm(d*a*sB)-pnorm((d*a*b-r)/sB))
#   # if n>1
#   H[2L] <- r * H[1L]
#   M[2L] <- b*(d*a*M[1L] + a*dnorm(d*sB)*(dnorm(d*a*sB)-dnorm((d*a*b-r)/sB)))
#   if(n>2){
#     A <- numeric(n)
#     L <- numeric(n-2L)
#     A[1L:2L] <- 1
#     L[1L] <- a * b * r * dnorm(r) * dnorm(a*r-d) / 2
#     for(k in 3L:n){
#       A[k] <- 1/(k-1L)/A[k-1L]
#     }
#     if(n>3){
#       for(k in 2L:(n-2L)){
#         L[k] <- A[k+2L] * r * L[k-1L]
#       }
#     }
#     for(k in 3L:n){
#       H[k] <- A[k] * r * H[k-1L]
#       M[k] <- (k-2L)/(k-1L) * b * (A[k-2L] * d * a * M[k-1L] + M[k-2L]) - L[k-2L]
#     }
#   }
#   return(cbind(H,M))
# }


#' @title First Owen Q-function
#' @description Evaluates the first Owen Q-function (integral from \eqn{0} to \eqn{R})
#' for an integer value of the degrees of freedom.
#' @param nu integer greater than \eqn{1}, the number of degrees of freedom
#' @param t finite number, positive or negative
#' @param delta finite number, positive or negative
#' @param R finite positive number, the upper bound of the integral
#' @return A number between \eqn{0} and \eqn{1}, the value of the integral from \eqn{0} to \eqn{R}.
#' @export
#' @examples
#' # OwenQ1(nu, t, delta, Inf) = pt(t, nu, delta)
#' OwenQ1(4, 3, 2, 100)
#' ptOwen(3, 4, 2)
OwenQ1 <- function(nu, t, delta, R){
  if(R<0){
    stop("R must be positive.")
  }
  if(isNotPositiveInteger(nu)){
    stop("`nu` must be an integer >=1.")
  }
  if(is.infinite(t) || is.infinite(delta) || is.infinite(R)){
    stop("Parameters must be finite.")
  }
  a <- t/sqrt(nu)
  b <- nu/(nu+t*t)
  sB <- sqrt(b)
  if(nu==1){
    C <- pnorm(R) - 2*OwenT(R, (a*R-delta)/R) -
      2*OwenT(delta*sB, (delta*a*b-R)/b/delta) + 2*OwenT(delta*sB, a) -
      (delta>=0)
    return(C)
  }
  nu <- as.integer(nu)
  n <- nu-1L
  H <- M <- numeric(n)
  H[1L] <- -dnorm(R) * pnorm (a*R-delta)
  M[1L] <- a*sB*dnorm(delta*sB)*(pnorm(delta*a*sB)-pnorm((delta*a*b-R)/sB))
  if(nu>=3L){
    H[2L] <- R * H[1L]
    M[2L] <- b*(delta*a*M[1L] + a*dnorm(delta*sB)*(dnorm(delta*a*sB)-dnorm((delta*a*b-R)/sB)))
    if(nu>=4L){
      A <- numeric(n)
      L <- numeric(n-2L)
      A[1L:2L] <- 1
      L[1L] <- a * b * R * dnorm(R) * dnorm(a*R-delta) / 2
      for(k in 3L:n){
        A[k] <- 1/(k-1L)/A[k-1L]
      }
      if(nu>=5L){
        for(k in 2L:(n-2L)){
          L[k] <- A[k+2L] * R * L[k-1L]
        }
      }
      for(k in 3L:n){
        H[k] <- A[k] * R * H[k-1L]
        M[k] <- (k-2L)/(k-1L) * b * (A[k-2L] * delta * a * M[k-1L] + M[k-2L]) - L[k-2L]
      }
    }
  }
  if(nu%%2L==1L){
    C <- pnorm(R) - 2*OwenT(R, (a*R-delta)/R) -
      2*OwenT(delta*sB, (delta*a*b-R)/b/delta) + 2*OwenT(delta*sB, a) -
      (delta>=0)
    indices <- seq(2L, n, by=2L)
    return(C + 2*sum(H[indices]+M[indices]))
  }else{
    indices <- seq(1L, n, by=2L)
    return(pnorm(-delta) + sqrt(2*pi) * sum(H[indices]+M[indices]))
  }
}

