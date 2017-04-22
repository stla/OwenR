Ssequences <- function(n, a, b, d){
  A <- M <- numeric(n)
  sb <- sqrt(b)
  M[1L] <- a * sb * dnorm(d*sb) * pnorm(d*a*sb)
  A[2L] <- 1
  M[2L] <- b * (d * a * M[1L] + a * dnorm(d) / sqrt(2*pi))
  for(k in 3L:n){
    A[k] <- 1 /(k-2)/A[k-1L]
    M[k] <- (k-2)/(k-1) * b * (A[k-1L] * d * a * M[k-1L] + M[k-2L])
  }
  return(M)
}

#' @title Student CDF with integer number of degrees of freedom
#' @description Cumulative distribution function of the noncentrel Student
#' distribution with an integer number of degrees of freedom.
#' @param q quantile
#' @param nu integer greater than \eqn{1}, the number of degrees of freedom
#' @param delta noncentrality parameter
#' @return Numeric value, the CDF evaluated at \code{q}.
#' @export
#' @examples
#' ptOwen(2, 3) - pt(2, 3)
#' ptOwen(2, 3, delta=1) - pt(2, 3, ncp=1)
ptOwen <- function(q, nu, delta=0){
  if(isNotPositiveInteger(nu)){
    stop("`nu` must be an integer >=1.")
  }
  if(is.infinite(q) || is.infinite(delta)){
    stop("Parameters must be finite.")
  }
  a <- q/sqrt(nu)
  b <- nu/(nu+q*q)
  nu <- as.integer(nu)
  if(nu==2L){
    sB <- sqrt(b)
    return(pnorm(-delta) + sqrt(2*pi) *
             (a * sB * dnorm(delta*sB) * pnorm(delta*a*sB)))
  }
  if(nu%%2L==1L){
    sB <- sqrt(b)
    C <- pnorm(-delta*sB) + 2*OwenT(delta*sB,a)
    if(nu==1L){
      return(C)
    }
    if(nu==3L){
      return(C + 2 * b * (delta * a * a * sB * dnorm(delta*sB) * pnorm(delta*a*sB)
                     + a * dnorm(delta) / sqrt(2*pi)))
    }
    return(C + 2*sum(Ssequences(nu-1, a, b, delta)[seq(2L, nu-1L, by=2L)]))
  }
  return(pnorm(-delta) + sqrt(2*pi) * sum(Ssequences(nu-1, a, b, delta)[seq(1L, nu-1L, by=2L)]))
}

