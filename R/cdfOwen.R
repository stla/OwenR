#' @title Owen's equality 11
#' @description Evaluates the Owen cumulative distribution function in the 4th case.
#' @export
pOwen4 <- function(nu, t1, t2, delta1, delta2, R){
  J <- length(R)
  if(J != length(delta1) || J != length(delta2)){
    stop("`delta1`, `delta2` and `R` must have the same length.")
  }
  if(any(R<0)){
    stop("R must be positive.")
  }
  if(isNotPositiveInteger(nu)){
    stop("`nu` must be an integer >=1.")
  }
  if(any(is.infinite(R))){
    stop("R must be finite.")
  }
  a1 <- sign(t1)*sqrt(t1*t1/nu)
  b1 <- nu/(nu+t1*t1)
  sB1 <- sqrt(b1)
  ab1 <- ifelse(is.infinite(t1), 0, a1*b1)
  asB1 <- ifelse(is.infinite(t1), sign(t1), sign(t1)*sqrt(t1*t1/(nu+t1*t1)))
  a2 <- sign(t2)*sqrt(t2*t2/nu)
  b2 <- nu/(nu+t2*t2)
  sB2 <- sqrt(b2)
  ab2 <- ifelse(is.infinite(t2), 0, a2*b2)
  asB2 <- ifelse(is.infinite(t2), sign(t2), sign(t2)*sqrt(t2*t2/(nu+t2*t2)))
  if(nu==1){
    C1 <- pnorm(R) - (delta1>=0) + 2*OwenT(delta1*sB1, a1) -
      vapply(seq_len(J), function(i){
        2*OwenT(R[i], (a1*R[i]-delta1[i])/R[i]) +
          2*OwenT(delta1[i]*sB1, (delta1[i]*ab1-R[i])/b1/delta1[i])
      }, FUN.VALUE=numeric(1L))
    C2 <- pnorm(R) - (delta2>=0) + 2*OwenT(delta2*sB2, a2) -
      vapply(seq_len(J), function(i){
        2*OwenT(R[i], (a2*R[i]-delta2[i])/R[i]) +
          2*OwenT(delta2[i]*sB2, (delta2[i]*ab2-R[i])/b2/delta2[i])
      }, FUN.VALUE=numeric(1L))
    return(C2-C1)
  }
  nu <- as.integer(nu)
  n <- nu-1L
  H1 <- H2 <- M1 <- M2 <- matrix(numeric(n*J), nrow=n)
  H1[1L,] <- -dnorm(R) * pnorm(a1*R-delta1)
  M1[1L,] <- asB1*dnorm(delta1*sB1)*(pnorm(delta1*asB1)-pnorm((delta1*ab1-R)/sB1))
  H2[1L,] <- -dnorm(R) * pnorm(a2*R-delta2)
  M2[1L,] <- asB2*dnorm(delta2*sB2)*(pnorm(delta2*asB2)-pnorm((delta2*ab2-R)/sB2))
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
    C1 <- pnorm(R) - (delta1>=0) + 2*OwenT(delta1*sB1, a1) -
      vapply(seq_len(J), function(i){
        2*OwenT(R[i], (a1*R[i]-delta1[i])/R[i]) +
          2*OwenT(delta1[i]*sB1, (delta1[i]*ab1-R[i])/b1/delta1[i])
      }, FUN.VALUE=numeric(1L))
    C2 <- pnorm(R) - (delta2>=0) + 2*OwenT(delta2*sB2, a2) -
      vapply(seq_len(J), function(i){
        2*OwenT(R[i], (a2*R[i]-delta2[i])/R[i]) +
          2*OwenT(delta2[i]*sB2, (delta2[i]*ab2-R[i])/b2/delta2[i])
      }, FUN.VALUE=numeric(1L))
    indices <- seq(2L, n, by=2L)
    return(C2-C1 +
             2*.colSums(H2[indices,]+M2[indices,]-H1[indices,]-M1[indices,],
                        m=length(indices), n=J))
  }else{
    indices <- seq(1L, n, by=2L)
    return(pnorm(-delta2)-pnorm(-delta1) + sqrt(2*pi) *
             .colSums(H2[indices,]+M2[indices,]-H1[indices,]-M1[indices,],
                      m=length(indices), n=J))
  }
}
