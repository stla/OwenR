context("OwenQ1")

test_that("OwenQ1 for large R equals ptOwen", {
  expect_true(OwenQ1(4, 3, 2, 100) == ptOwen(3, 4, 2))
  expect_equal(OwenQ1(5, 3, 2, 100), ptOwen(3, 5, 2), tolerance=1e-15)
  expect_equal(OwenQ1(5, 3, 2, 100), pt(3, 5, 2), tolerance=1e-12)
})

test_that("OwenQ1 for t=+Inf does not depend on delta", {
  expect_true(OwenQ1(5, Inf, 2, 2) == OwenQ1(5, Inf, 3, 2))
  expect_true(OwenQ1(6, Inf, 2, 2) == OwenQ1(6, Inf, 3, 2))
  # nu=2 => moment truncated normal
  expect_equal(OwenQ1(2, Inf, 1, 2), sqrt(2*pi)*(dnorm(0)-dnorm(2)), tolerance=1e-15)
  # nu >=1 => incomplete Gamma
  R <- 2; nu <- 6
  expect_equal(OwenQ1(nu, Inf, 3, R), 1-gsl::gamma_inc(nu/2,R^2/2)/gamma(nu/2),
               tolerance=1e-14)
  # does not depend on t for delta=-Inf and the same result
  expect_true(OwenQ1(5, Inf, 2, 2) == OwenQ1(5, 1, -100, 2))
  expect_equal(OwenQ1(6, Inf, 2, 2), OwenQ1(6, 1, -100, 2), tolerance=1e-16)
})

test_that("OwenQ1 for t=-Inf equals 0", {
  expect_true(OwenQ1(5, -Inf, 2, 100) == 0)
  expect_true(OwenQ1(6, -Inf, 2, 100) == 0)
})

test_that("OwenQ1 - bivariate Student", {
  t1 <- 2; t2 <- 1; delta1 <- 3; delta2 <- 2
  nu <- 6
  R <- sqrt(nu)*(delta1 - delta2)/(t1-t2)
  owen <- OwenQ1(nu, t2, delta2, R) - OwenQ1(nu, t1, delta1, R)
  pmvt <- mvtnorm::pmvt(lower=c(t1,-Inf), upper=c(Inf,t2), delta=c(delta1, delta2),
                df=nu, corr=cbind(c(1,1),c(1,1)))
  expect_equal(owen, pmvt[1], tolerance=1e-4, check.attributes=FALSE)
  wolfram <- 0.01785518085912236
  expect_equal(owen, wolfram, tolerance=1e-11)
  #
  nu <- 5
  R <- sqrt(nu)*(delta1 - delta2)/(t1-t2)
  owen <- OwenQ1(nu, t2, delta2, R) - OwenQ1(nu, t1, delta1, R)
  wolfram <- 0.01868982415809893
  expect_equal(owen, wolfram, tolerance=1e-9)
  # wolfram input:
  # With[{t1 = 2, t2=1,delta1=3, delta2=2, nu=6},
  #      NProbability[
  #        (x+delta1)/Sqrt[y/nu] >= t1 && (x+delta2)/Sqrt[y/nu] <=t2, {x \[Distributed]
  #          NormalDistribution[], y\[Distributed]ChiSquareDistribution[nu]}]]
})

test_that("OwenQ2 - bivariate Student", {
  t1 <- 2; t2 <- 1; delta1 <- 3; delta2 <- 2
  nu <- 6
  R <- sqrt(nu)*(delta1 - delta2)/(t1-t2)
  owen <- - (ptOwen(t2, nu, delta2) - OwenQ1(nu, t2, delta2, R)) +
             (ptOwen(t1, nu, delta1) - OwenQ1(nu, t1, delta1, R))
  wolfram <- 0.03257737810540227
  expect_equal(owen, wolfram, tolerance=1e-10)
  #
  nu <- 5
  R <- sqrt(nu)*(delta1 - delta2)/(t1-t2)
  owen <- - (ptOwen(t2, nu, delta2) - OwenQ1(nu, t2, delta2, R)) +
    (ptOwen(t1, nu, delta1) - OwenQ1(nu, t1, delta1, R))
  wolfram <- 0.0353568969628651
  expect_equal(owen, wolfram, tolerance=1e-9)
})



