context("OwenQ1")

test_that("OwenQ1 for large R equals ptOwen", {
  expect_true(OwenQ1(4, 3, 2, 100) == ptOwen(3, 4, 2))
  expect_equal(OwenQ1(5, 3, 2, 100), ptOwen(3, 5, 2), tolerance=1e-15)
  expect_equal(OwenQ1(5, 3, 2, 100), pt(3, 5, 2), tolerance=1e-12)
})
