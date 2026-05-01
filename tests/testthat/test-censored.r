set.seed(42)

x1 <- c(0,1,2,2,5,7)
x2 <- c(0,0,4,15,7,99)

x <- matrix(c(x1, x2), ncol = 2)

test_that("Estimation of Poisson correct in lambda_col", {
  res1 <- pois_rcens(x, omit_zero = FALSE)
  res2 <- pois_rcens(x, omit_zero = TRUE)
  expect_equal(res1[1], 6.06)
  expect_equal(res1[2], 51.70)
  expect_equal(res2[1], 6.37)
  expect_equal(res2[2], 51.76)
})






