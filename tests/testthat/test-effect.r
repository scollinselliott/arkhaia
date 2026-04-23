
# sequences

x1 <- c(2, 3, 10, 11, 5)
x2 <- c(1, 1, 17, 23, 3)
x3 <- c(0, 0, 2, 81, 11)
 
X_dep <- matrix(c(x1, x1), ncol = 2)
X_indep <- matrix(c(x1, x2, x3), ncol = 3)


# Cressie-Read

test_that("CR works", {
  expect_true(CR(X_dep) == 0)
  expect_true(CR(X_indep) > 0)
})

# Bias-Corrected Cramer's V

test_that("VB works", {
  expect_true(VB(X_dep) == 0)
  expect_true(VB(X_indep) > 0)
})




