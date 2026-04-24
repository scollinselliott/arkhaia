set.seed(42)

t_ <- sort(runif(50,0,100))
x1 <- 2 + 5 * cos(2 * pi * 0.03 * t_) - 0.25 * sin(2 * pi * 0.03 * t_) + rnorm(50)
x2 <- x1 + rnorm(50)
x3 <- 1 + 3 * cos(2 * pi * 0.001 * t_) + 7 * sin(2 * pi * 0.001 * t_) + rnorm(50)

dat1 <- data.frame(t = t_, y = x1)
dat2 <- data.frame(t = t_, y = x2)
dat3 <- data.frame(t = t_, y = x3)

Y1 <- list(dat1, dat2)
Y2 <- list(dat1, dat3)
Y3 <- list(dat1, dat2, dat3)


test_that("LSSA functions: S3 class check", {
  expect_s3_class(LSSA(dat1), "data.frame")
  expect_s3_class(LSSA_LFI(dat1,), "LSSA_LFI")
  expect_s3_class(LSSA_LFI_candidates(Y1), "LSSA_LFI_AIC")
  expect_s3_class(LSSA_LFI_model(dat1), "data.frame")
})

test_that("LSSA detects frequency", {
  res <- LSSA(dat1)
  freq_ <- res$freq[which.max(res$power)]
  expect_equal(freq_, 0.03)
})

test_that("Correct candidate model selected", {
  res1 <- LSSA_LFI_candidates(Y1, n_iter = 50)
  res2 <- LSSA_LFI_candidates(Y2, n_iter = 50)
  mod1 <- model_select(res1)
  mod2 <- model_select(res2)
  expect_equal(mod1, 1)
  expect_equal(mod2, 2)
})

test_that("trim_epoch works", {
  out <- trim_epoch(dat1, c(30,70))
  expect_true(out[,1] > 30 & out[,1] < 70)
})

test_that("Correct validated candidate model selected", {
  res <- LSSA_LFI_validated(Y3, pair = c(1,2), n_iter = 50) 
  expect_equal(res, 1)
})

test_that("Pairwise validated candidate model selection works", {
  expected <- matrix(NA, 3, 3)
  expected[1,2] <- 1
  expected[1:2,3] <- 0
  res <- LSSA_LFI_pairwise(Y3, n_iter = 50) 
  expect_equal(res, expected)
})






