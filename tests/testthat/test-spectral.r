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
  expect_equal(res1, 1)
  expect_equal(res2, 2)
})

test_that("trim_epoch works", {
  out <- trim_epoch(dat1, c(30,70))
  expect_false(FALSE %in% c(out[,1] > 30 & out[,1] < 70))
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

test_that("repartition works", {
  out <- repartition(Y3, list(c(1,2), 3))
  expect_equal(nrow(out[[1]]), 100)
  expect_equal(nrow(out[[2]]), 50)
})

test_that("LSSA_LFI_model fits results below 1e-12 threshold", {
  k <- 10
  out <- LSSA_LFI(dat1, n_iter = k, intercept = TRUE)
  freqs <- out$freqs
  coefs <- out$coefs

  y <- coefs[1,3]
  rss_R <- numeric(k)
  for (i in 1:k) {
    theta <- dat1[,1] * 2 * pi * freqs[i]
    y <- y + coefs[i,1] * cos(theta) + coefs[i,2] * sin(theta) 
    rss_R[i] <- sum((y - dat1[,2])^2)
  }

  mod <- LSSA_LFI_model(dat1, t_ = dat1[,1], n_iter = k, intercept = TRUE)

  expect_true(all(abs(y - mod$y) < 1e-12))
})

test_that("AIC values correct in LSSA_LFI_multi - zero intercept", {
  k <- 1:10
  aic_out <- numeric(length(k))
  aic_out_2 <- numeric(length(k))

  for (i in 1:length(k)) {
    out <- LSSA_LFI_multi(Y3, n_iter = k[i], intercept = FALSE, AIC = FALSE)
    out_2 <- LSSA_LFI_multi(Y1, n_iter = k[i], intercept = FALSE, AIC = FALSE)

    mod1 <- LSSA_LFI_model(dat1, t_ = dat1[,1], intercept = FALSE, n_iter = k[i])
    mod2 <- LSSA_LFI_model(dat2, t_ = dat2[,1], intercept = FALSE, n_iter = k[i])
    mod3 <- LSSA_LFI_model(dat3, t_ = dat3[,1], intercept = FALSE, n_iter = k[i])
    rss1 <- sum((dat1[,2] - mod1[,2])^2)
    rss2 <- sum((dat2[,2] - mod2[,2])^2)
    rss3 <- sum((dat3[,2] - mod3[,2])^2)
    n <- nrow(dat1) + nrow(dat2) + nrow(dat3)
    n_2 <- nrow(dat1) + nrow(dat2)

    aic <-   2 * 3 * (k[i] * 3 + 1) + n * log((rss1 + rss2 + rss3)/n) + n + n * log(2 * pi)
    aic_2 <- 2 * 2 * (k[i] * 3 + 1) + n_2 * log((rss1 + rss2)/n_2) + n_2 + n_2 * log(2 * pi)

    aic_out[i] <- out$aic[k[i]] - aic
    aic_out_2[i] <- out_2$aic[k[i]] - aic_2
  }

  expect_true(all(abs(aic_out) < 1e-12))
  expect_true(all(abs(aic_out_2) < 1e-12))
})

test_that("AIC values correct in LSSA_LFI_multi - including intercept", {
  k <- 1:10
  aic_out <- numeric(length(k))
  aic_out_2 <- numeric(length(k))

  for (i in 1:length(k)) {
    out <- LSSA_LFI_multi(Y3, n_iter = k[i], intercept = TRUE, AIC = FALSE)
    out_2 <- LSSA_LFI_multi(Y1, n_iter = k[i], intercept = TRUE, AIC = FALSE)

    mod1 <- LSSA_LFI_model(dat1, t_ = dat1[,1], intercept = TRUE, n_iter = k[i])
    mod2 <- LSSA_LFI_model(dat2, t_ = dat2[,1], intercept = TRUE, n_iter = k[i])
    mod3 <- LSSA_LFI_model(dat3, t_ = dat3[,1], intercept = TRUE, n_iter = k[i])
    rss1 <- sum((dat1[,2] - mod1[,2])^2)
    rss2 <- sum((dat2[,2] - mod2[,2])^2)
    rss3 <- sum((dat3[,2] - mod3[,2])^2)

    n <- nrow(dat1) + nrow(dat2) + nrow(dat3)
    n_2 <- nrow(dat1) + nrow(dat2)

    aic <- 2 * 3 * (k[i] * 3 + 2) + n * log((rss1 + rss2 + rss3)/n) + n + n * log(2 * pi)
    aic_2 <- 2 * 2 * (k[i] * 3 + 2) + n_2 * log((rss1 + rss2)/n_2) + n_2 + n_2 * log(2 * pi)

    aic_out[i] <- out$aic[k[i]] - aic
    aic_out_2[i] <- out_2$aic[k[i]] - aic_2
  }

  expect_true(all(abs(aic_out) < 1e-12))
  expect_true(all(abs(aic_out_2) < 1e-12))
})

test_that("LSSA_LFI_comp chooses better fitting model than LSSA_LFI_multi", {

  Y3_repar <- repartition(Y3, list(c(1,2), 3))
  Y3_repar_multi <- LSSA_LFI_multi(Y3_repar, n_iter = 100, intercept = FALSE, AIC = TRUE)
  aic1 <- Y3_repar_multi["aic"]

  Y3_repar_multi <- LSSA_LFI_multi(Y3_repar, n_iter = 100, intercept = FALSE, AIC = TRUE)
  out <- LSSA_LFI_comp(Y3_repar, intercept = FALSE) 
  aic2 <- out$aic

  expect_true(aic2 < aic1)
})

