#' Empirical Wasserstein R1 Distance betwen Two e.c.d.f.
#'
#' For two or more one-dimensional point processes, computes their Wasserstein distance via their empirical cumuluative density functions (e.c.d.f.), adjusting for sample size and the scale of their durations \insertCite{vasershtein_markov_1969,anderson_distribution_1962,ramdas_wasserstein_2017}{arkhaia}, using numerical quadrature: \eqn{W_{1s} = \sqrt{\frac{n_1 n_2}{n_1 + n_2}} \int_\mathbf{R} \vert F_1(x) - F_2(x) \vert \, dx}. More than two processes result in pairwise estimation of processes.
#' 
#' @param x1 Either a list of vectors or a vector containing observations of sample 1.
#' @param x2 If \code{x1} is vector, \code{x2} is vector containing observations of sample 2.
#' @param scaled Whether to scale the resulting Wasserstein R1 distance to the distributions of \code{x1}, \code{x2}, ensuring a maximum of 1. Default is \code{TRUE}.
#' @param by The increment of sampling for performing quadrature. Default is \code{1e-3}.
#' @param perm Whether to perform a Monte Carlo permutation across all samples in the input. Default is \code{FALSE}.
#' @param n_perm Number of Monte Carlo permutations to perform. Default is \code{10^5}.
#' @returns For two vectors \code{x1} and \code{x2}, the empirical scaled Wassertein R1 distance between the two e.c.d.f. of their samples. For a list of vectors input in \code{x1}, a matrix of the empirical scaled Wasserstein R1 distances pairwise between all vectors.
#' 
#' @references
#'   \insertAllCited{}
#' 
#' @importFrom Rdpack reprompt
#' @export
wass_R1 <- function(x1, x2 = NULL, scaled = TRUE, by = 1e-3, perm = FALSE, n_perm = 10^5) {
    UseMethod("wass_R1")
}

#' @rdname wass_R1
#' @export
wass_R1.numeric <- function(x1, x2 = NULL, scaled = TRUE, by = 1e-3, perm = FALSE, n_perm = 10^5) {
  if (!is.vector(x2)) {
    stop("x2 must be provided as a vector.")
  }
  if (scaled == TRUE) {
    sc <- 1
  } else {
    sc <- 0
  }
  x <- c(x1,x2)
  n <- c(length(x1),length(x2))
  W <- wass_R1_arma(x, n, sc, by)
  out <- W[2,1]
  if (perm == TRUE) {
    W_perm <- wass_R1_perm_arma(x, n, sc, by, n_perm)
    out <- list(W = W[2,1], W_perm = W_perm[2,1,])
    class(out) <- c("wass_perm", "list")
  } else {
    class(out) <- c("wass", "numeric")
  }

  return(out)
}

#' @rdname wass_R1
#' @export
wass_R1.list <- function(x1, x2 = NULL, scaled = TRUE, by = 1e-3, perm = FALSE, n_perm = 10^5) {
  if (scaled == TRUE) {
    sc <- 1
  } else {
    sc <- 0
  }
  x <- unlist(x1)
  n <- sapply(x1, length)
  W <- wass_R1_arma(x, n, sc, by)

  if (length(x1) == 2) {
    out <- W[2,1]
  } else {
    W[upper.tri(W, diag = TRUE)] <- NA
    out <- W
  }

  if (perm == TRUE) {
    W_perm <- wass_R1_perm_arma(x, n, sc, by, n_perm)
    out <- list(W = W , W_perm = W_perm)
    class(out) <- c("wass_perm", "list")
  } else {
    if (length(x1) == 2) {
      class(out) <- c("wass", "numeric")
    } else {
      class(out) <- c("wass", "matrix", "array")
    }
  }

  return(out)
}

#' @export 
print.wass <- function(x, ...) {
  if (is.vector(x)) {
    cat("Empirical Scaled Wasserstein Distance: ", "\n")
    cat(format(x, digits = 2),"\n")
  } else if (is.matrix(x)) {
    cat("Empirical Scaled Wasserstein Distances: ", "\n")
    print(format(x, digits = 2), quote = FALSE)
    cat("\n")
  }
}


#' @export 
print.wass_perm <- function(x, ...) {
  if (is.vector(x$W)) {
    cat("Empirical Scaled Wasserstein Distance: ", "\n")
    cat(format(x$W, digits = 2),"\n")
    cat("Mean of Permutation Samples: ", "\n")
    cat(format(mean(x$W_perm), digits = 2),"\n")
  } else if (is.matrix(x$W)) {
    cat("Empirical Scaled Wasserstein Distances:", "\n")
    print(format(x$W, digits = 2), quote = FALSE)
    cat("\n")
    cat("Means of Permutation Samples: ", "\n")
    tmp <- apply(x$W_perm, c(1,2), mean)
    tmp[upper.tri(tmp, diag = TRUE)] <- NA
    print(format(tmp, digits = 2), quote = FALSE)
    cat("\n")
  }
}





#' Integration in Distribution over Time
#'
#' For two or more one-dimensional point processes, this function estimates that the probability of processes are integrated in distribution over time (i.d.t.) by computing the probability of Monte Carlo-permuted empirical scaled Wasserstein R1 distances greater than the distance obtained (see \code{\link[arkhaia]{wass_R1}}).
#' 
#' @param x1 Either a list of vectors or a vector containing observations of sample 1.
#' @param x2 If \code{x1} is vector, \code{x2} is vector containing observations of sample 2.
#' @param by The increment of sampling for performing quadrature. Default is \code{1e-3}.
#' @param n_perm Number of Monte Carlo permutations to perform. Default is \code{10^5}.
#' @returns For two vectors \code{x1} and \code{x2}, the probability that Monte Carlo permutation of samples results in a lower  empirical scaled Wassertein R1 distance than the one obtained. For a list of vectors input in \code{x1}, a matrix of probabilities pairwise between all vectors.
#' 
#' @importFrom Rdpack reprompt
#' @export
idt <- function(x1, x2 = NULL, by = 1e-3, n_perm = 10^5) {
    UseMethod("idt")
}

#' @rdname idt
#' @export
idt.numeric <- function(x1, x2 = NULL, by = 1e-3, n_perm = 10^5) {
  res <- wass_R1(x1, x2, scaled = TRUE, by = by, perm = TRUE, n_perm = n_perm)
  w <- res$W
  w_perm <- res$W_perm
  if (is.vector(w)) {
    pr <- mean(w_perm > w)
    class(pr) <- c("pr_idt", "numeric")
  } else {
    w0 <- replicate(n_perm, w, simplify="array")
    w1 <- w_perm > w0
    pr <- apply(w1, c(1,2), mean)
    class(pr) <- c("pr_idt", "matrix", "array")
  }
  return(pr)
}

#' @rdname wass_R1
#' @export
idt.list <- function(x1, x2 = NULL, by = 1e-3, n_perm = 10^5) {
  idt.numeric(x1, x2 = NULL, by = 1e-3, n_perm = 10^5)
}

#' @export 
print.pr_idt <- function(x, ...) {
  if (is.vector(x)) {
    cat("Probability of i.d.t.: ", "\n")
    cat(format(x, digits = 2),"\n")
  } else if (is.matrix(x)) {
    cat("Probabilities of i.d.t.: ", "\n")
    print(format(x, digits = 2), quote = FALSE)
    cat("\n") 
  }
}






#' Thresholds of Integration in Distribution over Time
#'
#' For two one-dimensional point processes, this function estimates the condition of the two processes becoming i.d.t. at a moment of time \eqn{t} considering a window of time around that moment, \eqn{t-h,t+h} (see \code{\link[arkhaia]{idt}}). This is achieved by computing the difference of the empirical scaled Wasserstein R1 distances between the two processes before and after \eqn{t}. Monte Carlo permutation of samples is used to estimate the probability that the obtain distances is of practical signficance.
#' 
#' @param x1 Either a list of vectors or a vector containing observations of sample 1.
#' @param x2 If \code{x1} is vector, \code{x2} is vector containing observations of sample 2.
#' @param by The increment of sampling for performing quadrature. Default is \code{1e-3}.
#' @param incr The increment at which to sample values of \eqn{t} and \eqn{h}. Default is 0.1. 
#' @param perm Whether to perform a Monte Carlo permutation across all samples in the input, which will return the probability of an integrating threshold at time \eqn{t}. Default is \code{FALSE}.
#' @param n_perm Number of Monte Carlo permutations to perform. Default is \code{10^5}.
#' @returns For two vectors \code{x1} and \code{x2}, the probability that Monte Carlo permutation of samples results in a lower  empirical scaled Wassertein R1 distance than the one obtained. For a list of vectors input in \code{x1}, a matrix of probabilities pairwise between all vectors.
#' 
#' @importFrom Rdpack reprompt
#' @export
th_idt <- function(x1, x2 = NULL, by = 1e-3, incr = 0.1, perm = FALSE, n_perm = 10^5) {
    UseMethod("th_idt")
}

#' @rdname th_idt
#' @export
th_idt.numeric <- function(x1, x2 = NULL, by = 1e-3, incr = 0.1, perm = FALSE, n_perm = 10^5) {
  if (is.vector(x2)) {
    a0 <- min(c(x1,x2))
    a1 <- max(c(x1,x2))
  } else if (is.list(x1)) {
    x1 <- x1[[1]]
    x2 <- x2[[2]]
    a0 <- min(unlist(x1))
    a1 <- max(unlist(x1))
  } else {
    stop("x1 must be list, or x1 and x2 must be vectors.")
  }
  h_max <- (a1-a0)/2
  t_ <- seq(a0, a1, by = incr)
  h_ <- seq(incr, h_max, by = incr)
  if (perm == FALSE) {
    res <- matrix(NA, length(h_), length(t_))
    for (i in 1:length(t_)) {
      threshold <- t_[i]
      for (j in 1:length(h_)) {
        h <- h_[j]

        if ((threshold - h) >= min(t_) & (threshold + h) <= max(t_)) {

          x_11 <- x1[x1 < threshold & x1 > threshold - h]
          x_12 <- x1[x1 >= threshold & x1 < threshold + h]

          x_21 <- x2[x2 < threshold & x2 > threshold - h]
          x_22 <- x2[x2 >= threshold & x2 < threshold + h]

          if ((length(x_11) >= 1 & length(x_12) >= 1) & (length(x_21) >= 1 & length(x_22) >= 1)) {

          bef <- wass_R1(x_11, x_21, scaled = TRUE, by = by, perm = FALSE)
          aft <- wass_R1(x_12, x_22, scaled = TRUE, by = by, perm = FALSE)
          res[j,i] <- bef - aft
          }
        }
      }
    }
    hh <- rep(h_, length(t_))
    tt <- rep(t_, each=length(h_))
    out <- data.frame(t = tt, h = hh, theta = c(res))
    class(out) <- c("th_idt", "data.frame")
    return(out)

  } else {
    res <- matrix(NA, length(h_), length(t_))
    perm <- matrix(NA, length(h_), length(t_))

    for (i in 1:length(t_)) {
      threshold <- t_[i]
      for (j in 1:length(h_)) {
        h <- h_[j]

        if ((threshold - h) >= min(t_) & (threshold + h) <= max(t_)) {

          x_11 <- x1[x1 < threshold & x1 > threshold - h]
          x_12 <- x1[x1 > threshold & x1 < threshold + h]

          x_21 <- x2[x2 < threshold & x2 > threshold - h]
          x_22 <- x2[x2 > threshold & x2 < threshold + h]

          if ((length(x_11) >= 1 & length(x_12) >= 1) & (length(x_21) >= 1 & length(x_22) >= 1)) {

          bef <- wass_R1(x_11, x_21, scaled = TRUE, by = by, perm = TRUE, n_perm = n_perm)
          aft <- wass_R1(x_12, x_22, scaled = TRUE, by = by, perm = TRUE, n_perm = n_perm)
          d_ <- bef$W - aft$W
          res[j,i] <- d_
          d_perm <- bef$W_perm - aft$W_perm
          perm[j,i] <- mean(d_perm < d_)
          
          }
        }
      }
    }

    hh <- rep(h_, length(t_))
    tt <- rep(t_, each=length(h_))
    out1 <- data.frame(t = tt, h = hh, theta = c(res))
    out2 <- data.frame(t = tt, h = hh, theta = c(perm))
    class(out1) <- c("th_idt", "data.frame")
    class(out2) <- c("th_idt", "data.frame")

    out <- list(theta = out1, pr_theta = out2)

    class(out) <- c("th_idt", "list")
    return(out)
  }
}


#' @rdname th_idt
#' @export
th_idt.list <- function(x1, x2 = NULL, by = 1e-3, incr = 0.1, perm = FALSE, n_perm = 10^5) {
  th_idt.numeric(x1, x2, by, incr, perm, n_perm)
}


















