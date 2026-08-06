#' Wasserstein R1 Distance betwen Two e.c.d.f.
#'
#' For two or more one-dimensional point processes, computes their Wasserstein distance via their empirical cumuluative density functions (e.c.d.f.), adjusting for sample size and the scale of their durations \insertCite{vasershtein_markov_1969,anderson_distribution_1962,ramdas_wasserstein_2017}{arkhaia}, using numerical quadrature:
#' \deqn{\displaystyle{W_{1s} = \frac{1}{s} \int_\mathbf{R} \vert F_1(x) - F_2(x) \vert \, dx}}
#' where \eqn{s} is a scaling coefficient, namely either the range of the point process or its support. More than two processes result in pairwise estimation of processes.
#' 
#' @param x1 Either a list of vectors or a vector containing observations of sample 1.
#' @param x2 If \code{x1} is vector, \code{x2} is vector containing observations of sample 2.
#' @param scaled Whether to scale the resulting Wasserstein R1 distance. Default is \code{FALSE}.
#' @param scaled_by The magnitude by which to scale the Wasserstein R1 distance, if \code{scaled} is \code{TRUE}. Entering \code{0} (the default) will use the range of the maximum and minimum of the pair of point processes.
#' @param quad_by The increment of sampling for performing quadrature. Default is \code{1e-3}.
#' @param perm Whether to perform a Monte Carlo permutation across all samples in the input. Default is \code{FALSE}.
#' @param n_perm Number of Monte Carlo permutations to perform. Default is \code{10^5}.
#' @returns For two vectors \code{x1} and \code{x2}, the empirical scaled Wassertein R1 distance between the two e.c.d.f. of their samples. For a list of vectors input in \code{x1}, a matrix of the empirical scaled Wasserstein R1 distances pairwise between all vectors.
#' 
#' @references
#'   \insertAllCited{}
#' 
#' @importFrom Rdpack reprompt
#' @export
wass_R1 <- function(x1, x2 = NULL, scaled = FALSE, scaled_by = 0, quad_by = 1e-3, perm = FALSE, n_perm = 10^5) {
    UseMethod("wass_R1")
}

#' @rdname wass_R1
#' @export
wass_R1.numeric <- function(x1, x2 = NULL, scaled = FALSE, scaled_by = 0, quad_by = 1e-3, perm = FALSE, n_perm = 10^5) {
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
  
  W <- wass_R1_arma(x, n, sc, scaled_by, quad_by)
  if (perm == FALSE) {
    out <- list(W = W[2,1], scaled = scaled, perm = perm)
  } else {
    W_perm <- wass_R1_perm_arma(x, n, sc, scaled_by, quad_by, n_perm)
    out <- list(W = W[2,1], W_perm = W_perm[2,1,], scaled = scaled, perm = perm)
  } 
  class(out) <- c("wass", "list")
  return(out)
}

#' @rdname wass_R1
#' @export
wass_R1.list <- function(x1, x2 = NULL, scaled = FALSE, scaled_by = 0, quad_by = 1e-3, perm = FALSE, n_perm = 10^5) {
  if (!is.null(x2)) {
    message("x1 is list; ignoring x2 input.")
  }
  if (scaled == TRUE) {
    sc <- 1
  } else {
    sc <- 0
  }
  x <- unlist(x1)
  n <- sapply(x1, length)
  W <- wass_R1_arma(x, n, sc, scaled_by, by)
  if (length(x1) == 2) {
    W <- W[2,1]
  } else {
    W[upper.tri(W, diag = TRUE)] <- NA
  }
  if (perm == TRUE) {
    W_perm <- wass_R1_perm_arma(x, n, sc, scaled_by, by, n_perm)
    if (length(x1) == 2) {
      W_perm <- W_perm[2,1,]
    } else {
      tmp <- matrix(0, length(x1), length(x1))
      tmp[upper.tri(tmp, diag = TRUE)] <- NA
      tmp0 <- replicate(n_perm, tmp, simplify="array")
      W_perm <- W_perm + tmp0
    }
    out <- list(W = W, W_perm = W_perm, scaled = scaled, perm = perm)
  } else {
    out <- list(W = W, scaled = scaled, perm = perm)
  }
  class(out) <- c("wass", "list")
  return(out)
}

#' @export 
print.wass <- function(x, ...) {
  if (x$scaled == TRUE & !is.matrix(x$W)) {
    txt <- "Empirical Scaled Wasserstein Distance: "
  } else if (x$scaled == TRUE & is.matrix(x$W)) {
    txt <- "Empirical Scaled Wasserstein Distances: "
  } else if (x$scaled == FALSE & !is.matrix(x$W)) {
    txt <- "Empirical Wasserstein Distance: "
  } else if (x$scaled == FALSE & is.matrix(x$W)) {
    txt <- "Empirical Wasserstein Distances: "
  }
  cat(txt, "\n")
  if (!is.matrix(x$W)) {
    cat(format(x$W, digits = 2),"\n")
  } else if (is.matrix(x$W)) {
    print(format(x$W, digits = 2), quote = FALSE)
    cat("\n")
  }
  if (x$perm == TRUE) {
    if (!is.matrix(x$W)) {
    cat("Mean of Permutation Samples: ", "\n")
    cat(format(mean(x$W_perm), digits = 2))
    cat("\n")
    } else if (is.matrix(x$W)) {
    cat("Means of Permutation Samples: ", "\n")
    tmp <- apply(x$W_perm, c(1,2), mean)
    tmp[upper.tri(tmp, diag = TRUE)] <- NA
    print(format(tmp, digits = 2), quote = FALSE)
    cat("\n")
    }
  }
}



#' Kolmogorov-Smirnov Distance
#'
#' This function computes the Kolmogorov-Smirnov (KS) distance between two empirical cumulative distribution functions. For more functionality, use \code{\link[stats]{ks.test}}.
#' 
#' The Rcpp-based function \code{\link[arkhaia]{ks_arkhaia}} is used for faster estimation of the KS distance than \code{\link[stats]{ks.test}}.
#' 
#' @param x1 A vector of observations, or a list containing both vectors.
#' @param x2 A vector of observations, if \code{x1} is not a list. 
#' @returns The KS distance between the two vectors.
#' 
#' @importFrom Rdpack reprompt
#' @export
ks_arkhaia <- function(x1, x2 = NULL) {
    UseMethod("ks_arkhaia")
}

#' @rdname ks_arkhaia
#' @export
ks_arkhaia.numeric <- function(x1, x2 = NULL) {
  if (is.null(x2)) {
    stop("x1 and x2 must be vectors")
  }
  res <- ks_Rcpp(x1,x2)
  return(res)
}

#' @rdname ks_arkhaia
#' @export
ks_arkhaia.list <- function(x1, x2 = NULL) {
  y1 <- x1[[1]]
  y2 <- x1[[2]]
  res <- ks_Rcpp(y1,y2)
  return(res)
}




#' Thresholds of Integration in Distribution over Time (KS)
#'
#' For two one-dimensional temporal point processes \eqn{X_1} and \eqn{X_2}, this function estimates the condition of the two processes becoming "integrated in distribution over time" (i.d.t.) at a moment in time \eqn{t}. This is achieved by computing the Kolmogorov-Smirnov distance \eqn{D} between their empirical culumative distribution functions (e.c.d.f.), before and after \eqn{t} for some duration \eqn{h}, and taking the difference between the distances \eqn{\theta = D_{before}} - \eqn{D_{after}}. For practical significance in the estimation of \eqn{\theta}, Monte Carlo permutation of samples is used to estimate the probability that the observed value of \eqn{\theta} is greater than the distribution of permuted samples (i.e., that observations are interchangable between \eqn{X_1} and \eqn{X_2}). 
#' 
#' @param x1 Either a list of vectors or a vector containing observations of sample 1.
#' @param x2 If \code{x1} is vector, \code{x2} is vector containing observations of sample 2.
#' @param stat Whether to multiply the distance by \eqn{\sqrt{\frac{n_1 n_2}{n_1 + n_2}}}. Default is \code{FALSE}.
#' @param incr The increment at which to sample values of \eqn{t} and \eqn{h}. Default is 0.5. 
#' @param perm Whether to perform a Monte Carlo permutation across all samples in the input, which will return the probability of an integrating threshold at time \eqn{t}. Default is \code{FALSE}.
#' @param n_perm Number of Monte Carlo permutations to perform. Default is \code{1e5}.
#' @returns A list containing a data frame  two vectors \code{x1} and \code{x2}, the probability that Monte Carlo permutation of samples results in a lower distance than the one obtained. For a list of vectors input in \code{x1}, a matrix of probabilities pairwise between all vectors.
#' 
#' @importFrom Rdpack reprompt
#' @export
th_idt_ks <- function(x1, x2 = NULL, stat = FALSE, incr = 0.5, perm = FALSE, n_perm = 1e5) {
    UseMethod("th_idt_ks")
}

#' @rdname th_idt_ks
#' @export
th_idt_ks.numeric <- function(x1, x2 = NULL, stat = FALSE, incr = 0.5, perm = FALSE, n_perm = 1e5) {
  if (is.null(x2)) {
    stop("x1 and x2 must be vectors")
  }
  if (stat == FALSE) {
    res <- th_idt_ks_arma(x1,x2,incr)
  } else {
    res <- th_idt_ks_stat_arma(x1,x2,incr)
  }
  mat <- res[[1]]
  mat[mat > 8000] <- NA
  t_ <- c(res[[2]])
  h_ <- c(res[[3]])
  hh <- rep(h_, length(t_))
  tt <- rep(t_, each=length(h_))
  out1 <- data.frame(t = tt, h = hh, theta = c(mat))
  if (perm == TRUE) {
    if (stat == FALSE) {
      res_perm <- th_idt_ks_perm_arma(x1,x2,incr,n_perm)
    } else {
      res_perm <- th_idt_ks_stat_perm_arma(x1,x2,incr,n_perm)
    }
    res_perm[res_perm > 8] <- NA
    res0 <- replicate(n_perm, mat, simplify="array")
    res_perm0 <- (res0^2 < res_perm^2)
    res2 <- apply(res_perm0, c(1,2), mean, na.rm = TRUE)
    out2 <- data.frame(t = tt, h = hh, pr_theta = c(res2))
    out <- list(theta = out1, pr_theta = out2)
  } else {
    out <- list(theta = out1)
  }
  class(out) <- c("th_idt", "list")
  return(out)
}

#' @rdname th_idt_ks
#' @export
th_idt_ks.list <- function(x1, x2 = NULL, stat = FALSE, incr = 0.5, perm = FALSE, n_perm = 10e5) {
  y1 <- x1[[1]]
  y2 <- x1[[2]]
  if (stat == FALSE) {
    res <- th_idt_ks_arma(y1,y2,incr)
  } else {
    res <- th_idt_ks_stat_arma(y1,y2,incr)
  }
  mat <- res[[1]]
  mat[mat > 8000] <- NA
  t_ <- c(res[[2]])
  h_ <- c(res[[3]])
  hh <- rep(h_, length(t_))
  tt <- rep(t_, each=length(h_))
  out1 <- data.frame(t = tt, h = hh, theta = c(mat))
  if (perm == TRUE) {
    if (stat == FALSE) {
      res_perm <- th_idt_ks_perm_arma(y1,y2,incr,n_perm)
    } else {
      res_perm <- th_idt_ks_stat_perm_arma(y1,y2,incr,n_perm)
    }
    res_perm[res_perm > 8] <- NA
    res0 <- replicate(n_perm, mat, simplify="array")
    res_perm0 <- (res0^2 < res_perm^2)
    res2 <- apply(res_perm0, c(1,2), mean, na.rm = TRUE)
    out2 <- data.frame(t = tt, h = hh, pr_theta = c(res2))
    out <- list(theta = out1, pr_theta = out2)
  } else {
    out <- list(theta = out1)
  }
  class(out) <- c("th_idt", "list")
  return(out)
}






#' Thresholds of Integration in Distribution over Time (Wasserstein)
#'
#' For two one-dimensional point processes, this function estimates the condition of the two processes becoming i.d.t. at a moment of time \eqn{t} considering a window of time around that moment, \eqn{(t-h,t+h)}. This is achieved by computing the difference of the empirical scaled Wasserstein R1 distances between the two processes before and after \eqn{t}. Monte Carlo permutation of samples is used to estimate the probability that the obtain distances is of practical signficance.
#' 
#' @param x1 Either a list of vectors or a vector containing observations of sample 1.
#' @param x2 If \code{x1} is vector, \code{x2} is vector containing observations of sample 2.
#' @param stat Whether to multiply the distance by \eqn{\sqrt{\frac{n_1 n_2}{n_1 + n_2}}}. Default is \code{FALSE}.
#' @param scaled Whether to scale the Wasserstein R1 distance used in the course of estimation. Default is \code{FALSE}.
#' @param scaled_by The magnitude by which to scale the Wasserstein R1 distance, if \code{scaled} is \code{TRUE}. Entering \code{0} (the default) will use the range of the maximum and minimum of the pair of point processes, while entering \code{-1} will use the value of \eqn{h}.
#' @param quad_by The increment of sampling for performing quadrature. Default is \code{1e-3}.
#' @param incr The increment at which to sample values of \eqn{t} and \eqn{h}. Default is 0.5. 
#' @param perm Whether to perform a Monte Carlo permutation across all samples in the input, which will return the probability of an integrating threshold at time \eqn{t}. Default is \code{FALSE}.
#' @param n_perm Number of Monte Carlo permutations to perform. Default is \code{10^5}.
#' @returns For two vectors \code{x1} and \code{x2}, the probability that Monte Carlo permutation of samples results in a lower  empirical scaled Wassertein R1 distance than the one obtained. For a list of vectors input in \code{x1}, a matrix of probabilities pairwise between all vectors.
#' 
#' @importFrom Rdpack reprompt
#' @export
th_idt_wass <- function(x1, x2 = NULL, stat = FALSE, scaled = FALSE, scaled_by = 0, quad_by = 1e-3, incr = 0.5, perm = FALSE, n_perm = 10^5) {
    UseMethod("th_idt_wass")
}

#' @rdname th_idt_wass
#' @export
th_idt_wass.numeric <- function(x1, x2 = NULL, stat = FALSE, scaled = FALSE, scaled_by = 0, quad_by = 1e-3, incr = 0.5, perm = FALSE, n_perm = 10^5) {
  if (is.null(x2)) {
    stop("x1 and x2 must be vectors")
  }
  if (!(scaled_by == -1 | scaled_by >= 0 )) {
    stop("scaled_by must be -1, 0, or greater than 0.")
  }
  if (stat == FALSE) {
    res <- th_idt_wass_arma(x1,x2,scaled,scaled_by,quad_by,incr)
  } else {
    res <- th_idt_wass_stat_arma(x1,x2,scaled,scaled_by,quad_by,incr)
  }
  mat <- res[[1]]
  mat[mat > 8000] <- NA
  t_ <- c(res[[2]])
  h_ <- c(res[[3]])
  hh <- rep(h_, length(t_))
  tt <- rep(t_, each=length(h_))
  out1 <- data.frame(t = tt, h = hh, theta = c(mat))

  if (perm == TRUE) {
    if (stat == FALSE) {
      res_perm <- th_idt_wass_perm_arma(x1,x2,scaled,scaled_by,quad_by,incr,n_perm)
    } else {
      res_perm <- th_idt_wass_stat_perm_arma(x1,x2,scaled,scaled_by,quad_by,incr,n_perm)
    }

    res_perm[res_perm > 8000] <- NA
    res0 <- replicate(n_perm, mat, simplify="array")
    res_perm0 <- (res0^2 < res_perm^2)
    res2 <- apply(res_perm0, c(1,2), mean, na.rm = TRUE)
    out2 <- data.frame(t = tt, h = hh, pr_theta = c(res2))
    out <- list(theta = out1, pr_theta = out2)
  } else {
    out <- list(theta = out1)
  }
  class(out) <- c("th_idt", "list")
  return(out)

  #####################


  # if (is.vector(x2)) {
  #   a0 <- min(c(x1,x2))
  #   a1 <- max(c(x1,x2))
  # } else {
  #   stop("x1 must be list, or x1 and x2 must be vectors.")
  # }
  # h_max <- (a1-a0)
  # t_ <- seq(a0, a1, by = incr)
  # h_ <- seq(incr, h_max, by = incr)
  # if (perm == FALSE) {
  #   res <- matrix(NA, length(h_), length(t_))
  #   for (i in 1:length(t_)) {
  #     threshold <- t_[i]
  #     for (j in 1:length(h_)) {
  #       h <- h_[j]
  #       if (scaled_by == -1) {
  #         sc_by <- h
  #       } else {
  #         sc_by <- scaled_by
  #       }

  #       if ((threshold - h) >= min(t_) & (threshold + h) <= max(t_)) {
  #         x_11 <- x1[x1 > threshold -h & x1 < threshold ]
  #         x_21 <- x2[x2 > threshold -h & x2 < threshold ]
  #         x_12 <- x1[x1 >= threshold & x1 < threshold + h]
  #         x_22 <- x2[x2 >= threshold & x2 < threshold + h]

  #         if ((length(x_11) >= 1 & length(x_12) >= 1) & (length(x_21) >= 1 & length(x_22) >= 1)) {

  #           bef <- wass_R1(x_11, x_21, scaled, sc_by, by, perm = FALSE)
  #           aft <- wass_R1(x_12, x_22, scaled, sc_by, by, perm = FALSE)
  #           res[j,i] <- bef$W -  aft$W
  #         }
  #       }
  #     }
  #   }
  #   hh <- rep(h_, length(t_))
  #   tt <- rep(t_, each=length(h_))
  #   res <- data.frame(t = tt, h = hh, theta = c(res))
  #   out <- list(theta = res)
  #   class(out) <- c("th_idt", "list")
  #   return(out)

  # } else {
  #   res <- matrix(NA, length(h_), length(t_))
  #   perm <- matrix(NA, length(h_), length(t_))

  #   for (i in 1:length(t_)) {
  #     threshold <- t_[i]
  #     for (j in 1:length(h_)) {
  #       h <- h_[j]
  #       if (scaled_by == -1) {
  #         sc_by <- h
  #       } else {
  #         sc_by <- scaled_by
  #       }
  #         if ((threshold - h) >= min(t_) & (threshold + h) <= max(t_)) {
  #           x_11 <- x1[x1 > threshold -h & x1 < threshold ]
  #           x_21 <- x2[x2 > threshold -h & x2 < threshold ]
  #           x_12 <- x1[x1 > threshold & x1 < threshold + h]
  #           x_22 <- x2[x2 > threshold & x2 < threshold + h]

  #           if ((length(x_11) >= 1 & length(x_12) >= 1) & (length(x_21) >= 1 & length(x_22) >= 1)) {

  #             bef <- wass_R1(x_11, x_21, scaled, sc_by, by, perm = TRUE, n_perm)
  #             aft <- wass_R1(x_12, x_22, scaled, sc_by, by, perm = TRUE, n_perm)

  #             d_ <- bef$W - aft$W
  #             res[j,i] <- d_
  #             d_perm <- aft$W_perm
  #             perm[j,i] <- mean(d_perm > d_, na.rm = TRUE)
            
  #           }
  #       }
  #     }
  #   }

    # hh <- rep(h_, length(t_))
    # tt <- rep(t_, each=length(h_))
    # out1 <- data.frame(t = tt, h = hh, theta = c(res))
    # out2 <- data.frame(t = tt, h = hh, pr_theta = c(perm))

    # out <- list(theta = out1, pr_theta = out2)

    # class(out) <- c("th_idt", "list")
    # return(out)
  #}
}

#' @rdname th_idt_wass
#' @export
th_idt_wass.list <- function(x1, x2 = NULL, stat = FALSE, scaled = FALSE, scaled_by = 0, quad_by = 1e-3, incr = 0.5, perm = FALSE, n_perm = 10^5) {
  y1 <- x1[[1]]
  y2 <- x2[[2]]
  if (!is.null(x2)) {
    message("x1 is list; ignoring x2 input")
  } 
  th_idt_wass.numeric(y1, y2, stat, scaled, scaled_by, by, incr, perm, n_perm)
}


















###################



#' Integration in Distribution over Time (KS)
#'
#' For two one-dimensional temporal point processes \eqn{X_1} and \eqn{X_2}, this function estimates the condition of the two processes becoming "integrated in distribution over time" (i.d.t.) at a moment in time \eqn{t}. This is achieved by computing the Kolmogorov-Smirnov distance \eqn{D} between their empirical culumative distribution functions (e.c.d.f.), as a function \eqn{d(t,h)} that encompasses observations on the support \eqn{(t-h,t+h)}. For practical significance in the estimation of \eqn{d}, Monte Carlo permutation of samples is used to estimate the probability that the observed value of \eqn{d} is greater than the distribution of permuted samples (i.e., that observations are interchangable between \eqn{X_1} and \eqn{X_2}). 
#' 
#' The Rcpp-based function \code{\link[arkhaia]{ks_arkhaia}} relied upon for faster estimation of the KS distance than \code{\link[stats]{ks.test}}.
#' 
#' @param x1 Either a list of vectors or a vector containing observations of sample 1.
#' @param x2 If \code{x1} is vector, \code{x2} is vector containing observations of sample 2.
#' @param stat Whether to multiply the distance by \eqn{\sqrt{\frac{n_1 n_2}{n_1 + n_2}}}. Default is \code{FALSE}.
#' @param incr The increment at which to sample values of \eqn{t} and \eqn{h}. Default is 0.5. 
#' @param perm Whether to perform a Monte Carlo permutation across all samples in the input, which will return the probability of an integrating threshold at time \eqn{t}. Default is \code{FALSE}.
#' @param n_perm Number of Monte Carlo permutations to perform. Default is \code{1e5}.
#' @returns A list containing a data frame  two vectors \code{x1} and \code{x2}, the probability that Monte Carlo permutation of samples results in a lower distance than the one obtained. For a list of vectors input in \code{x1}, a matrix of probabilities pairwise between all vectors.
#' 
#' @importFrom Rdpack reprompt
#' @export
idt_ks <- function(x1, x2 = NULL, stat = FALSE, incr = 0.5, perm = FALSE, n_perm = 1e5) {
    UseMethod("idt_ks")
}

#' @rdname idt_ks
#' @export
idt_ks.numeric <- function(x1, x2 = NULL, stat = FALSE, incr = 0.5, perm = FALSE, n_perm = 1e5) {
  if (is.null(x2)) {
    stop("x1 and x2 must be vectors")
  }
  if (stat == FALSE) {
    res <- idt_ks_arma(x1,x2,incr)
  } else {
    res <- idt_ks_stat_arma(x1,x2,incr)
  }
  mat <- res[[1]]
  mat[mat > 8000] <- NA
  t_ <- c(res[[2]])
  h_ <- c(res[[3]])
  hh <- rep(h_, length(t_))
  tt <- rep(t_, each=length(h_))
  out1 <- data.frame(t = tt, h = hh, d = c(mat))
  if (perm == TRUE) {
    if (stat == FALSE) {
      res_perm <- idt_ks_perm_arma(x1,x2,incr,n_perm)
    } else {
      res_perm <- idt_ks_stat_perm_arma(x1,x2,incr,n_perm)
    }
    res_perm[res_perm > 8] <- NA
    res0 <- replicate(n_perm, mat, simplify="array")
    res_perm0 <- (res0 > res_perm)
    res2 <- apply(res_perm0, c(1,2), mean, na.rm = TRUE)
    out2 <- data.frame(t = tt, h = hh, pr_d = c(res2))
    out <- list(d = out1, pr_d = out2)
  } else {
    out <- list(d = out1)
  }
  class(out) <- c("idt", "list")
  return(out)
}

#' @rdname idt_ks
#' @export
idt_ks.list <- function(x1, x2 = NULL, stat = FALSE, incr = 0.5, perm = FALSE, n_perm = 10e5) {
  y1 <- x1[[1]]
  y2 <- x1[[2]]
  if (stat == FALSE) {
    res <- idt_ks_arma(y1,y2,incr)
  } else {
    res <- idt_ks_stat_arma(y1,y2,incr)
  }
  mat <- res[[1]]
  mat[mat > 8000] <- NA
  t_ <- c(res[[2]])
  h_ <- c(res[[3]])
  hh <- rep(h_, length(t_))
  tt <- rep(t_, each=length(h_))
  out1 <- data.frame(t = tt, h = hh, d = c(mat))
  if (perm == TRUE) {
    if (stat == FALSE) {
      res_perm <- idt_ks_perm_arma(y1,y2,incr,n_perm)
    } else {
      res_perm <- idt_ks_stat_perm_arma(y1,y2,incr,n_perm)
    }
    res_perm[res_perm > 8] <- NA
    res0 <- replicate(n_perm, mat, simplify="array")
    res_perm0 <- (res0 > res_perm)
    res2 <- apply(res_perm0, c(1,2), mean, na.rm = TRUE)
    out2 <- data.frame(t = tt, h = hh, pr_d = c(res2))
    out <- list(d = out1, pr_d = out2)
  } else {
    out <- list(d = out1)
  }
  class(out) <- c("idt", "list")
  return(out)
}






#' Integration in Distribution over Time (Wasserstein)
#'
#' For two one-dimensional point processes, this function computes their Wasserstein R1 distances as a function \eqn{d(t,h)}, which encompasses observations on the support \eqn{(t-h,t+h)}. Monte Carlo permutation of samples is used to estimate the probability that the obtain distances is of practical signficance.
#' 
#' @param x1 Either a list of vectors or a vector containing observations of sample 1.
#' @param x2 If \code{x1} is vector, \code{x2} is vector containing observations of sample 2.
#' @param stat Whether to multiply the distance by \eqn{\sqrt{\frac{n_1 n_2}{n_1 + n_2}}}. Default is \code{FALSE}.
#' @param scaled Whether to scale the Wasserstein R1 distance used in the course of estimation. Default is \code{FALSE}.
#' @param scaled_by The magnitude by which to scale the Wasserstein R1 distance, if \code{scaled} is \code{TRUE}. Entering \code{0} (the default) will use the range of the maximum and minimum of the pair of point processes, while entering \code{-1} will use the value of \eqn{h}.
#' @param quad_by The increment of sampling for performing quadrature. Default is \code{1e-3}.
#' @param incr The increment at which to sample values of \eqn{t} and \eqn{h}. Default is 0.5. 
#' @param perm Whether to perform a Monte Carlo permutation across all samples in the input, which will return the probability of an integrating threshold at time \eqn{t}. Default is \code{FALSE}.
#' @param n_perm Number of Monte Carlo permutations to perform. Default is \code{10^5}.
#' @returns For two vectors \code{x1} and \code{x2}, the probability that Monte Carlo permutation of samples results in a lower  empirical scaled Wassertein R1 distance than the one obtained. For a list of vectors input in \code{x1}, a matrix of probabilities pairwise between all vectors.
#' 
#' @importFrom Rdpack reprompt
#' @export
idt_wass <- function(x1, x2 = NULL, stat = FALSE, scaled = FALSE, scaled_by = 0, quad_by = 1e-3, incr = 0.5, perm = FALSE, n_perm = 10^5) {
    UseMethod("idt_wass")
}

#' @rdname idt_wass
#' @export
idt_wass.numeric <- function(x1, x2 = NULL, stat = FALSE, scaled = FALSE, scaled_by = 0, quad_by = 1e-3, incr = 0.5, perm = FALSE, n_perm = 10^5) {
  if (is.null(x2)) {
    stop("x1 and x2 must be vectors")
  }
  if (!(scaled_by == -1 | scaled_by >= 0 )) {
    stop("scaled_by must be -1, 0, or greater than 0.")
  }
  if (stat == FALSE) {
    res <- idt_wass_arma(x1,x2,scaled,scaled_by,quad_by,incr)
  } else {
    res <- idt_wass_stat_arma(x1,x2,scaled,scaled_by,quad_by,incr)
  }
  mat <- res[[1]]
  mat[mat > 8000] <- NA
  t_ <- c(res[[2]])
  h_ <- c(res[[3]])
  hh <- rep(h_, length(t_))
  tt <- rep(t_, each=length(h_))
  out1 <- data.frame(t = tt, h = hh, d = c(mat))

  if (perm == TRUE) {
    if (stat == FALSE) {
      res_perm <- idt_wass_perm_arma(x1,x2,scaled,scaled_by,quad_by,incr,n_perm)
    } else {
      res_perm <- idt_wass_stat_perm_arma(x1,x2,scaled,scaled_by,quad_by,incr,n_perm)
    }

    res_perm[res_perm > 8000] <- NA
    res0 <- replicate(n_perm, mat, simplify="array")
    res_perm0 <- (res0 > res_perm)
    res2 <- apply(res_perm0, c(1,2), mean, na.rm = TRUE)
    out2 <- data.frame(t = tt, h = hh, pr_d = c(res2))
    out <- list(d = out1, pr_d = out2)
  } else {
    out <- list(d = out1)
  }
  class(out) <- c("idt", "list")
  return(out)
}

#' @rdname idt_wass
#' @export
idt_wass.list <- function(x1, x2 = NULL, stat = FALSE, scaled = FALSE, scaled_by = 0, quad_by = 1e-3, incr = 0.5, perm = FALSE, n_perm = 10^5) {
  y1 <- x1[[1]]
  y2 <- x2[[2]]
  if (!is.null(x2)) {
    message("x1 is list; ignoring x2 input")
  } 
  idt_wass.numeric(y1, y2, stat, scaled, scaled_by, by, incr, perm, n_perm)
}
