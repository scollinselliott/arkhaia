#' Least Squares Spectral Analysis (LSSA)
#'
#' Peforms a simple least squares fitting to time indexed data of the form \eqn{x(t) = \beta_{0} + \sum_{i =1}^n \beta_{1i} \cos (2 \pi f t) + \beta_{2i} \sin(2 \pi f t)}, using a range of potential frequencies \insertCite{vanicek_approximate_1969,vanicek_further_1971}{arkhaia}. Intercept may be ommitted. 
#' 
#' @param x A data frame of two columns, the first containing time indices and the second containing values.
#' @param freqs A vector of frequencies to evaluate. By default a grid from 0.001 to 0.5 is tested at an interval of 0.005.
#' @param intercept Whether to include the intercept. Default is \code{TRUE}.
#' @param type Can be either \code{"frequency"} (the default) or \code{"period"}.
#' @returns The a data frame containing the power and residual sum of squares for each frequency, as well as coefficients.
#' 
#' @references
#'   \insertAllCited{}
#' 
#' @importFrom Rdpack reprompt
#' @export
LSSA <- function(x, freqs = seq(0.001, 0.5, by = 0.0005), intercept = TRUE, type = "frequency") {
    UseMethod("LSSA")
}

#' @rdname LSSA
#' @export
LSSA.matrix <- function(x, freqs = seq(0.001, 0.5, by = 0.0005), intercept = TRUE, type = "frequency") {
  if (intercept == TRUE) {
      intrcpt <- 1
  } else if (intercept == FALSE) {
      intrcpt <- 0
  }
  P <- LSSA_arma(x, freqs, intrcpt)

  if (type == "frequency") {
    out <- data.frame(freq = freqs, power = P[,1], rss = P[,2], coef_cos = P[,3], coef_sin = P[,4], coef_const = P[,3])
  } else if (type == "period") {
    out <- data.frame(period = 1/freqs, power = P[,1], rss = P[,2], coef_cos = P[,3], coef_sin = P[,4], coef_const = P[,3])
  }
  return(out)
}

#' @rdname LSSA
#' @export
LSSA.data.frame <- function(x, freqs = seq(0.001, 0.5, by = 0.0005), intercept = TRUE, type = "frequency") {
  x <- as.matrix(x)
  out <- LSSA.matrix(x, freqs, intercept, type) 
  return(out)
}





#' Least Squares Spectral Analysis via Lowest Frequency Iteration (LSSA-LFI)
#'
#' Peforms a simple least squares fitting to time indexed data of the form y = SIGMA B cos(2 pi f t) + C sin(2 pi f t) + D, using an input of frequencies; the lowest frequency peak (not the highest power frequency) is chosen for regression, up to a chosen number of iterations. Intercept may be ommitted. The lowest frequency is equivalent to the longest period.
#' 
#' @param dat A data frame of two columns, with the first column containing time indices and the second containing values.
#' @param n_iter The number of iterations to run. Default is 1.
#' @param intercept Whether to include the intercept. Default is \code{TRUE}.
#' @param AIC If \code{TRUE}, only the result that has yeilded the lowest AIC (Aikake Information Criterion) is given. Default is \code{FALSE}.
#' @returns A list containing:
#'   * A list of the coefficients for each iteration (the intercept is included in the first iteration).
#'   * A vector of the frequencies.
#'   * The residual sum of squares (RSS) after each iteration (decreasing).
#'   * The AIC upon each iteration. (If the paraemter AIC is \code{TRUE}, this will stop at the lowest AIC value produced by the frequencies tested).
#' 
#' @export
LSSA_LFI <- function(dat, n_iter = 1, intercept = TRUE, AIC = FALSE) {
    UseMethod("LSSA_LFI")
}

#' @rdname LSSA_LFI
#' @export
LSSA_LFI.matrix <- function(dat, n_iter = 1, intercept = TRUE, AIC = FALSE) {
  coefs <- list()
  freqs <- numeric(n_iter)
  rss <- numeric(n_iter)
  aic <- numeric(n_iter)

  first <- LSSA(dat, type = "period", intercept = intercept)
  idx <- which(diff(first$power) < 0)[1]
  freq1 <- 1/first[idx,1]
  epsilon2 <- first[idx,3]

  freqs[1] <- freq1
  rss[1] <- epsilon2
  coefs[[1]] <- first[idx, 4:6]
  
  dat <- as.matrix(dat)
  resid <- LSSA_resid_arma(dat, freq1, intercept = intercept)

  prev_freq <- freq1
  dat0 <- data.frame(dat[,1], resid)
  dat0 <- as.matrix(dat0)

  if (intercept == TRUE) {
    aic[1] <- 2 * (3 + 1) + log(epsilon2/nrow(dat)) * nrow(dat)
  } else {
    aic[1] <- 2 * (2 + 1) + log(epsilon2/nrow(dat)) * nrow(dat)
  }

  if (n_iter > 1) {
    for (j in 2:n_iter) {
      dat2 <- LSSA(dat0, intercept = FALSE, type = "period")
      dat2a <- dat2[dat2$period < 1/ prev_freq,]
      idx <- which(diff(dat2a$power) < 0)[1]
      freq_ <- 1/dat2a[idx,1]
      epsilon2 <- dat2a[idx,3]
      
      rss[j] <- epsilon2
      coefs[[j]] <- dat2a[idx, 4:5]
      freqs[j] <- freq_

      resid <- LSSA_resid_arma(dat0, freq_, intercept = FALSE)
      dat0 <- data.frame(dat[,1], resid)
      dat0 <- as.matrix(dat0)

      prev_freq <- freq_

      if (intercept == TRUE) {
          aic[j] <- 2 * (3 * j + 1) + log(epsilon2/nrow(dat)) * nrow(dat)
      } else {
          aic[j] <- 2 * (3 * j) + log(epsilon2/nrow(dat)) * nrow(dat)
      }
      
    }
  }

  if (AIC == FALSE) {
    out <- list(coefs = coefs, freqs = freqs, rss = rss, aic = aic)
  } else {
    idx <- which.min(aic)
    out <- list(coefs = coefs[1:idx], freqs = freqs[1:idx], rss = rss[1:idx], aic = aic[1:idx])
  }
  class(out) <- c("LSSA_LFR", "list")

  return(out)
}

#' @rdname LSSA_LFI
#' @export
LSSA_LFI.data.frame <- function(dat, n_iter = 1, intercept = TRUE, AIC = FALSE) {
  dat <- as.matrix(dat)
  out <- LSSA_LFI.matrix(dat, n_iter, intercept, AIC) 
  return(out)
}









#' Linear Dependence of LSSA-LFI Candidate Models via AIC
#'
#' For a set of time series (namely partitions of set of series) contained in a \code{list}, will compute the Akaike Information Criterion (AIC) for each candidate set.  
#' 
#' @param x A \code{list} containing the time series, each of which should be a matrix or data frame with time index in the first column and value in the second.
#' @param sets Candidate sets to evaluate; must be a \code{list} of \code{lists} containing the indices of the sets in \code{x}. If left \code{NULL}, two sets are evaluated: all series pooled together [1] and all series kept separate [2].
#' @param n_iter The number of iterations to run for the least squares spectral analysis via lowest frequency iteration (LSSA-LFI). Default is 1.
#' @param intercept Whether to include the intercept in the least squares spectral analysis via lowest frequency iteration (LSSA-LFI). Default is \code{TRUE}.
#' @returns A \code{list} containing the AIC for each candidate set, for each iteration.
#' 
#' @export
LSSA_LFI_candidates <- function(x, sets = NULL, n_iter = 1, intercept = TRUE) {
    UseMethod("LSSA_LFI_candidates")
}

#' @rdname LSSA_LFI_candidates
#' @export
LSSA_LFI_candidates.list <- function(x, sets = NULL, n_iter = 1, intercept = TRUE) {
  out <- list()

  n <- length(x)
  all <- data.frame(t_ = c(), y = c())

  for (i in 1:n) {
    all <- rbind(all, x[[i]])
  }

  n_obs <- nrow(all)

  # if sets is null; first set is pooled, second separate
  if ( is.null(sets) ) {
    sets <- list()
    sets[[1]] <- list(1:n) 
    set2 <- list()
    for (i in 1:n) {
      set2[[i]] <- i
    }
    sets[[2]] <- set2
  }

  # loop throught the partitions
  for (P in 1:length(sets)) {

    # number of partitions (sets of terms)
    n_m <- length(sets[[P]])

    epsilon <- numeric(n_iter)
    # loop through each set in the partition
    for (j in 1:length(sets[[P]])) {

      # construct the pooled observations for each set
      X <- c()
      for (jj in 1:length(sets[[P]][[j]])) {
        idx <- sets[[P]][[j]][[jj]]
        X <- rbind(X, x[[ idx ]])
      }

      W <- LSSA_LFI(X, n_iter = n_iter, intercept = intercept, AIC = FALSE)

      epsilon <- epsilon + W$rss
    }

    AIC_ <- 2 * n_m * (1:n_iter * 3 + 1) + n_obs * log(epsilon / n_obs)
    out[[P]] <- AIC_
  }
  
  class(out) <- c("LSSA_LFI_AIC", "list")

  return(out)
}





#' Trim to Epoch
#'
#' For a data frame in which the first column contains a time index and the second colum observations, trim the data frame to include observations only with a given epoch (time period).
#' 
#' @param x A \code{data frame} or \code{matrix} containing time-indexed data, with time index in the first column and value in the second.
#' @param epoch A numeric vector given the start and end time indices of the epoch.
#' @returns A \code{data frame} containing only those observations with time indices within the epoch.
#' 
#' @export
trim_epoch <- function(x, epoch = NULL) {
    UseMethod("trim_epoch")
}

#' @rdname trim_epoch
#' @export
trim_epoch.matrix <- function(x, epoch = NULL) {
  if (!(is.numeric(epoch) & length(epoch) ==2)  ) {
    stop("Epoch must a numeric vector of length 2 indicating start and end.")
  }
  if (epoch[2] < epoch[1]) {
    stop("Start of epoch must be earlier than end.")
  }

  out <- x[x[,1] > epoch[1] & x[,1] < epoch[2]  ,]

  class(out) <- "data.frame"
  return(out)
}

#' @rdname trim_epoch
#' @export
trim_epoch.data.frame <- function(x, epoch = NULL) {
  x <- as.matrix(x)
  trim_epoch.matrix(x, epoch)
}





#' Model Selection of LSSA-LFI Candidates
#'
#' For an \code{LSSA_LFI_AIC} object see \code{\link[arkhaia]{LSSA_LFI_candidates}}), returns the index of the candidate model. If \code{sets} is \code{NULL} in the \code{LSSA_LFI_candidates} function, an index of [1] refers to the model of a pooled (joint) grouping, [2] refers to the model of discrete, separate groupings. 
#' 
#' @param x An \code{LSSA_LFI_AIC} object. 
#' @returns The index of the grouping which contains the lowest AIC score.
#' 
#' @export
model_select <- function(x) {
    UseMethod("model_select")
}

#' @rdname model_select
#' @export
model_select.LSSA_LFI_AIC <- function(x) {
  x_min <- numeric(length(x))
  for (j in 1:length(x)) {
    y <- x[[j]]
    x_min[j] <- min(y, na.rm = TRUE)
  }
  out <- which.min(x_min)
  return(out)
}




#' LSSA-LFI Model
#'
#' Generates a data frame of values \eqn{f(t)} of the model generated by LSSA-LFI (see \code{\link[arkhaia]{LSSA_LFI}})).
#' 
#' @param dat A data frame where time indices are in the first column and values are in the second.
#' @param t_  A vector giving samples range of \eqn{t} for computing \eqn{f(t)}. Default is from the minimum to maximum time index sampled at 0.01 intervals.
#' @param n_iter The number of iterations to run. Default is 1.
#' @param intercept Whether to include the intercept in the least squares spectral analysis via lowest frequency iteration (LSSA-LFI). Default is \code{TRUE}.
#' @param label Default is \code{"model"}.
#' @returns A data frame containing \eqn{t, f(t)}.
#' 
#' @export
LSSA_LFI_model <- function(dat, t_ = NULL, n_iter = 1, intercept = TRUE, label = "model") {
    UseMethod("LSSA_LFI_model")
}

#' @rdname LSSA_LFI_model
#' @export
LSSA_LFI_model.matrix <- function(dat, t_ = NULL, n_iter = 1, intercept = TRUE, label = "model") {  
  dat_LSSA <- LSSA_LFI(dat[,1:2], n_iter = n_iter, intercept = intercept)

  if (is.null(t_)) {
    t_ <- seq(min(dat[,1]), max(dat[,1]), by = 0.01)
  }

  theta <- t_ * 2 * pi * dat_LSSA$freqs[1]
  if (intercept == TRUE) {
    Z <- dat_LSSA$coefs[[1]]$coef_cos * cos(theta) + dat_LSSA$coefs[[1]]$coef_sin * sin(theta) + dat_LSSA$coefs[[1]]$coef_const
  } else if (intercept == FALSE) {
    Z <- dat_LSSA$coefs[[1]]$coef_cos * cos(theta) + dat_LSSA$coefs[[1]]$coef_sin * sin(theta)
  }

  if (n_iter > 1) {
    for (i in 2:length(dat_LSSA$freqs)) {
      theta <- t_ * 2 * pi * dat_LSSA$freqs[i]
      Z <- Z + dat_LSSA$coefs[[i]]$coef_cos * cos(theta) + dat_LSSA$coefs[[i]]$coef_sin * sin(theta) 
    }
  }

  out <- data.frame(t = t_, y = Z, label =  rep(label, length(t_)))
  return(out)
}

#' @rdname LSSA_LFI_model
#' @export
LSSA_LFI_model.data.frame <- function(dat, t_ = NULL, n_iter = 1, intercept = TRUE, label = "model") {  
  dat <- as.matrix(dat[,1:2])
  LSSA_LFI_model.matrix(dat, t_ = t_, n_iter = n_iter, intercept = intercept, label = label)
}







#' Validated Linear Dependence via LSSA-LFI
#'
#' Probability of linear dependence between two groups of time series observations using a LSSA-LFI model selection (see \code{\link[arkhaia]{LSSA_LFI_candidates}})), given a list of at least three time series. Confounding variate is selected from the remaining time series in the list.
#' 
#' @param x A list of data frames.
#' @param pair The pair of commodities to evaluate in the list \code{x}, either names or indices.
#' @param n_iter The number of iterations to run. Default is 1.
#' @param intercept Whether to include the intercept in the least squares spectral analysis via lowest frequency iteration (LSSA-LFI). Default is \code{TRUE}.

#' @returns An upper-triangular matrix, in which 1 indicates linear dependence between two variates and 0 indicates independence.
#' 
#' @export
LSSA_LFI_validated <- function(x, pair = NULL, n_iter = 1, intercept = TRUE) {
    UseMethod("LSSA_LFI_validated")
}

#' @rdname LSSA_LFI_validated
#' @export
LSSA_LFI_validated.list <- function(x, pair = NULL, n_iter = 1, intercept = TRUE) {
  if (is.null(pair)) {
    stop("Must supply the pair of data frames in x to be evaluated.")
  }
  if (typeof(pair) == "character") {
    pair <- which(names(x) %in% pair)
  }
  nx <- 1:length(x)
  conf <- nx[-pair]
  res <- numeric(length(conf))

  partition <- list(list(c(1,2,3)),
                    list(c(1,2), c(3)),
                    list(c(1,3), c(2)),
                    list(c(2,3), c(1)),
                    list(c(1),c(2), c(3)))

  for (i in 1:length(conf)) {
    idx <- c(pair, conf[i])
    lssa_i <- LSSA_LFI_candidates(x[idx], sets = partition, n_iter = n_iter, intercept = intercept )
    mod <- model_select(lssa_i)
    if (mod %in% c(1,2)) {
      res[i] <- 1
    } else {
      res[i] <- 0
    }
  }
  return(mean(res))
}




#' Pairwise Selection of Linearlly Dependent LSSA-LFI Models
#'
#' Evaluate pairwise linear dependence between observations using a LSSA-LFI valdiated model selection (see \code{\link[arkhaia]{LSSA_LFI_validated}})). 
#' 
#' @param x A list of data frames.
#' @param n_iter The number of iterations to run. Default is 1.
#' @param intercept Whether to include the intercept in the least squares spectral analysis via lowest frequency iteration (LSSA-LFI). Default is \code{TRUE}.

#' @returns An upper-triangular matrix, containing the probability of linear dependence between series.
#' 
#' @export
LSSA_LFI_pairwise <- function(x, n_iter = 1, intercept = TRUE) {
    UseMethod("LSSA_LFI_pairwise")
}

#' @rdname LSSA_LFI_pairwise
#' @export
LSSA_LFI_pairwise.list <- function(x, n_iter = 1, intercept = TRUE) {
  idx_list <- list()
  mat <- matrix(NA, nrow = length(x), ncol = length(x))
  for (i in 1:(length(x)-1)) {
    for (j in (i+1):length(x)) {
      pair_ <- c(i,j)
      mat[i,j] <- LSSA_LFI_validated(x, pair = pair_, n_iter = n_iter, intercept = intercept)
    }
  }
  rownames(mat) <- names(x)
  colnames(mat) <- names(x)
  return(mat)
}







#' Period-Variable Pairwise Selection of Linearlly Dependent LSSA-LFI Models
#'
#' Evaluate pairwise linear dependence between two time series using a LSSA-LFI valdiated model selection (see \code{\link[arkhaia]{LSSA_LFI_validated}})), with variable length time period. 
#' 
#' @param x A list of data frames.
#' @param pair The pair of commodities to evaluate in the list \code{x}, either names or indices.
#' @param n_iter The number of iterations to run. Default is 1.
#' @param intercept Whether to include the intercept in the least squares spectral analysis via lowest frequency iteration (LSSA-LFI). Default is \code{TRUE}.
#' @param t_range The range of the time period. Default are the minimum and maximum dates spanned by the data in \code{x}.
#' @param h_range The range of the potential windows of \eqn{h}. Default is from 0 to entire span of \code{t_range}. 
#' @param t_grid The length of interval along which to sample \eqn{t}. Default is 1.
#' @param h_grid The length of interval along which to sample \eqn{h}. Default is 1.
#' @returns A matrix giving the probability of linear dependence, predicated upon time \eqn{t} and window length \eqn{h}.
#'
#' @export
LSSA_LFI_epoch <- function(x, n_iter = 1, intercept = TRUE) {
    UseMethod("LSSA_LFI_epoch")
}

#' @rdname LSSA_LFI_epoch
#' @export
LSSA_LFI_epoch <- function(x, pair = NULL, n_iter = 1, intercept = TRUE, t_range = NULL, h_range = NULL, t_grid = 1, h_grid = 1) {
  if (is.null(pair)) {
    stop("Must supply the pair of data frames in x to be evaluated.")
  }
  all <- data.frame()
  for (i in 1:length(x)) {
    all <- rbind(all, x[[i]])
  }
  start <- min(all$t)
  end <- max(all$t)
  if (is.null(t_range)) {
    t_ <- seq(start, end, by = t_grid)
  } else {
    t_ <- seq(t_range[1], t_range[2], by = t_grid)
  }
  if (is.null(h_range)) {
    h_ <- seq(h_grid, (end - start) + h_grid, by = h_grid)
  } else {
    h_ <- seq(h_range[1], h_range[2], by = h_grid)
  }
  out <- matrix(NA, nrow = length(t_), ncol = length(h_))

  for (ti in 1:length(t_)) {
    for (hi in 1:length(h_)) {
      tmp <- list()
      do_lssa <- TRUE
      n_obs_ <- c()
      for (i in 1:length(x)) {
        trimmed <- trim_epoch(x[[i]], epoch = c(t_[ti], t_[ti] + h_[hi]))
        n_obs <- nrow(trimmed)
        n_obs_ <- c(n_obs_, n_obs)
        if (n_obs < 4) {
          do_lssa <- FALSE
        }
        tmp[[names(x)[i]]] <- trimmed
      } 
      if (do_lssa == TRUE) {
        n_iter_ <- min(c(n_obs_, n_iter))

        out[(length(h_) - hi)+1,ti] <- LSSA_LFI_validated(tmp, pair = pair, n_iter = n_iter_, intercept = intercept)
      }
    }
    cat("\n Percent complete: ", round(ti/length(t_), 2)*100, "%")
  }



  return(out)

}



