#' Least Squares Spectral Analysis (LSSA)
#'
#' Performs a simple least squares fitting to time indexed data of the form \eqn{x(t) = \beta_{0} + \sum_{i =1}^n \beta_{1i} \cos (2 \pi f t) + \beta_{2i} \sin(2 \pi f t)}, using a range of potential frequencies \insertCite{vanicek_approximate_1969,vanicek_further_1971}{arkhaia}. Intercept may be ommitted. 
#' 
#' @param x A time_series object (see \code{\link[arkhaia]{time_series}}) with only one variate.
#' @param freqs A vector of frequencies to evaluate. By default a grid from 0.001 to 0.5 is tested at an interval of 0.005.
#' @param intercept Whether to include the intercept. Default is \code{TRUE}.
#' @param type Type of output. Can be either \code{"frequency"} (the default) or \code{"period"}.
#' 
#' @examples
#' t_ <- runif(100, 0, 100)
#' y <- 5 * sin(0.1 * t_) - 2.3 * cos(0.3 * t_ + 1)
#' dat <- time_series(t = t_, x = y)
#' LSSA(dat)
#' 
#' @returns A data frame containing the power and residual sum of squares for each frequency, as well as coefficients.
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
LSSA.time_series <- function(x, freqs = seq(0.001, 0.5, by = 0.0005), intercept = TRUE, type = "frequency") {
  if (nlevels(x[,3]) > 1) {
    warning("More than one variate in time series. Pooling observations.")
  }

  if (intercept == TRUE) {
      intrcpt <- 1
  } else if (intercept == FALSE) {
      intrcpt <- 0
  }
  P <- LSSA_arma(as.matrix(x[,1:2]), freqs, intrcpt)

  if (type == "frequency") {
    out <- data.frame(freq = freqs, power = P[,1], rss = P[,2], coef_cos = P[,3], coef_sin = P[,4], coef_const = P[,5])
  } else if (type == "period") {
    out <- data.frame(period = 1/freqs, power = P[,1], rss = P[,2], coef_cos = P[,3], coef_sin = P[,4], coef_const = P[,5])
  }
  return(out)
}

# #' @rdname LSSA
# #' @export
# LSSA.data.frame <- function(x, freqs = seq(0.001, 0.5, by = 0.0005), intercept = TRUE, type = "frequency") {
#   x <- as.matrix(x)
#   out <- LSSA.matrix(x, freqs, intercept, type) 
#   return(out)
# }





#' Least Squares Spectral Analysis via Lowest Frequency Iteration (LSSA-LFI)
#'
#' Performs a simple least squares fitting to time indexed data of the form  \eqn{x(t) = \beta_{0} + \sum_{i =1}^n \beta_{1i} \cos (2 \pi f t) + \beta_{2i} \sin(2 \pi f t)}, using an input of frequencies; the lowest frequency peak (not the highest power frequency) is chosen for regression, up to a chosen number of iterations. Intercept may be ommitted. The lowest frequency is equivalent to the longest period.
#' 
#' @param x A time_series object (see \code{\link[arkhaia]{time_series}}) with only one variate.
#' @param n_iter The number of iterations to run. Default is 1.
#' @param intercept Whether to include the intercept. Default is \code{TRUE}.
#' @param AIC If \code{TRUE}, only the result that has yeilded the lowest AIC is given \insertCite{akaike_information_1973}{arkhaia}. Default is \code{FALSE}.
#' @returns A list containing:
#' \itemize{
#'  \item A list of the coefficients for each iteration (the intercept is included in the first iteration).
#'  \item A vector of the frequencies.
#'  \item The residual sum of squares (RSS) after each iteration (decreasing).
#'  \item The AIC upon each iteration. (If the parameter AIC is \code{TRUE}, this will stop at the lowest AIC value produced by the frequencies tested).
#' }
#' 
#' @references
#'   \insertAllCited{}
#'
#' @importFrom Rdpack reprompt
#' @export
LSSA_LFI <- function(x, n_iter = 1, intercept = TRUE, AIC = FALSE) {
    UseMethod("LSSA_LFI")
}

#' @rdname LSSA_LFI
#' @export
LSSA_LFI.time_series <- function(x, n_iter = 1, intercept = TRUE, AIC = FALSE) {
  if (nlevels(x[,3]) > 1) {
    warning("more than one variate in time series. Pooling observations.")
  }

  if (intercept == TRUE) {
    intrcpt = 1
  } else {
    intrcpt = 0
  }
  out0 <- LSSA_LFI_arma(as.matrix(x[,1:2]), n_iter, intrcpt)

  out <- list(coefs = out0[,1:3], freqs = out0[,4], rss = out0[,5], aic = out0[,6])

  if (AIC == TRUE) {
    idx <- which.min(out$aic)
    out <- list(coefs = out0[1:idx,1:3], freqs = out0[1:idx, 4], rss = out0[1:idx, 5], aic = out0[1:idx, 6])
  }
  class(out) <- c("LSSA_LFI", "list")

  return(out)
}

# #' @rdname LSSA_LFI
# #' @export
# LSSA_LFI.data.frame <- function(x, n_iter = 1, intercept = TRUE, AIC = FALSE) {
#   x <- as.matrix(x)
#   out <- LSSA_LFI.matrix(x, n_iter, intercept, AIC) 
#   return(out)
# }




#' Linear Dependence of LSSA-LFI Candidate Models via AIC
#'
#' For a set of time series (namely partitions of set of series) contained in a list, will compute the Akaike Information Criterion (AIC) for each candidate set \insertCite{akaike_information_1973}{arkhaia}.  
#' 
#' @param x A time_series object (see \code{\link[arkhaia]{time_series}}) with more than one variate.
#' @param sets Candidate sets of variates to evaluate; must be a \code{list} of \code{lists} containing the indices of the levels  in \code{x$variate}. If left \code{NULL}, two sets are evaluated: all series pooled together [1] and all series kept separate [2].
#' @param n_iter The number of iterations to run for the least squares spectral analysis via lowest frequency iteration (LSSA-LFI). Default is 1.
#' @param intercept Whether to include the intercept in the least squares spectral analysis via lowest frequency iteration (LSSA-LFI). Default is \code{TRUE}.
#' @returns The index of the set yielding the lowest AIC. If \code{sets} is \code{NULL}, the output of \code{[1]} indicates linear dependence; if the output is \code{[2]}, linear independence.
#'  
#' @references
#'   \insertAllCited{}
#'
#' @importFrom Rdpack reprompt
#' @export
LSSA_LFI_candidates <- function(x, sets = NULL, n_iter = 1, intercept = TRUE) {
    UseMethod("LSSA_LFI_candidates")
}

#' @rdname LSSA_LFI_candidates
#' @export
LSSA_LFI_candidates.time_series <- function(x, sets = NULL, n_iter = 1, intercept = TRUE) {
  if (nlevels(x$variate) < 2) {
    stop("time_series input must contain more than one variate.")
  }
  n <- nlevels(x$variate)
  # y <- as.matrix(do.call(rbind, x))
  x_rows <- table(x$variate) # unlist(lapply(x, nrow))
  idx <- c(0, cumsum(x_rows)[1:(n-1)]) 

  # if sets is null; first set is pooled, second separate
  if ( is.null(sets) ) {
    pld <- rep(0, n)
    tog <- 0:(n-1)
    sets0 <- matrix(c(pld, tog), nrow = 2, byrow = TRUE)
  } else {
    lev <- levels(x$variate)
    for (i in 1:n) {
      sets1 <- unlist(sets)
      sets1[sets1 == lev[i]] <- i
      sets <- relist(sets1, sets)
    }
    sets1 <- unlist(sets)
    sets1 <- as.numeric(sets1)
    sets <- relist(sets1, sets)
    sets0 <- partition_list_to_matrix_arma(sets)
  }

  if (intercept == TRUE) {
    intrcpt <- 1
  } else {
    intrcpt <- 0
  }

  # add 1 for R indexing
  out <- LSSA_LFI_candidates_arma(as.matrix(x[,1:2]), idx, x_rows, sets0, n_iter, intrcpt) + 1
  
  return(out)
}



#' Trim to Epoch
#'
#' For a data frame in which the first column contains a time index and the second colum observations, trim the data frame to include observations only with a given epoch (time period).
#' 
#' @param x A time_series object (see \code{\link[arkhaia]{time_series}}).
#' @param epoch A numeric vector given the start and end time indices of the epoch.
#' @returns A data frame containing only those observations with time indices within the epoch.
#' 
#' @export
trim_epoch <- function(x, epoch = NULL) {
    UseMethod("trim_epoch")
}

#' @rdname trim_epoch
#' @export
trim_epoch.time_series <- function(x, epoch = NULL) {
  if (!(is.numeric(epoch) & length(epoch) == 2)  ) {
    stop("epoch must a numeric vector of length 2 indicating start and end.")
  }
  if (epoch[2] < epoch[1]) {
    stop("start of epoch must be earlier than end.")
  }

  out <- x[x[,1] > epoch[1] & x[,1] < epoch[2]  ,]
  out <- as.data.frame(out)
  class(out) <- c("time_series", "data.frame")
  return(out)
}





#' Repartition (Regroup) Variates
#'
#' Given multiple time series in a \code{\link[arkhaia]{time_series}} object (\code{x}), group the variates into a new set of variates, expressed as a list of the variate names (\code{set}). Returns a \code{time_series} object using the new partitioning of variates.
#' 
#' @param x A \code{time_series} object containing multiple  time series, each of which should be a matrix or data frame with time index in the first column and value in the second.
#' @param set A \code{list} of vectors containing the names of variates in each of the new groups.
#' 
#' @examples 
#' idx <- runif(40, 0, 100)
#' x1 <- runif(10, 0, 1) + 3
#' x2 <- runif(10, 0, 1) + 2
#' x3 <- runif(10, 0, 1) + 1
#' x4 <- runif(10, 0, 1)
#' x <- time_series(t = idx, x = c(x1,x2,x3,x4),
#'   variate = c(rep("x1", 10), rep("x2", 10), 
#'   rep("x3", 10), rep("x4", 10)))
#' 
#' repartition(x, set = list(c("x1", "x4"), "x2", "x3"))
#' 
#' @returns A \code{list} containing the data according to the partition structure described by \code{set}.
#' 
#' @export
repartition <- function(x, set = NULL) {
    UseMethod("repartition")
}

#' @rdname repartition
#' @export
repartition.time_series <- function(x, set = NULL) {
  if (!is.list(set)) {
    stop("repartitioning set must be provided.")
  }
  vars <- unlist(set)
  if (TRUE %in% is.na(match(vars, levels(x$variate)))) {
    stop("set contains variates not in x")
  }
  if (TRUE %in% is.na(match(levels(x$variate), vars))) {
    stop("variates in x are not designated in set")
  }
  x$variate <- as.character(x$variate)
  for (i in 1:length(set)) {
    x[x$variate %in% set[[i]], 3] <- i
  }
  out <- x
  out$variate <- factor(out$variate)
  class(out) <- c("time_series", "data.frame")
  return(out)
}





# #' Model Selection of LSSA-LFI Candidates
# #'
# #' For an \code{LSSA_LFI_AIC} object (see \code{\link[arkhaia]{LSSA_LFI_candidates}}), returns the index of the candidate model. If \code{sets} is \code{NULL} in the \code{LSSA_LFI_candidates} function, an index of [1] refers to the model of a pooled (joint) grouping, [2] refers to the model of discrete, separate groupings. 
# #' 
# #' @param x An \code{LSSA_LFI_AIC} object. 
# #' @returns The index of the grouping which contains the lowest AIC score.
# #' 
# #' @export
# model_select <- function(x) {
#     UseMethod("model_select")
# }

# #' @rdname model_select
# #' @export
# model_select.LSSA_LFI_AIC <- function(x) {
#   x_min <- numeric(length(x))
#   for (j in 1:length(x)) {
#     y <- x[[j]]
#     x_min[j] <- min(y, na.rm = TRUE)
#   }
#   out <- which.min(x_min)
#   return(out)
# }




#' LSSA-LFI Model
#'
#' Generates a data frame of values \eqn{x(t)} of the model generated by LSSA-LFI (see \code{\link[arkhaia]{LSSA_LFI}}).
#' 
#' @param x A time_series object (see \code{\link[arkhaia]{time_series}}) with only one variate.
#' @param t_  A vector giving samples range of \eqn{t} for computing \eqn{x(t)}. Default is from the minimum to maximum time index sampled at 0.01 intervals.
#' @param n_iter The number of iterations to run. Default is 1.
#' @param intercept Whether to include the intercept in the least squares spectral analysis via lowest frequency iteration (LSSA-LFI). Default is \code{TRUE}.
#' @param label Default is \code{"model"}.
#' @returns A data frame containing \eqn{t, x(t)}.
#' 
#' @export
LSSA_LFI_model <- function(x, t_ = NULL, n_iter = 1, intercept = TRUE, label = "model") {
    UseMethod("LSSA_LFI_model")
}

#' @rdname LSSA_LFI_model
#' @export
LSSA_LFI_model.time_series <- function(x, t_ = NULL, n_iter = 1, intercept = TRUE, label = "model") {  
  if (intercept == TRUE) {
    intrcpt <- 1
  } else {
    intrcpt <- 0
  }

  if (is.null(t_)) {
    t_ <- seq(min(x[,1]), max(x[,1]), by = 0.01)
  }

  Y <- LSSA_LFI_model_arma(as.matrix(x[,1:2]), t_, n_iter, intrcpt) 

  out <- data.frame(t = Y[,1], x = Y[,2], variate =  rep(label, length(t_)))
  colnames(out) <- colnames(x)
  class(out) <- c("time_series", "data.frame")
  return(out)
}

# #' @rdname LSSA_LFI_model
# #' @export
# LSSA_LFI_model.data.frame <- function(x, t_ = NULL, n_iter = 1, intercept = TRUE, label = "model") {  
#   x <- as.matrix(x[,1:2])
#   LSSA_LFI_model.matrix(x, t_ = t_, n_iter = n_iter, intercept = intercept, label = label)
# }



#' Validated Linear Dependence via LSSA-LFI
#'
#' Probability of linear dependence between two groups of time series observations using a LSSA-LFI model selection (see \code{\link[arkhaia]{LSSA_LFI_candidates}}), on the basis of the inclusion of a third "attendant" variate. The third attendant variate is selected from the remaining time series in the list, hence the input list must include at least three times series.
#' 
#' @param x A time_series object (see \code{\link[arkhaia]{time_series}}) with at least three variates.
#' @param pair The pair of series to evaluate in the \code{x}, given levels in \code{x$variate}.
#' @param n_iter The number of iterations to run. Default is 1.
#' @param intercept Whether to include the intercept in the least squares spectral analysis via lowest frequency iteration (LSSA-LFI). Default is \code{TRUE}.
#'
#' @returns The probability of linear dependence between two time series, in which 1 indicates linear dependence almost surely and 0 indicates independence.
#' 
#' @export
LSSA_LFI_validated <- function(x, pair = NULL, n_iter = 1, intercept = TRUE) {
    UseMethod("LSSA_LFI_validated")
}

#' @rdname LSSA_LFI_validated
#' @export
LSSA_LFI_validated.time_series <- function(x, pair = NULL, n_iter = 1, intercept = TRUE) {
  if (nlevels(x$variate) < 3) {
    stop("must have at least three time series in x.")
  }
  if (is.null(pair)) {
    stop("must supply the pair of data frames in x to be evaluated.")
  }
  if (!is.character(pair)) {
    stop("invalid input type for pair.")
  }
  if (length(pair) != 2) {
    stop("pair must contain two elements")
  }
  if (NA %in% match(pair, levels(x$variate))) {
    stop("pair element not in variate")
  }

  base <- x[x$variate %in% pair,]
  base$variate <- as.numeric(factor(base$variate))
  lev <- levels(x$variate)

  idx_p <- match(pair, lev)
  attend <- lev[-idx_p]
  res <- numeric(length(attend))

  if (intercept == TRUE) {
    intrcpt <- 1
  } else {
    intrcpt <- 0
  }

  partition <-  matrix(c(0,0,0,
                         0,0,1,
                         0,1,0,
                         0,1,1,
                         0,1,2), ncol = 3, byrow = TRUE)

  for (i in 1:length(attend)) {
    y0 <- x[x$variate == attend[i], ]
    y0$variate <- rep(3, nrow(y0))
    y <- rbind(base, y0)
    y_rows <- table(y$variate)
    idx <- c(0, cumsum(y_rows)[1:2]) 

    mod <- LSSA_LFI_candidates_arma(as.matrix(y[,1:2]), idx, y_rows, partition, n_iter, intrcpt) + 1 # add 1 for R indexing
    if (mod %in% c(1,2)) {
      res[i] <- 1
    } # else {
    #  res[i] <- 0
    # }
  }
  return(mean(res))
}



#' Pairwise Selection of Linearly Dependent LSSA-LFI Models
#'
#' Evaluate pairwise linear dependence between observations using a LSSA-LFI validated model selection (see \code{\link[arkhaia]{LSSA_LFI_validated}}). 
#' 
#' @param x A time_series object (see \code{\link[arkhaia]{time_series}}) with at least three variates.
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
LSSA_LFI_pairwise.time_series <- function(x, n_iter = 1, intercept = TRUE) {
  n_factors <- nlevels(x$variate)
  if (n_factors < 3) {
    stop("need at least three variates to evaluate pairwise.")
  }
  lev <- levels(x$variate)
  mat <- matrix(NA, nrow = n_factors, ncol = n_factors)
  for (i in 1:(n_factors-1)) {
    for (j in (i+1):n_factors) {
      pair <- c(lev[i],lev[j])
      mat[i,j] <- LSSA_LFI_validated(x, pair, n_iter = n_iter, intercept = intercept)
    }
  }
  rownames(mat) <- levels(x$variate)
  colnames(mat) <- levels(x$variate)
  return(mat)
}



#' Period-Variable Pairwise Selection of Linearly Dependent LSSA-LFI Models
#'
#' Evaluate pairwise linear dependence between two time series using a LSSA-LFI valdiated model selection (see \code{\link[arkhaia]{LSSA_LFI_validated}})), with variable length time period. 
#' 
#' @param x A list of data frames.
#' @param pair The pair of commodities to evaluate in the list \code{x}, given as indices.
#' @param n_iter The number of iterations to run. Default is 1.
#' @param intercept Whether to include the intercept in the least squares spectral analysis via lowest frequency iteration (LSSA-LFI). Default is \code{TRUE}.
#' @param t_range The range of the time period. Default are the minimum and maximum dates spanned by the data in \code{x}.
#' @param h_range The range of the potential windows of \eqn{h}. Default is from 0 to entire span of \code{t_range}. 
#' @param t_grid The length of interval along which to sample \eqn{t}. Default is 1.
#' @param h_grid The length of interval along which to sample \eqn{h}. Default is 1.
#' @returns A matrix giving the probability of linear dependence, predicated upon time \eqn{t} and window length \eqn{h}.
#'
#' @export
LSSA_LFI_epoch <- function(x, pair = NULL, n_iter = 1, intercept = TRUE, t_range = NULL, h_range = NULL, t_grid = 1, h_grid = 1) {
    UseMethod("LSSA_LFI_epoch")
}

#' @rdname LSSA_LFI_epoch
#' @export
LSSA_LFI_epoch.list <- function(x, pair = NULL, n_iter = 1, intercept = TRUE, t_range = NULL, h_range = NULL, t_grid = 1, h_grid = 1) {
  if (is.null(pair)) {
    stop("must supply the pair of data frames in x to be evaluated.")
  }
  if (intercept == TRUE) {
    intrcpt <- 1
  } else {
    intrcpt <- 0
  }
  all <- data.frame()
  for (i in 1:length(x)) {
    all <- rbind(all, x[[i]])
  }
  start_ <- min(all$t)
  end_ <- max(all$t)
  if (is.null(t_range)) {
    t_ <- seq(start_, end_, by = t_grid)
  } else {
    t_ <- seq(t_range[1], t_range[2], by = t_grid)
  }
  if (is.null(h_range)) {
    h_ <- seq(h_grid, (end_ - start_) + h_grid, by = h_grid)
  } else {
    h_ <- seq(h_range[1], h_range[2], by = h_grid)
  }
  out <- matrix(NA, nrow = length(t_), ncol = length(h_))

  for (ti in 1:length(t_)) {
    for (hi in 1:length(h_)) {
      if (t_[ti] + h_[hi] <= end_) {
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
          # tmp[[names(x)[i]]] <- trimmed
          tmp[[i]] <- trimmed
        } 
        if (do_lssa == TRUE) {
          n_iter_ <- min(c(n_obs_, n_iter))

          out[(length(h_) - hi)+1,ti] <- LSSA_LFI_validated(tmp, pair, n_iter_, intrcpt)
        }
      }
    }
    message("Percent complete: ", round(ti/length(t_), 2)*100, "%")
  }

  colnames(out) <- t_
  rownames(out) <- rev(h_)
  return(out)

}





#' Partitioned LSSA-LFI Model Terms and AIC
#'
#' For multiple time series contained in a \code{\link[arkhaia]{time_series}} object, this function will compute the Akaike Information Criterion (AIC) as the number of terms (\code{n_iter}) is increased. This function can be performed after \code{\link[arkhaia]{LSSA_LFI_candidates}} and \code{\link[arkhaia]{repartition}} to determine the appropriate number of terms.
#' 
#' @param x A time_series object (see \code{\link[arkhaia]{time_series}}).
#' @param n_iter The number of iterations to run for the least squares spectral analysis via lowest frequency iteration (LSSA-LFI). Default is 1.
#' @param intercept Whether to include the intercept in the least squares spectral analysis via lowest frequency iteration (LSSA-LFI). Default is \code{TRUE}.
#' @param AIC If \code{TRUE}, only the result that has yeilded the lowest AIC is given \insertCite{akaike_information_1973}{arkhaia}. Default is \code{FALSE}.
#' 
#' @returns A data frame giving the number of iterations and the AIC for each. 
#'  
#' @references
#'   \insertAllCited{}
#'
#' @importFrom Rdpack reprompt
#' @export
LSSA_LFI_multi <- function(x, n_iter = 1, intercept = TRUE, AIC = FALSE) {
    UseMethod("LSSA_LFI_multi")
}

#' @rdname LSSA_LFI_multi
#' @export
LSSA_LFI_multi.time_series <- function(x, n_iter = 1, intercept = TRUE, AIC = FALSE) {
  n <- nlevels(x$variate)
  x_rows <- table(x$variate) # unlist(lapply(x, nrow))
  idx <- c(0, cumsum(x_rows)[1:(n-1)]) 
  # n <- length(x)
  # y <- as.matrix(do.call(rbind, x))
  # x_rows <- unlist(lapply(x, nrow))
  # idx <- c(0, cumsum(x_rows)[1:(n-1)]) 
  if (length(x_rows) == 1) {
    idx <- idx[1]
  }

  if (intercept == TRUE) {
    intrcpt <- 1
  } else {
    intrcpt <- 0
  }

  out <- LSSA_LFI_multi_arma(as.matrix(x[,1:2]), idx, x_rows, n_iter, intrcpt)

  out <- as.data.frame(out)
  colnames(out) <- c("n_iter", "rss", "aic")

  if (AIC == TRUE) {
    m <- which.min(out$aic)
    out <- as.numeric(out[m,])
    names(out) <- c("n_iter", "rss", "aic")
  }

  return(out)
}



#' Partitioned LSSA-LFI Model Terms and AIC: Composite Number of Parameters
#'
#' For multiple time series contained in a list, this function will compute the Akaike Information Criterion (AIC) for the best-fitting models particular to each series, including the case of an intercept-only model. This function can be performed after \code{\link[arkhaia]{LSSA_LFI_candidates}} and \code{\link[arkhaia]{repartition}} to determine the appropriate number of terms. For the situation where all time series must have the same number of terms, \code{\link[arkhaia]{LSSA_LFI_multi}} should be used, but this function should chosose a better fitting model, since it adjusts the number of parameters for each dataset.
#' 
#' @param x A time_series object (see \code{\link[arkhaia]{time_series}}).
#' @param intercept Whether to include the intercept in the least squares spectral analysis via lowest frequency iteration (LSSA-LFI). Default is \code{TRUE}.
#' 
#' @returns A list giving the number of parameters (\code{k}) and residual sum of squares (\code{rss}) for each dataset in \code{x}, along with the AIC for the paritioned model. 
#' 
#' @export
LSSA_LFI_comp <- function(x, intercept = TRUE) {
    UseMethod("LSSA_LFI_comp")
}

#' @rdname LSSA_LFI_comp
#' @export
LSSA_LFI_comp.time_series <- function(x, intercept = TRUE) {  
  # n <- length(x)
  # y <- as.matrix(do.call(rbind, x))
  # nx <- nrow(y)
  nx <- nrow(x)
  n <- nlevels(x$variate)
  x_rows <- table(x$variate) # unlist(lapply(x, nrow))
  idx <- c(0, cumsum(x_rows)[1:(n-1)]) 

  if (intercept == TRUE) {
    intrcpt <- 1
  } else {
    intrcpt <- 0
  }

  rss_i <- numeric(n)
  K <- numeric(n)

  lev <- levels(x$variate)

  for (i in 1:n) {
    xsub <- x[x$variate == lev[i], ]
    xsub$variate <- factor(xsub$variate)
    nxi <- nrow(xsub)
    lssa_lfi_i <- LSSA_LFI(xsub, n_iter = nxi-1, intercept = intrcpt, AIC = TRUE)
    k <- length(lssa_lfi_i$rss)
    rss0 <- sum((xsub[,2] - mean(xsub[,2]))^2)
    aic0 <- 2 * (1) + nxi * log(rss0 / nxi) + nxi + nxi * log(2 * pi)
    if (aic0 < lssa_lfi_i$aic[k]) {
      rss_i[i] <- rss0
      K[i] <- 0
    } else {
      rss_i[i] <- lssa_lfi_i$rss[k]
      K[i] <- k
    }
  }

  aic <- 2 * (3 *(sum(K)) + n) + nx * log(sum(rss_i) / nx) + nx + nx * log(2 * pi) 

  out <- list(k = K, rss = rss_i, aic = aic)

  return(out)
}



#' Decaying Trend-after-Event Model
#'
#' For multiple time-indexed observations and a number of events, \eqn{r = 1, ldots, n_r}, that occur at a given time index, this function performs a regression of the form \eqn{x(t) = \beta_0 + \sum_{r = 1}^{n_r} \gamma t I(t \leq t_r ) e^{-\lambda t} }, where \eqn{I(t \geq t_r )} is an indicator function equaling 1 if \eqn{t \geq t_r}. Estimation of \eqn{\gamma}, \eqn{\lambda}, and \eqn{\beta_0} is performed numerically on specified ranges of potential values. The intercept term is not used where cases where a zero intercept is expected (e.g., standardized observations).
#' 
#' @param x A data frame of two columns, with the first column containing time indices and the second containing values.
#' @param events A vector of the time indices of the events, \eqn{t_1, \ldots, t_{n_r}}.
#' @param gamma The range of values of \eqn{\gamma} to evaluate. Default is \code{c(-10,10)}.
#' @param lambda The range of values of \eqn{\lambda} to evaluate. Default is \code{c(0,5)}.
#' @param intercept The range of values of the intercept \eqn{\beta_0} to evaluate. If  \code{FALSE} is given instead (the default) no intercept term is included.
#' @param incr Increment by which to sample values of \code{gamma} and \code{lambda}. Default is \code{0.01}.
#' 
#' @examples
#' # generate example of positive trend (gamma = 2) 
#' # after event with slow decay (lambda = 0.05)
#' t_ <- runif(100, 0, 100)
#' events <- c(0, 54, 89)
#' gam <- 2
#' lam <- 0.05
#' x_ev1 <-  (t_ - events[1]) * gam  * exp(- lam * (t_ - events[1]))
#' x_ev2 <-  (t_ - events[2]) * gam  * exp(- lam * (t_ - events[2]))
#' x_ev3 <-  (t_ - events[3]) * gam  * exp(- lam * (t_ - events[3]))
#' x_ev1[t_ < events[1]] <- 0
#' x_ev2[t_ < events[2]] <- 0
#' x_ev3[t_ < events[3]] <- 0
#' y <- x_ev1 + x_ev2 + x_ev3 + rnorm(length(t_), 0, 0.25)
#' dat <- time_series(t = t_, x = y)
#' 
#' # fit model
#' res <- dtem(dat, events, gamma = c(-2.5,2.5), lam = c(0,0.1))
#' 
#' # coefficients
#' res$coefs
#' 
#' @returns A list containing
#' \itemize{
#'  \item A vector of coefficients and the residual sum of squares.
#'  \item A matrix of the residual sum of squares for samples within the ranges of gamma and lambda.
#' }
#' 
#' @references
#'   \insertAllCited{}
#' 
#' @importFrom Rdpack reprompt
#' @export
dtem <- function(x, events, gamma = c(-10,10), lambda = c(0,5), intercept = FALSE, incr = 0.01) {
    UseMethod("dtem")
}

#' @rdname dtem
#' @export
dtem.time_series <- function(x, events, gamma = c(-10,10), lambda = c(0,5), intercept = FALSE, incr = 0.01) {
  if (!(is.numeric(events))) {
    stop("events must be numeric")
  }
  if (TRUE %in% c(max(x$t) < events)) {
    stop("at least one event is later than maximum in x")
  }
  if (!(is.numeric(gamma) & length(gamma) == 2)) {
    stop("invalid range for gamma")
  }
  if (!(is.numeric(lambda) & length(lambda) == 2)) {
    stop("invalid range for lambda")
  }
  if (is.numeric(intercept)) {
    if (length(intercept) != 2) {
      stop("invalid range for intercept")
    }
    if (intercept[1] > intercept[2]) {
      stop("invalid range for intercept")
    }
    if ((intercept[2]-intercept[1]) < incr) {
      stop("invalid increment used")
    }
  } else {
    if (intercept != FALSE) {
      stop("invalid input for intercept")
    }
  }
  if (gamma[1] > gamma[2]) {
    stop("invalid range for gamma")
  }  
  if (lambda[1] > lambda[2]) {
    stop("invalid range for lambda")
  }
  if (lambda[1] < 0 | lambda[2] < 0) {
    stop("lambda must be non-negative")
  }
  if ((gamma[2]-gamma[1]) < incr | (lambda[2]-lambda[1]) < incr) {
    stop("invalid increment used")
  }

  gammas <- seq(gamma[1], gamma[2], by = incr)
  lambdas <- seq(lambda[1], lambda[2], by = incr)
  if (nlevels(x$variate) > 1) {
    total <- table(x$variate)
    weights <- 1/total[x$variate]
  } else {
    weights <- rep(1, nrow(x))
  }

  if (!is.numeric(intercept)) {
    out0 <- dtem_arma(as.matrix(x[,1:2]), events, gammas, lambdas, weights) 
    idx <- which(out0 == min(out0), arr.ind = TRUE)
    coefs <- c(gammas[idx[1]], lambdas[idx[2]], out0[idx])
    names(coefs) <- c("gamma", "lambda", "rss")
  } else {
    intercepts <- seq(intercept[1], intercept[2], by = incr)
    out0 <- dtem_intrcpt_arma(as.matrix(x[,1:2]), events, gammas, lambdas, intercepts, weights) 
    idx <- which(out0 == min(out0), arr.ind = TRUE)
    coefs <- c(gammas[idx[1]], lambdas[idx[2]], intercepts[idx[3]], out0[idx])
    names(coefs) <- c("gamma", "lambda", "intercept", "rss")
  }

  out <- list(coefs = coefs,  rss = out0)
  class(out) <- c("dtem", "list")
  return(out)
}


#' @export 
print.dtem <- function(x, ...) {
  cat("Decaying Trend-after-Event Regression: \n")
  print(x$coefs)
}







#' Create a Time Series Object
#'
#' Creates an object of class \code{time_series} from a matrix or data frame, or by giving a filename for a \code{csv} file in the working directory. The first column contains time indices (\code{t}), the second the value of observations (\code{x}), both of which must be numeric, and the third contains a label as a factor column (\code{variate}). Multiple time series are contained in the same \code{time_series} object, distinguished by different labels.
#' 
#' @param input A matrix, data frame, or a filename given as character.
#' @param t If used, ignores ime indices of the observations; requires also input  \code{x}
#' @param x The values of the 
#' @param variate The label of the variate if 
#' @param log Whether to take the logarithm of observations for \code{x} (default is \code{FALSE}).
#' @param inverse Whether to use the inverse of observations (default is \code{FALSE}). If both \code{log} and \code{inverse} are true, the log inverse (\eqn{\log(1/x)}) is used.
#' @param standard Whether to standardize the observations, using \eqn{z = (x -\bar{x})/s}, grouped by variate (default is \code{FALSE}). This operation is performed after \code{log} and/or \code{inverse}.
#' @param ... Additional arguments for importing from csv.
#' \itemize{
#'  \item \code{header} If a filename, whether the csv file contains a header (default \code{TRUE}).
#'  \item \code{tabular} If a filename, whether the data are presented in a tabular fonrmat, as opposed to "long" (tidyverse) format (default \code{TRUE}). 
#' } 
#' 
#' @examples
#' idx <- runif(100,0,100)
#' x <- cumsum(rnorm(100,0,1))
#' label <- rep("x", 100)
#' mat <- matrix(c(idx, x, label), ncol = 3)
#' time_series(mat)
#' 
#' dat <- data.frame(idx, x, label)
#' time_series(dat)
#' 
#' time_series(t = idx, x = x, variate = label)
#'  
#' @returns A time_series object.
#' 
#' @export
time_series <- function(input, t = NULL, x = NULL, variate = NULL, log = FALSE, inverse = FALSE, standard = FALSE, ...) {
  UseMethod("time_series")
}
#' 
#' @rdname time_series
#' @export
time_series.matrix <- function(input, t = NULL, x = NULL, variate = NULL, log = FALSE, inverse = FALSE, standard = FALSE, ...) {
  out <- time_series.data.frame(input, log = log, inverse = inverse, standard = standard)
  return(out)
}
#' 
#' @rdname time_series
#' @export
time_series.data.frame <- function(input, t = NULL, x = NULL, variate = NULL, log = FALSE, inverse = FALSE, standard = FALSE, ...) {
  out <- input
  if (NA %in% out) {
    stop('input cannot contain NA')
  }
  if (ncol(out) == 2) {
    out <- data.frame(t = as.numeric(out[,1]), x = as.numeric(out[,2]), variate = factor(rep("x", nrow(out))))
  } else if (ncol(out) == 3) {
    out <- data.frame(t = as.numeric(out[,1]), x = as.numeric(out[,2]), variate = factor(as.character(out[,3])))
  } else {
    stop('invalid dimensions. If reading from csv check whether "tabular = TRUE" should be used')
  }
  if (inverse == TRUE) {
    out$x <- 1/(out$x)
  }
  if (log == TRUE) {
    out$x <- log(out$x)
  }
  if (standard == TRUE) {
    out$x <- array2DF(tapply(out$x, out$variate, scale))$Value
    colnames(out) <- c("t", "z", "variate")
  }
  class(out) <- c("time_series", "data.frame")
  return(out)
}
#' 
#' @rdname time_series
#' @export
time_series.character <- function(input, t = NULL, x = NULL, variate = NULL, log = FALSE, inverse = FALSE, standard = FALSE, ...) {
  out <- input
  other <- list(...)
  if (!"header" %in% names(other)) {
    header <- FALSE
  } else {
    header <- other$header
  }
  if (!"tabular" %in% names(other)) {
    tabular <- FALSE
  } else {
    tabular <- other$tabular
  }
  out <- utils::read.csv(out, header = header)
  if (tabular == TRUE) {
    out <- data.frame(out[1], utils::stack(out[2:ncol(out)]))
  }
  out <- out[!is.na(out[,2]),]
  rownames(out) <- NULL
  out <- time_series.data.frame(input = out, log = log, inverse = inverse, standard = standard)
  return(out)
}
#' 
#' @rdname time_series
#' @export
time_series.numeric <- function(input = NULL, t = NULL, x = NULL, variate = NULL, log = FALSE, inverse = FALSE, standard = FALSE, ...) {
  if (is.null(t) | is.null(x)) {
    stop("no input given; t and/or x are NULL")
  }
  if (length(t) != length(x)) {
    stop("t and x are of different length")
  }
  if (is.null(variate)) {
    variate <- rep("x", length(t))
  } else {
    if (length(variate) == 1) {
      variate <- rep(variate, length(t))
    } else {
      if (length(variate) != length(t)) {
        stop("t and variate are of different length")
      }
    }
  }
  input <- data.frame(t = t, x = x, variate = variate)
  out <- time_series.data.frame(input = input, log = log, inverse = inverse, standard = standard)
  return(out)
}

# #' @export 
# print.time_series <- function(...) {
# }

