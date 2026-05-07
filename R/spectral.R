#' Least Squares Spectral Analysis (LSSA)
#'
#' Peforms a simple least squares fitting to time indexed data of the form \eqn{x(t) = \beta_{0} + \sum_{i =1}^n \beta_{1i} \cos (2 \pi f t) + \beta_{2i} \sin(2 \pi f t)}, using a range of potential frequencies \insertCite{vanicek_approximate_1969,vanicek_further_1971}{arkhaia}. Intercept may be ommitted. 
#' 
#' @param x A data frame of two columns, the first containing time indices and the second containing values.
#' @param freqs A vector of frequencies to evaluate. By default a grid from 0.001 to 0.5 is tested at an interval of 0.005.
#' @param intercept Whether to include the intercept. Default is \code{TRUE}.
#' @param type Type of output. Can be either \code{"frequency"} (the default) or \code{"period"}.
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
    out <- data.frame(freq = freqs, power = P[,1], rss = P[,2], coef_cos = P[,3], coef_sin = P[,4], coef_const = P[,5])
  } else if (type == "period") {
    out <- data.frame(period = 1/freqs, power = P[,1], rss = P[,2], coef_cos = P[,3], coef_sin = P[,4], coef_const = P[,5])
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
#' Peforms a simple least squares fitting to time indexed data of the form  \eqn{x(t) = \beta_{0} + \sum_{i =1}^n \beta_{1i} \cos (2 \pi f t) + \beta_{2i} \sin(2 \pi f t)}, using an input of frequencies; the lowest frequency peak (not the highest power frequency) is chosen for regression, up to a chosen number of iterations. Intercept may be ommitted. The lowest frequency is equivalent to the longest period.
#' 
#' @param x A data frame of two columns, with the first column containing time indices and the second containing values.
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
LSSA_LFI <- function(x, n_iter = 1, intercept = TRUE, AIC = FALSE) {
    UseMethod("LSSA_LFI")
}

#' @rdname LSSA_LFI
#' @export
LSSA_LFI.matrix <- function(x, n_iter = 1, intercept = TRUE, AIC = FALSE) {
  if (intercept == TRUE) {
    intrcpt = 1
  } else {
    intrcpt = 0
  }
  out0 <- LSSA_LFI_arma(x, n_iter, intrcpt)

  out <- list(coefs = out0[,1:3], freqs = out0[,4], rss = out0[,5], aic = out0[,6])

  if (AIC == TRUE) {
    idx <- which.min(out$aic)
    out <- list(coefs = out0[1:idx,1:3], freqs = out0[1:idx, 4], rss = out0[1:idx, 5], aic = out0[1:idx, 6])
  }
  class(out) <- c("LSSA_LFI", "list")

  return(out)
}

#' @rdname LSSA_LFI
#' @export
LSSA_LFI.data.frame <- function(x, n_iter = 1, intercept = TRUE, AIC = FALSE) {
  x <- as.matrix(x)
  out <- LSSA_LFI.matrix(x, n_iter, intercept, AIC) 
  return(out)
}




#' Linear Dependence of LSSA-LFI Candidate Models via AIC
#'
#' For a set of time series (namely partitions of set of series) contained in a list, will compute the Akaike Information Criterion (AIC) for each candidate set.  
#' 
#' @param x A \code{list} containing the time series, each of which should be a matrix or data frame with time index in the first column and value in the second.
#' @param sets Candidate sets to evaluate; must be a \code{list} of \code{lists} containing the indices of the sets in \code{x}. If left \code{NULL}, two sets are evaluated: all series pooled together [1] and all series kept separate [2].
#' @param n_iter The number of iterations to run for the least squares spectral analysis via lowest frequency iteration (LSSA-LFI). Default is 1.
#' @param intercept Whether to include the intercept in the least squares spectral analysis via lowest frequency iteration (LSSA-LFI). Default is \code{TRUE}.
#' @returns The index of the set yielding the lowest AIC. If \code{sets} is \code{NULL}, the output of \code{[1]} indicates linear dependence; if the output is \code{[2]}, linear independence.
#' 
#' @export
LSSA_LFI_candidates <- function(x, sets = NULL, n_iter = 1, intercept = TRUE) {
    UseMethod("LSSA_LFI_candidates")
}

#' @rdname LSSA_LFI_candidates
#' @export
LSSA_LFI_candidates.list <- function(x, sets = NULL, n_iter = 1, intercept = TRUE) {
  n <- length(x)
  y <- as.matrix(do.call(rbind, x))
  x_rows <- unlist(lapply(x, nrow))
  idx <- c(0, cumsum(x_rows)[1:(n-1)]) 

  # if sets is null; first set is pooled, second separate
  if ( is.null(sets) ) {
    pld <- rep(0, length(x))
    tog <- 0:(length(x)-1)
    sets0 <- matrix(c(pld, tog), nrow = 2, byrow = TRUE)
  } else {
    sets0 <- partition_list_to_matrix_arma(sets)
  }

  if (intercept == TRUE) {
    intrcpt <- 1
  } else {
    intrcpt <- 0
  }

  # add 1 for R indexing
  out <- LSSA_LFI_candidates_arma(y, idx, x_rows, sets0, n_iter, intrcpt) + 1
  
  return(out)
}



#' Trim to Epoch
#'
#' For a data frame in which the first column contains a time index and the second colum observations, trim the data frame to include observations only with a given epoch (time period).
#' 
#' @param x A data frame or matrix containing time-indexed data, with time index in the first column and value in the second.
#' @param epoch A numeric vector given the start and end time indices of the epoch.
#' @returns A data frame containing only those observations with time indices within the epoch.
#' 
#' @export
trim_epoch <- function(x, epoch = NULL) {
    UseMethod("trim_epoch")
}

#' @rdname trim_epoch
#' @export
trim_epoch.matrix <- function(x, epoch = NULL) {
  if (!(is.numeric(epoch) & length(epoch) == 2)  ) {
    stop("Epoch must a numeric vector of length 2 indicating start and end.")
  }
  if (epoch[2] < epoch[1]) {
    stop("Start of epoch must be earlier than end.")
  }

  out <- x[x[,1] > epoch[1] & x[,1] < epoch[2]  ,]
  out <- as.data.frame(out)
  return(out)
}

#' @rdname trim_epoch
#' @export
trim_epoch.data.frame <- function(x, epoch = NULL) {
  x <- as.matrix(x)
  trim_epoch.matrix(x, epoch)
}





#' Repartition (Group) Datasets
#'
#' Given a list of data frames (or matrices), labeled \code{x}, and another list which contains indices of \code{x} as vectors, labeled \code{set}, returns a new list which pools (row-binds) the data frames together  according to the grouprings in \code{set}.
#' 
#' @param x A \code{list} containing the time series, each of which should be a matrix or data frame with time index in the first column and value in the second.
#' @param set A \code{list} of vectors containing the indices of the sets in \code{x}, according to which to pool the datasets into a new list. 
#' @returns A \code{list} containing the data according to the partition structure described by \code{set}.
#' 
#' @export
repartition <- function(x, set = NULL) {
    UseMethod("repartition")
}

#' @rdname repartition
#' @export
repartition.list <- function(x, set = NULL) {
  if (!is.list(set)) {
    stop("Repartitioning set must be provided.")
  }
  out <- vector("list", length(set))
  for (i in 1:length(set)) {
    idx <- set[[i]]
    sub <- x[idx]
    y <- as.matrix(do.call(rbind, sub))
    out[[i]] <- y
  }
  return(out)
}




# #' Model Selection of LSSA-LFI Candidates
# #'
# #' For an \code{LSSA_LFI_AIC} object see \code{\link[arkhaia]{LSSA_LFI_candidates}}), returns the index of the candidate model. If \code{sets} is \code{NULL} in the \code{LSSA_LFI_candidates} function, an index of [1] refers to the model of a pooled (joint) grouping, [2] refers to the model of discrete, separate groupings. 
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
#' Generates a data frame of values \eqn{f(t)} of the model generated by LSSA-LFI (see \code{\link[arkhaia]{LSSA_LFI}})).
#' 
#' @param x A data frame where time indices are in the first column and values are in the second.
#' @param t_  A vector giving samples range of \eqn{t} for computing \eqn{f(t)}. Default is from the minimum to maximum time index sampled at 0.01 intervals.
#' @param n_iter The number of iterations to run. Default is 1.
#' @param intercept Whether to include the intercept in the least squares spectral analysis via lowest frequency iteration (LSSA-LFI). Default is \code{TRUE}.
#' @param label Default is \code{"model"}.
#' @returns A data frame containing \eqn{t, f(t)}.
#' 
#' @export
LSSA_LFI_model <- function(x, t_ = NULL, n_iter = 1, intercept = TRUE, label = "model") {
    UseMethod("LSSA_LFI_model")
}

#' @rdname LSSA_LFI_model
#' @export
LSSA_LFI_model.matrix <- function(x, t_ = NULL, n_iter = 1, intercept = TRUE, label = "model") {  
  if (intercept == TRUE) {
    intrcpt <- 1
  } else {
    intrcpt <- 0
  }

  if (is.null(t_)) {
    t_ <- seq(min(x[,1]), max(x[,1]), by = 0.01)
  }

  Y <- LSSA_LFI_model_arma(x[,1:2], t_, n_iter, intrcpt) 

  out <- data.frame(t = Y[,1], y = Y[,2], label =  rep(label, length(t_)))
  return(out)
}

#' @rdname LSSA_LFI_model
#' @export
LSSA_LFI_model.data.frame <- function(x, t_ = NULL, n_iter = 1, intercept = TRUE, label = "model") {  
  x <- as.matrix(x[,1:2])
  LSSA_LFI_model.matrix(x, t_ = t_, n_iter = n_iter, intercept = intercept, label = label)
}



#' Validated Linear Dependence via LSSA-LFI
#'
#' Probability of linear dependence between two groups of time series observations using a LSSA-LFI model selection (see \code{\link[arkhaia]{LSSA_LFI_candidates}}), on the basis of the inclusion of a third "attendant" variate. The third attendant variate is selected from the remaining time series in the list, hence the input list must include at least three times series.
#' 
#' @param x A list of data frames.
#' @param pair The pair of series to evaluate in the list \code{x}, either names or indices.
#' @param n_iter The number of iterations to run. Default is 1.
#' @param intercept Whether to include the intercept in the least squares spectral analysis via lowest frequency iteration (LSSA-LFI). Default is \code{TRUE}.

#' @returns The probability of linear dependence between two time series, in which 1 indicates linear dependence almost surely and 0 indicates independence.
#' 
#' @export
LSSA_LFI_validated <- function(x, pair = NULL, n_iter = 1, intercept = TRUE) {
    UseMethod("LSSA_LFI_validated")
}

#' @rdname LSSA_LFI_validated
#' @export
LSSA_LFI_validated.list <- function(x, pair = NULL, n_iter = 1, intercept = TRUE) {
  if (length(x) < 3) {
    stop("Must have at least three time series in x.")
  }
  if (is.null(pair)) {
    stop("Must supply the pair of data frames in x to be evaluated.")
  }
  if (typeof(pair) == "character") {
    pair <- which(names(x) %in% pair)
  }
  n <- length(x)
  nx <- 1:n
  attend <- nx[-pair]
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
    ii <- c(pair, attend[i])
    y <- x[ii]

    x_rows <- unlist(lapply(y, nrow))
    idx <- c(0, cumsum(x_rows)[1:(n-1)]) 
    y0 <- as.matrix(do.call(rbind, y))

    mod <- LSSA_LFI_candidates_arma(y0, idx, x_rows, partition, n_iter, intrcpt) + 1 # add 1 for R indexing
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
#' Evaluate pairwise linear dependence between observations using a LSSA-LFI valdiated model selection (see \code{\link[arkhaia]{LSSA_LFI_validated}}). 
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
      pair <- c(i,j)
      mat[i,j] <- LSSA_LFI_validated(x, pair, n_iter = n_iter, intercept = intercept)
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
  if (intercept == TRUE) {
    intrcpt <- 1
  } else {
    intrcpt <- 0
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
      if (t_[ti] + h_[hi] <= end) {
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

          out[(length(h_) - hi)+1,ti] <- LSSA_LFI_validated(tmp, pair, n_iter_, intrcpt)
        }
      }
    }
    message("Percent complete: ", round(ti/length(t_), 2)*100, "%")
  }


  rownames(out) <- t_
  colnames(out) <- h_
  return(out)

}



