#' Cressie-Read Power-Divergence Statistic
#'
#' For a matrix of cross-tabulated counts of observations, computes the Cressie-Read power-divergence statistic according to the selection of a parameter \eqn{\lambda} \insertCite{cressie_multinomial_1984,read_goodness--fit_1988}{arkhaia}.
#' 
#' @param x A matrix of cross-tabulated counts.
#' @param lambda The model parameter of the Cressie-Read power-divergence statistic. Default is the recommended value of 2/3.
#' @examples 
#' x1 <- c(2, 0, 10, 11, 5)
#' x2 <- c(1, 1, 17, 23, 3)
#' x3 <- c(0, 0, 2, 81, 11)
#' 
#' x <- matrix(c(x1, x2, x3), ncol = 3)
#' 
#' CR(x)
#' CR(x, lambda = 1)
#'
#' @returns The Cressie-Read power-divergence statistic.
#' 
#' @references
#'   \insertAllCited{}
#' 
#' @importFrom Rdpack reprompt
#' @export
CR <- function(x, lambda = 2/3) {
    UseMethod("CR")
}

#' @rdname CR
#' @export
CR.matrix <- function(x, lambda = 2/3) {
  x_e <- outer(rowSums(x), colSums(x)) / sum(x)
  out <- (2/(lambda * (lambda + 1))) * sum(x * ((x / x_e)^(lambda) - 1))
  return(out)
}

#' @rdname CR
#' @export
CR.data.frame <- function(x, lambda = 2/3) {
  x <- as.matrix(x)
  out <- CR.matrix(x, lambda = lambda)
  return(out)
}

#' @rdname CR
#' @export
CR.xtabs <- function(x, lambda = 2/3) {
  x <- as.matrix(x)
  out <- CR.matrix(x, lambda = lambda)
  return(out)
}

#' @rdname CR
#' @export
CR.table <- function(x, lambda = 2/3) {
  x <- as.matrix(x)
  out <- CR.matrix(x, lambda = lambda)
  return(out)
}




#' Bias-Corrected Cramer's V
#'
#' For a matrix of cross-tabulated counts of observations, estimates Cramer's V using Bergsma's bias correction \insertCite{bergsma_bias-correction_2013}{arkhaia}.
#' 
#' @param x A matrix or data frame of cross-tabulated counts.
#' @param method Method for estimating chi squared. Default is \code{CR}, the Cressie-Read method (see \code{\link[arkhaia]{CR}})). Other options include \code{Pearson}, which uses R's base \code{\link[stats]{chisq.test}}) function, and \code{G}, which computes the likelihood-ratio \eqn{G} statistic, which is the case for the Cressie-Read statistic as \eqn{\lambda} approaches 0.
#' @param lambda Only applied to the method \code{CR}. Default is the recommended value of 2/3.
#' @param simulate Only applied ot the method \code{Pearson}.  Whether to simulate \eqn{p} values in the \code{\link[stats]{chisq.test}}) function. Default is \code{TRUE}.
#' @examples 
#' x1 <- c(2, 0, 10, 11, 5)
#' x2 <- c(1, 1, 17, 23, 3)
#' x3 <- c(0, 0, 2, 81, 11)
#' 
#' x <- matrix(c(x1, x2, x3), ncol = 3)
#' 
#' VB_pair(x)
#' VB(x, lambda = 1)
#' 
#' @returns Bergsma's bias-corrected etimate of Cramer's \eqn{V}.
#' 
#' @references
#'  \insertAllCited{}
#' 
#' @importFrom Rdpack reprompt
#' @export
VB <- function(x, method = "CR", lambda = 2/3, simulate = TRUE) {
    UseMethod("VB")
}

#' @rdname VB
#' @export
VB.matrix <- function(x, method = "CR", lambda = 2/3, simulate = TRUE) {
  if (method == "CR") {
    chi2 <- CR(x, lambda)
  } else if (method == "Pearson") {
    chi2 <- stats::chisq.test(x, simulate.p.value = simulate)$statistic
  } else if (method == "G") {
    x_e <- outer(rowSums(x), colSums(x)) / sum(x)
    chi2 <- 2 * sum(x * log(x / x_e))
  } else {
    stop('method must be CR, Pearson, or G')
  }
  phi2 <- chi2 / sum(x)
  phi2_plus <- max(c(phi2 - (nrow(x) - 1) * (ncol(x) - 1) / (sum(x) - 1), 0))
  out <- sqrt(phi2_plus/(min(c(nrow(x), ncol(x))) - 1  )) 
  return(out)
}

#' @rdname VB
#' @export
VB.data.frame <- function(x, method = "CR", lambda = 2/3, simulate = TRUE) {
  x <- as.matrix(x)
  out <- VB.matrix(x, method = method, lambda = lambda, simulate = simulate)
  return(out)
}

#' @rdname VB
#' @export
VB.xtabs <- function(x, method = "CR", lambda = 2/3, simulate = TRUE) {
  x <- as.matrix(x)
  out <- VB.matrix(x, method = method, lambda = lambda, simulate = simulate)
  return(out)
}

#' @rdname VB
#' @export
VB.table <- function(x, method = "CR", lambda = 2/3, simulate = TRUE) {
  x <- as.matrix(x)
  out <- VB.matrix(x, method = method, lambda = lambda, simulate = simulate)
  return(out)
}


#' Bias-Corrected Cramer's V Pairwise between Columns
#'
#' For a matrix or data frame of cross-tabulated counts of observations, estimates Cramer's V using Bergsma's bias correction \insertCite{bergsma_bias-correction_2013}{arkhaia} by subsetting the table by pairs of columns. See \code{\link[arkhaia]{VB}}
#' 
#' @param x A matrix or data frame of cross-tabulated counts.
#' @param method Method for estimating chi squared. Default is \code{CR}, the Cressie-Read method (see \code{\link[arkhaia]{CR}})). Other options include \code{Pearson}, which uses R's base \code{\link[stats]{chisq.test}}) function, and \code{G}, which computes the likelihood-ratio \eqn{G} statistic, which is the case for the Cressie-Read statistic as \eqn{\lambda} approaches 0.
#' @param lambda Only applied to the method \code{CR}. Default is the recommended value of 2/3.
#' @param simulate Only applied ot the method \code{Pearson}.  Whether to simulate \eqn{p} values in the \code{\link[stats]{chisq.test}}) function. Default is \code{TRUE}.
#' @examples 
#' x1 <- c(2, 0, 10, 11, 5)
#' x2 <- c(1, 1, 17, 23, 3)
#' x3 <- c(0, 0, 2, 81, 11)
#' 
#' x <- matrix(c(x1, x2, x3), ncol = 3)
#' 
#' VB_pair(x)
#' VB_pair(x, lambda = 1)
#' 
#' @returns A matrix of Bergsma's bias-corrected etimate of Cramer's \eqn{V}, pairwise between columns of the input matrix, of class \code{effect_sizes}.
#' 
#' @references
#'  \insertAllCited{}
#' 
#' @importFrom Rdpack reprompt
#' @export
VB_pair <- function(x, method = "CR", lambda = 2/3, simulate = TRUE) {
    UseMethod("VB_pair")
}

#' @rdname VB_pair
#' @export
VB_pair.matrix <- function(x, method = "CR", lambda = 2/3, simulate = TRUE) {

  M_V <- matrix(NA, ncol(x), ncol(x))

  for (i in 1:(ncol(x))) {
    for (j in 1:ncol(x)) {
      
      Y <- x[,c(i,j)]
      Y <- t(Y[rowSums(Y) > 0  ,  ])

      V <- VB(Y, method = method, lambda = lambda, simulate = simulate)
      M_V[i,j] <- V
    }
  }

  names_ <- colnames(x)
  rownames(M_V) <- names_
  colnames(M_V) <- names_

  class(M_V) <- c("effect_sizes", "matrix")

  return(M_V)
}

#' @rdname VB_pair
#' @export
VB_pair.data.frame <- function(x, method = "CR", lambda = 2/3, simulate = TRUE) {
  x <- as.matrix(x)
  out <- VB_pair(x, method = method, lambda = lambda, simulate = simulate)
  return(out)
}

#' @rdname VB_pair
#' @export
VB_pair.xtabs <- function(x, method = "CR", lambda = 2/3, simulate = TRUE) {
  x <- as.matrix(x)
  out <- VB_pair.matrix(x, method = method, lambda = lambda, simulate = simulate)
  return(out)
}

#' @rdname VB_pair
#' @export
VB_pair.table <- function(x, method = "CR", lambda = 2/3, simulate = TRUE) {
  x <- as.matrix(x)
  out <- VB_pair.matrix(x, method = method, lambda = lambda, simulate = simulate)
  return(out)
}




#' Leave-One-Out Type Routine for Pairwise Effect Sizes of Archaeological Contexts
#'
#' To evaluate the stability of estimates of effect sizes between archaeological contexts in light of the inclusion or exclusion of a given type, this routine computes bias-corrected Cramer's V \code{\link[arkhaia]{VB}}) omitting a type on each iteration.
#' 
#' @param x A contingency table as a matrix or data frame expressing counts, with contexts along columns and types along rows.
#' @param method Method for estimating chi squared. Default is \code{CR}, the Cressie-Read method (see \code{\link[arkhaia]{CR}})). Other options include \code{Pearson}, which uses R's base \code{\link[stats]{chisq.test}}) function, and \code{G}, which computes the likelihood-ratio \eqn{G} statistic, which is the case for the Cressie-Read statistic as \eqn{\lambda} approaches 0.
#' @param lambda Only applied to the method \code{CR}. Default is the recommended value of 2/3.
#' @param simulate Only applied ot the method \code{Pearson}.  Whether to simulate \eqn{p} values in the \code{\link[stats]{chisq.test}}) function. Default is \code{TRUE}.
#' @examples 
#' x1 <- c(2, 0, 10, 11, 5)
#' x2 <- c(1, 1, 17, 23, 3)
#' x3 <- c(0, 0, 2, 81, 11)
#' 
#' x <- matrix(c(x1, x2, x3), ncol = 3)
#' 
#' rownames(x) <- LETTERS[1:nrow(x)]
#' colnames(x) <- c("S1", "S2", "S3")
#' 
#' VB_LOO_type(x)
#' 
#' @returns A three-dimensional array of the pairwise context-by-context effect sizes given the ommission of a type, for all types given in the input data frame.
#' 
#' @export
VB_LOO_type <- function(x, method = "CR", lambda = 2/3, simulate = TRUE) {
    UseMethod("VB_LOO_type")
}

#' @rdname VB_LOO_type
#' @export
VB_LOO_type.matrix <- function(x, method = "CR", lambda = 2/3, simulate = TRUE) {
  if (is.null(colnames(x)) | is.null(rownames(x))) {
    stop("Rows and and columns must be named.")
  }
  contexts <- colnames(x)
  types <- rownames(x)

  V_mat <- array(NA, dim = c(length(contexts), length(contexts), length(types)))
  for (k in 1:length(types)) {
    Y <- x[-k, ]
    Y <- Y[rowSums(Y) > 0, ]
    Y <- Y[, colSums(Y) > 0]

    if (nrow(Y) > 1 & ncol(Y) > 1) {
    for (i in 1:(ncol(Y))) {
      for (j in 1:ncol(Y)) {
        
        Z <- Y[,c(i,j)]
        Z <- t(Z[rowSums(Z) > 0  ,  ])

        if (nrow(Z) > 1 & ncol(Z) > 1) {
          eff <- VB(Z)#, method = method, lambda = lambda, simulate = simulate)
        } else {
          eff <- NA
        }

        V_mat[i,j,k] <- eff
      }
    }
    }
  }

  rownames(V_mat) <- colnames(Y)
  colnames(V_mat) <- colnames(Y)
  class(V_mat) <- c("effect_sizes3d", "array")

  return(V_mat)
}

#' @rdname VB_LOO_type
#' @export
VB_LOO_type.data.frame <- function(x, method = "CR", lambda = 2/3, simulate = TRUE) {
  x <- as.matrix(x)
  out <- VB_LOO_type.matrix(x, method = method, lambda = lambda, simulate = simulate)
  return(out)
}

#' @rdname VB_LOO_type
#' @export
VB_LOO_type.xtabs <- function(x, method = "CR", lambda = 2/3, simulate = TRUE) {
  x <- as.matrix(x)
  out <- VB_LOO_type.matrix(x, method = method, lambda = lambda, simulate = simulate)
  return(out)
}

#' @rdname VB_LOO_type
#' @export
VB_LOO_type.table <- function(x, method = "CR", lambda = 2/3, simulate = TRUE) {
  x <- as.matrix(x)
  out <- VB_LOO_type.matrix(x, method = method, lambda = lambda, simulate = simulate)
  return(out)
}





#' Log Odds Ratio with Haldane-Anscombe's Correction Pairwise between Columns
#'
#' For a contingency table with columns representing contexts and rows representing types, computes the log odds ratio of a 2 x 2 contingency table of the presence-absences across each pair of columns. The addition of 0.5 to cells of the contingency table is performed \insertCite{anscombe_estimating_1956,haldane_estimation_1956}{arkhaia}.
#' 
#' @param x A matrix or data frame with contexts along columns and types along rows.
#' @examples 
#' x1 <- c(2, 0, 0, 11, 5, 0, 2, 0, 4)
#' x2 <- c(1, 1, 0, 23, 3, 3, 0, 0, 0)
#' x3 <- c(1, 0, 0, 0, 10, 0, 4, 0, 1)
#' 
#' x <- data.frame(S1 = x1, S2 = x2, S3 = x3)
#' rownames(x) <- LETTERS[1:nrow(x)]
#' 
#' log_OR_pair(x)
#' 
#' @returns A matrix giving the log odds ratio between all columns of the input.
#' 
#' @export
log_OR_pair <- function(x) {
    UseMethod("log_OR_pair")
}

#' @rdname log_OR_pair
#' @export
log_OR_pair.matrix <- function(x) {

  X <- matrix(NA, ncol(x), ncol(x))

  for (i in 1:(ncol(x))) {
    for (j in 1:ncol(x)) {
      Y <- x[,c(i,j)]
      PA <- pa_matrix(Y)
      PA <- PA + 0.5
      X[i,j] <- log(PA[1,1] * PA[2,2] / (PA[1,2] * PA[2,1]))
    }
  }

  rownames(X) <- colnames(x)
  colnames(X) <- colnames(x)

  class(X) <- c("effect_sizes", "matrix")

  return(X)
}

#' @rdname log_OR_pair
#' @export
log_OR_pair.data.frame <- function(x) {
  x <- as.matrix(x)
  out <- log_OR_pair.matrix(x)
  return(out)
}

#' @rdname log_OR_pair
#' @export
log_OR_pair.xtabs <- function(x) {
  x <- as.matrix(x)
  out <- log_OR_pair.matrix(x)
  return(out)
}

#' @rdname log_OR_pair
#' @export
log_OR_pair.table <- function(x) {
  x <- as.matrix(x)
  out <- log_OR_pair.matrix(x)
  return(out)
}


#' Trace Coefficient Pairwise between Columns
#'
#' For a contingency table with columns representing contexts and rows representing types, computes the trace of a 2 x 2 contingency table of the presence-absences across each pair of columns.
#' 
#' @param x A matrix or data frame with contexts along columns and types along rows.
#' @examples 
#' x1 <- c(2, 0, 0, 11, 5, 0, 2, 0, 4)
#' x2 <- c(1, 1, 0, 23, 3, 3, 0, 0, 0)
#' x3 <- c(1, 0, 0, 0, 10, 0, 4, 0, 1)
#' 
#' x <- data.frame(S1 = x1, S2 = x2, S3 = x3)
#' rownames(x) <- LETTERS[1:nrow(x)]
#' 
#' trace_pair(x)
#' 
#' @returns A matrix giving the trace between all columns of the input.
#' 
#' @export
trace_pair <- function(x) {
    UseMethod("trace_pair")
}

#' @rdname trace_pair
#' @export
trace_pair.matrix <- function(x) {

  X <- matrix(NA, ncol(x), ncol(x))

  for (i in 1:(ncol(x))) {
    for (j in 1:ncol(x)) {
      Y <- x[,c(i,j)]
      PA <- pa_matrix(Y)
      X[i,j] <- PA[1,1] + PA[2,2]
    }
  }

  rownames(X) <- colnames(x)
  colnames(X) <- colnames(x)

  class(X) <- c("effect_sizes", "matrix")

  return(X)
}

#' @rdname trace_pair
#' @export
trace_pair.data.frame <- function(x) {
  x <- as.matrix(x)
  out <- trace_pair.matrix(x)
  return(out)
}

#' @rdname trace_pair
#' @export
trace_pair.xtabs <- function(x) {
  x <- as.matrix(x)
  out <- trace_pair.matrix(x)
  return(out)
}

#' @rdname trace_pair
#' @export
trace_pair.table <- function(x) {
  x <- as.matrix(x)
  out <- trace_pair.matrix(x)
  return(out)
}




#' Presence-Absence Matrix
#'
#' Create a 2 x 2 contingency table of the presence/absence of a given type 
#' 
#' @param x A matrix or data frame where the two columns indicate contexts and the rows indicate types, with each cell containing counts of types in a context. 
#' @examples 
#' x1 <- c(2, 0, 10, 11, 5)
#' x2 <- c(1, 1, 17, 23, 3)
#' 
#' x <- data.frame(S1 = x1, S2 = x2)
#' rownames(x) <- LETTERS[1:5]
#' 
#' pa_matrix(x)
#' 
#' @returns A 2 x 2 contingency table of the counts of types present in both, either, or neither context.
#' 
#' @export
pa_matrix <- function(x) {
    UseMethod("pa_matrix")
}

#' @rdname pa_matrix
#' @export
pa_matrix.matrix <- function(x) {
  pa <- matrix(0, nrow = 2, ncol = 2)
  for (i in 1:nrow(x)) {
    ii <- 2
    if (x[i,1] > 0 ) {ii <- 1}
    jj <- 2
    if (x[i,2] > 0) {jj <- 1}
    pa[ii,jj] <- pa[ii,jj] + 1
  }
  return(pa)
}

#' @rdname pa_matrix
#' @export
pa_matrix.data.frame <- function(x) {
  x <- as.matrix(x)
  out <- pa_matrix.matrix(x)
  return(out)
}

#' @rdname pa_matrix
#' @export
pa_matrix.table <- function(x) {
  x <- as.matrix(x)
  out <- pa_matrix.matrix(x)
  return(out)
}

#' @rdname pa_matrix
#' @export
pa_matrix.xtabs <- function(x) {
  x <- as.matrix(x)
  out <- pa_matrix.matrix(x)
  return(out)
}



#' Homogeneity of Related Assemblages via Effect Size
#'
#' Given a contingency table of counts, with contexts along the columns and types along rows, as well as identified pairs of contexts ass 
#' 
#' @param x An effect_sizes object as returned by \code{\link[arkhaia]{VB_pair}}.
#' @param related The related pairs of contexts as a two-column matrix or data frame (contexts between which one anticipates a meaningful relationship). Names must match \code{colnames(x)}.
#' @param unrelated The unrelated pairs of contexts as a two-column matrix or data frame (contexts between which one does not anticipate a meaningful relationship). Names must match \code{colnames(x)}. May be left \code{NULL}, in which event all pairs not expressed in \code{related} are created
#' @param direction Whether the related or unrelated effect size should come first. Default is \code{"UW"}; alternative is \code{"WU"}.
#' 
#' @examples 
#' x1 <- c(2, 0, 10, 11, 5)
#' x2 <- c(1, 1, 17, 23, 3)
#' x3 <- c(0, 0, 2, 81, 11)
#' x4 <- c(3, 18, 9, 0, 23)
#' x <- matrix(c(x1, x2, x3, x4), ncol = 4)
#' colnames(x) <- c("surface1", "subsurface1", "surface2", "subsurface2")
#' rownames(x) <- LETTERS[1:5]
#' 
#' # related pairs
#' W_contexts <- matrix(c("surface1", "surface2", "subsurface1", "subsurface2"), ncol = 2)
#' 
#' # unrelated pairs (will be automatically created if left NULL)
#' U_contexts <- matrix(c("surface1", "surface1", "surface2", "subsurface1",
#'                        "surface2","subsurface2", "surface1", "subsurface2"), ncol = 2)
#' 
#' effect_sizes <- VB_pair(x)
#' homogeneity(effect_sizes, related = W_contexts, unrelated = U_contexts) 
#' homogeneity(effect_sizes, related = W_contexts) 
#' 
#' 
#' @returns A list of results:
#'  * \code{U} : The effect sizes between unrelated pairs of contexts.
#'  * \code{W} : The effect sizes between related pairs of contexts.
#'  * \code{Q} : The quantile indicating the proportion of related contexts more homogenous than unrelated contexts (if direction is \code{"UW"}); less homogenous if direction is set to \code{"UW"}.
#'  * \code{D} : The distribution of differences, \eqn{D_{ij} = U_j - W_i}, if direction is set to \code{"UW"}. The proportion of \eqn{D > 0} is equivalent to the mean of \code{Q}.
#' 
#' @export
homogeneity <- function(x, related = NULL, unrelated = NULL, direction = "UW") {
    UseMethod("homogeneity")
}
#' 
#' @rdname homogeneity
#' @export
homogeneity.effect_sizes <- function(x, related = NULL, unrelated = NULL, direction = "UW") {
  if (is.null(related)) {
    stop("Names of related contexts must be specified.")
  }
  if (ncol(related) > 2) {
    message("Using only first two columns of related.")
  }
  if (is.null(unrelated)) {
    relpairs <- paste0(related[,1], related[,2])
    unrelated <- matrix(NA, nrow = 0, ncol = 2)
    all <- unique(c(related[,1], related[,2]))
    for (i in 1:(length(all)-1)) {
      for (j in (i+1):length(all)) {
        if (!(paste0(all[i],all[j]) %in% relpairs)) {
          unrelated <- rbind(unrelated, c(all[i], all[j]))
        }
      }
    }
  }

  check_colnames1 <- match(related, colnames(x))
  if (TRUE %in% is.na(check_colnames1)) {
    stop("Names of related contexts not found in colnames of input x.")
  }
  check_colnames2 <- match(unrelated, colnames(x))
  if (TRUE %in% is.na(check_colnames2)) {
    stop("Names of unrelated contexts not found in colnames of input x.")
  }

  #y <- VB_pair(x, method = method, lambda = lambda, simulate = simulate)

  U <- numeric(nrow(unrelated))
  for (i in 1:nrow(unrelated)) {
    U[i] <- x[unrelated[i, 1], unrelated[i, 2]]
  }

  W <- numeric(nrow(related))
  for (i in 1:nrow(related)) {
    W[i] <- x[related[i,1], related[i,2]]
  }

  U_names <- numeric(nrow(unrelated))
  W_names <- numeric(nrow(related))
  
  for (i in 1:nrow(unrelated)) {
    U_names[i] <- paste0(unrelated[i,1],'/', unrelated[i,2]) 
  }

  for (i in 1:nrow(related)) {
    W_names[i] <- paste0(related[i,1],'/', related[i,2]) 
  }
  names(U) <- U_names
  names(W) <- W_names

  if (direction == "UW") {
    Q <- c()
    for (i in 1:length(W)) {
      Q[i] <- mean(U > W[i])
    }
    names(Q) <- names(W)
    delta <- numeric(length(W) * length(U))
    k <- 1
    for (i in 1:length(U)) {
      for (j in 1:length(W)) {
        delta[k] <- U[i] - W[j]
        k <- k + 1
      }
    }
  } else if (direction == "WU") {
    Q <- c()
    for (i in 1:length(W)) {
      Q[i] <- mean(U < W[i])
    }
    names(Q) <- names(W)
    delta <- numeric(length(W) * length(U))
    k <- 1
    for (i in 1:length(U)) {
      for (j in 1:length(W)) {
        delta[k] <- W[j] - U[i]
        k <- k + 1
      }
    }

  }

  res <- list(U = U, W = W, Q = Q, D = delta)
  class(res) <- c("homogeneity", "list")
  return(res)
}




#' Leave-One-Out Routine for Homogeneity of Related Assemblages
#'
#' To evaluate the stability of estimates of effect sizes between archaeological contexts in light of the inclusion or exclusion of types and contexts, this routine computes the homogeneity of related assemblages iteratively in two ways, first by leaving out a type on each iteration and second by leaving out a context on each iteration (see \code{\link[arkhaia]{homogeneity}} and \code{\link[arkhaia]{VB}}). 
#' 
#' @param x A data frame or matrix representing a contingency table of counts, with contexts along the columns and types along rows.
#' @param related The related pairs of contexts as a two-column matrix or data frame (contexts between which one anticipates a meaningful relationship). Names must match \code{colnames(x)}.
#' @param unrelated The unrelated pairs of contexts as a two-column matrix or data frame (contexts between which one does not anticipate a meaningful relationship). Names must match \code{colnames(x)}.
#' @param direction Whether the related or unrelated effect size should come first. Default is \code{"UW"}; alternative is \code{"WU"}.
#' @param method Method for estimating chi squared. Default is \code{CR}, the Cressie-Read method (see \code{\link[arkhaia]{CR}})). Other options include \code{Pearson}, which uses R's base \code{\link[stats]{chisq.test}}) function, and \code{G}, which computes the likelihood-ratio \eqn{G} statistic, which is the case for the Cressie-Read statistic as \eqn{\lambda} approaches 0.
#' @param lambda Only applied to the method \code{CR}. Default is the recommended value of 2/3.
#' @param simulate Only applied ot the method \code{Pearson}.  Whether to simulate \eqn{p} values in the \code{\link[stats]{chisq.test}}) function. Default is \code{TRUE}.
#' 
#' @examples 
#' x1 <- c(2, 0, 10, 11, 5)
#' x2 <- c(1, 1, 17, 23, 3)
#' x3 <- c(0, 0, 2, 81, 11)
#' x4 <- c(3, 18, 9, 0, 23)
#' x <- matrix(c(x1, x2, x3, x4), ncol = 4)
#' colnames(x) <- c("surface1", "subsurface1", "surface2", "subsurface2")
#' rownames(x) <- LETTERS[1:5]
#' 
#' # related pairs
#' W_contexts <- matrix(c("surface1", "surface2", "subsurface1", "subsurface2"), ncol = 2)
#' 
#' # unrelated pairs (will be automatically created if left NULL)
#' U_contexts <- matrix(c("surface1", "surface1", "surface2", "subsurface1",
#'                        "surface2","subsurface2", "surface1", "subsurface2"), ncol = 2)
#' 
#' homogeneity_LOO(x, related = W_contexts, unrelated = U_contexts) 
#' 
#' @returns A list containing:
#'   * \code{EQ} : The mean quantile expressing the degree of homogeneity among related contexts as computed by the \code{\link[arkhaia]{homogeneity}} function.
#'   * \code{EQ_T_mean} : The mean quantile upon iterating over the omission of each type (of the LOO samples). 
#'   * \code{EQ_T_var} : The variance of the LOO type samples.
#'   * \code{EQ_C_mean} : The mean quantile upon iterating over the omission of each context.
#'   * \code{EQ_C_var} : The variance of the LOO context samples.
#' 
#' @export
homogeneity_LOO <- function(x, related = NULL, unrelated = NULL, direction = "UW", method = "CR", lambda = 2/3, simulate = TRUE) {
    UseMethod("homogeneity_LOO")
}

#' @rdname homogeneity_LOO
#' @export
homogeneity_LOO.matrix <- function(x, related = NULL, unrelated = NULL, direction = "UW", method = "CR", lambda = 2/3, simulate = TRUE) {
  # if (!("Count" %in% colnames(x) & ("Context" %in% colnames(x) & "Type" %in% colnames(x))) ) {
  #   x <- data.frame(Context = x[,1],  Type = x[,2], Count = x[,3])
  #   message("Context identified as column 1, Type as column 2, and Count as column 3.")
  # } 
  # Y <- stats::xtabs(Count ~ Type + Context, x) # Context along columns  if 

  if (is.null(colnames(x)) | is.null(rownames(x))) {
    stop("Rows and and columns must be named.")
  }

  context_names <- colnames(x)
  type_names <- rownames(x)

  V_B <- VB_pair(x, method = method, lambda = lambda, simulate = simulate)
  eff <- homogeneity(V_B, related = related, unrelated = unrelated, direction = direction)
  EQ <- mean(eff$D > 0)

  # LOO - type
  V_B_jt_3d <- VB_LOO_type(x, method = method, lambda = lambda, simulate = simulate)
  EQ_jt0 <- numeric(length(type_names))
  EQ_jt0[] <- NA
  for (k in 1:(dim(V_B_jt_3d)[3])) {
    Z <- V_B_jt_3d[,,k]
    class(Z) <- "effect_sizes"
    eff_jt00 <- homogeneity(Z, related, unrelated)
    EQ_jt0[k] <- mean(eff_jt00$D > 0)
  }

  eff_jt <- mean(EQ_jt0)
  eff_jt_var <- stats::var(EQ_jt0) # sum((EQ_jt0  - eff_jt)^2) * ((length(EQ_jt0) - 1)/length(EQ_jt0))

  EQ_jc0 <- numeric(length(context_names))
  EQ_jc0[] <- NA

  # LOO - contexts
  for (i in 1:length(context_names)) {
    related_tmp <- matrix( related[rowSums(related == context_names[i]) == 0, ]   , ncol = 2)
    unrelated_tmp <- matrix( unrelated[rowSums(unrelated == context_names[i]) == 0, ] , ncol = 2)
    if (nrow(related_tmp) > 0 & nrow(unrelated_tmp) > 0) {
      eff_tmp <- homogeneity(V_B, related = related_tmp, unrelated = unrelated_tmp, direction = direction)
    }
    EQ_jc0[i] <- mean(eff_tmp$D >0)
  }

  eff_jc <- mean(EQ_jc0)
  eff_jc_var <-  stats::var(EQ_jc0) # sum((EQ_jc0  - eff_jc)^2)  * ((length(EQ_jc0) - 1)/length(EQ_jc0))

  out <- list(EQ = EQ, EQ_T_mean = eff_jt, EQ_T_var = eff_jt_var, EQ_C_mean = eff_jc, EQ_C_var = eff_jc_var)

  return(out)
}

#' @rdname homogeneity_LOO
#' @export
homogeneity_LOO.data.frame <- function(x, related = NULL, unrelated = NULL, direction = "UW", method = "CR", lambda = 2/3, simulate = TRUE) {
  x <- as.matrix(x)
  out <- homogeneity_LOO.matrix(x, related = related, unrelated = unrelated, direction = direction, method = method, lambda = lambda, simulate = TRUE)
  return(out)
}

#' @rdname homogeneity_LOO
#' @export
homogeneity_LOO.table <- function(x, related = NULL, unrelated = NULL, direction = "UW", method = "CR", lambda = 2/3, simulate = TRUE) {
  x <- as.matrix(x)
  out <- homogeneity_LOO.matrix(x, related = related, unrelated = unrelated, direction = direction, method = method, lambda = lambda, simulate = TRUE)
  return(out)
}

#' @rdname homogeneity_LOO
#' @export
homogeneity_LOO.xtabs <- function(x, related = NULL, unrelated = NULL, direction = "UW", method = "CR", lambda = 2/3, simulate = TRUE) {
  x <- as.matrix(x)
  out <- homogeneity_LOO.matrix(x, related = related, unrelated = unrelated, direction = direction, method = method, lambda = lambda, simulate = TRUE)
  return(out)
}

