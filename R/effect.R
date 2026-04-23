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




#' Bias-Corrected Cramer's V
#'
#' For a matrix of cross-tabulated counts of observations, estimates Cramer's V using Bergsma's bias correction \insertCite{bergsma_bias-correction_2013}{arkhaia}.
#' 
#' @param x A matrix of cross-tabulated counts.
#' @param method Method for estimating chi squared. Default is \code{CR}, the Cressie-Read method (see \code{\link[arkhaia]{CR}})). Other options include \code{Pearson}, which uses R's base \code{\link[stats]{chisq.test}}) function, and \code{G}, which computes the likelihood-ratio \eqn{G} statistic, which is the case for the Cressie-Read statistic as \eqn{\lambda} approaches 0.
#' @param lambda Only applied to the method \code{CR}. Default is the recommended value of 2/3.
#' @param simulate Only applied ot the method \code{Pearson}.  Whether to simulate \eqn{p} values in the \code{\link[stats]{chisq.test}}) function. Default is \code{TRUE}.
#' #' @examples 
#' x1 <- c(2, 0, 10, 11, 5)
#' x2 <- c(1, 1, 17, 23, 3)
#' x3 <- c(0, 0, 2, 81, 11)
#' 
#' x <- matrix(c(x1, x2, x3), ncol = 3)
#' 
#' VB(x)
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






#' Leave-One-Out Routine for Pairwise Effect Sizes of Archaeological Contexts
#'
#' To evaluate the stability of estimates of effect sizes between archaeological contexts in light of the inclusion or exclusion of a given type, this routine computes bias-corrected Cramer's V \code{\link[arkhaia]{VB}}) omitting a type on each iteration.
#' 
#' @param x A data frame in which the first column is labeled \code{Context} and the second column is labeled \code{Type}, with each row containing a type present in that context.
#' @param lambda The Cressie-Read power-divergence statistic parameter. Default is the recommended value of 2/3.
#' @returns A three-dimensional array of the pairwise context-by-context effect sizes given the ommission of a type, for all types given in the input data frame.
#' 
#' @export
VB_LOO_type <- function(x, lambda = 2/3) {
    UseMethod("VB_LOO_type")
}

#' @rdname VB_LOO_type
#' @export
VB_LOO_type.data.frame <- function(x, lambda = 2/3) {
  if (!("Context" %in% colnames(x) & "Type" %in% colnames(x)) ) {
    colnames(x) <- c("Context", "Type")
    warning("Context identified as first column and type identified as second.")
  } else {

    contexts <- unique(x$Context)
    types <- unique(x$Type)

    V_mat <- array(NA, dim = c(length(contexts), length(contexts), length(types)))
    for (k in 1:length(types)) {
      dat1 <- x[x$Type != types[k], ]
      Y <- stats::xtabs( Count ~ Type + Context, dat1) # Context along columns
      Y <- Y[rowSums(Y) > 0, ]
      Y <- Y[, colSums(Y) > 0]

      if (nrow(Y) > 1 & ncol(Y) > 1) {
      for (i in 1:(ncol(Y))) {
        for (j in 1:ncol(Y)) {
          
          Z <- Y[,c(i,j)]
          Z <- t(Z[rowSums(Z) > 0  ,  ])

          if (nrow(Z) > 1 & ncol(Z) > 1) {
            eff <- VB(Z, lambda)
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
  return(V_mat)
  }
}









# effect_sizes <- function(mat, related, unrelated) {

#   U <- numeric(nrow(unrelated))

#   for (i in 1:nrow(unrelated)) {
#     U[i] <- mat[unrelated[i, 1], unrelated[i, 2]]
#   }

#   W <- numeric(nrow(related))

#   for (i in 1:nrow(related)) {
#     W[i] <- mat[related[i,1], related[i,2]]
#   }

#   U_names <- numeric(nrow(unrelated))
#   W_names <- numeric(nrow(related))
  
#   for (i in 1:nrow(unrelated)) {
#     U_names[i] <- paste0(unrelated[i,1],'/', unrelated[i,2]) 
#   }

#   for (i in 1:nrow(related)) {
#     W_names[i] <- paste0(related[i,1],'/', related[i,2]) 
#   }
#   names(U) <- U_names
#   names(W) <- W_names

#   Q <- c()

#   for (i in 1:length(W)) {
#     Q[i] <- mean(U > W[i])
#   }

#   names(Q) <- names(W)


#   delta <- numeric(length(W) * length(U))
#   k <- 1

#   for (i in 1:length(U)) {
#     for (j in 1:length(W)) {
#       delta[k] <- U[i] - W[j]
#       k <- k + 1
#     }
#   }

#   res <- list(U = U, W = W, Q = Q, D = delta)

#   return(res)
# }





# effect_sizes0 <- function(dat, related, unrelated, lambda = 2/3) {
#   if (!("Count" %in% colnames(dat) & ("Context" %in% colnames(dat) & "Type" %in% colnames(dat))) ) {
#     message('colnames in data frame must contain labels "Count", "Context", and "Type"')
#   } else {
#     Y <- xtabs(Count ~ Type + Context, dat) # Context along columns
#     context_names <- colnames(Y)
#     type_names <- rownames(Y)

#     V_B <- V_B_mat(Y, lambda = lambda)
#     eff <- effect_sizes(V_B, related, unrelated)
#     EQ <- mean(eff$D > 0)

#     # jackknife - type
#     V_B_jt_3d <- V_B_jack_type(dat, lambda = lambda)
#     EQ_jt0 <- numeric(length(type_names))
#     EQ_jt0[] <- NA
#     for (k in 1:(dim(V_B_jt_3d)[3])) {
#       eff_jt00 <- effect_sizes(V_B_jt_3d[,,k], related, unrelated)
#       EQ_jt0[k] <- mean(eff_jt00$D > 0)
#     }

#     eff_jt <- mean(EQ_jt0)
#     eff_jt_var <- sum((EQ_jt0  - eff_jt)^2) * ((length(EQ_jt0) - 1)/length(EQ_jt0))

#     EQ_jc0 <- numeric(length(context_names))
#     EQ_jc0[] <- NA

#     # jackknife - contexts
#     for (i in 1:length(context_names)) {
#       related_tmp <- matrix( related[rowSums(related == context_names[i]) == 0, ]   , ncol = 2)
#       unrelated_tmp <- matrix( unrelated[rowSums(unrelated == context_names[i]) == 0, ] , ncol = 2)
#       if (nrow(related_tmp) > 0 & nrow(unrelated_tmp) > 0) {
#         eff_tmp <- effect_sizes(V_B, related_tmp, unrelated_tmp)
#       }
#       EQ_jc0[i] <- mean(eff_tmp$D >0)
#     }

#     eff_jc <- mean(EQ_jc0)
#     eff_jc_var <- sum((EQ_jc0  - eff_jc)^2)  * ((length(EQ_jc0) - 1)/length(EQ_jc0))


#   list(EQ = EQ, EQ_jt = eff_jt, EQ_jt_var = eff_jt_var, EQ_jc = eff_jc, EQ_jc_var = eff_jc_var)
#   }
# }




#' Trace Coefficient
#'
#' Computes the trace of the 
#' 
#' @param x A data frame in which the first column is labeled \code{Context} and the second column is labeled \code{Type}, with each row containing a type present in that context.
#' @returns A three-dimensional array of the pairwise context-by-context effect sizes given the ommission of a type, for all types given in the input data frame.
#' 
#' @export
log_OR <- function(x) {
    UseMethod("log_OR")
}

#' @rdname log_OR
#' @export
log_OR.matrix <- function(x) {
  x <- x + 0.5
  log_or <- log(x[1,1] * x[2,2] / (x[1,2] * x[2,1]))
  return(log_or)
}






#' Pairwise Log Odds Ratios with Haldane-Anscombe's Correction
#'
#' Selecting all pairs of columns (contexts) of a contingency table, computes the log odds ratio of a 2 x 2 contingency table that indicates the number of matching types (rows) that are present. The addition of 0.5 to cells of the contingency table is performed  \insertCite{anscombe_estimating_1956,haldane_estimation_1956}{arkhaia}.
#' 
#' @param x A data frame in which the first column is labeled \code{Context} and the second column is labeled \code{Type}, with each row containing a type present in that context.
#' @returns A matrix of the pairwise log odds ratios giving the overlap in types between contexts. 
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
      log_or <- log(PA[1,1] * PA[2,2] / (PA[1,2] * PA[2,1]))
      X[i,j] <- log_OR(PA)
    }
  }

  rownames(X) <- colnames(x)
  colnames(X) <- colnames(x)

  return(X)
}





#' Pairwise Trace Coefficient
#'
#' Selecting all pairs of columns (contexts) of a contingency table, computes the trace of a 2 x 2 contingency table that indicates the number of matching types (rows) that are present.
#' 
#' @param x A data frame in which the first column is labeled \code{Context} and the second column is labeled \code{Type}, with each row containing a type present in that context.
#' @returns A matrix of the pairwise trace giving the overlap in types between contexts. 
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

  return(X)
}




#' Effect Sizes for Related vs. Unrelated Contexts
#'
#' Selecting all pairs of columns (contexts) of a contingency table, computes the trace of a 2 x 2 contingency table that indicates the number of matching types (rows) that are present.
#' 
#' @param x A data frame in which the first column is labeled \code{Context} and the second column is labeled \code{Type}, with each row containing a type present in that context.
#' @param related The related pairs of contexts.
#' @param unrelated The unrelated pairs of contexts.
#' @param direction Whether the related or unrelated effect size should come first. Default is \code{"UW"}; alternative is \code{"WU"}.
#' @returns A list of results.
#' 
#' @export
effect_sizes <- function(x, related = NULL, unrelated = NULL, direction = "UW") {
    UseMethod("effect_sizes")
}
#' 
#' @rdname effect_sizes
#' @export
effect_sizes <- function(x, related = NULL, unrelated = NULL, direction = "UW") {

  if (is.null(related) | is.null(unrelated)) {
    stop("Names of ")
  }

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

  return(res)
}








#' Presence-Absence Matrix
#'
#' Create a 2 x 2 contingency table of the presence/absence of a given type 
#' 
#' @param x A matrix or data frame where the two columns indicate contexts and the rows indicate types, with each cell containing counts of types in a context. 
#' #' @examples 
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
pa_matrix.data.fame <- function(x) {
  x <- as.matrix(x)
  out <- pa_matrix.matrix(x)
  return(out)
}
