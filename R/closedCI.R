#' Multi-visit closed population mark-recapture estimates
#'
#' @description Calculate adjusted Schnabel and Schumacher-Eschmeyer estimates with confidence intervals.
#'
#' @param marked number of animals marked on first visit (\code{M2})
#' @param caught vector of catch on subsequent visits (\code{nk})
#' @param recaptured vector of recaptures on subsequent visits (\code{mk})
#' @param newmarks vector of newly marked animals on subsequent visits (default: \code{nk - mk})
#' @param alpha type I error rate for confidence intervals (default: 0.05)
#' @param ndraws number of bootstrap draws (default: 100,000)
#'
#' @details Bias adjusted estimators are based on Dettloff (2023).
#' Bootstrap confidence intervals are computed using a beta-binomial distribution with
#' \code{n = nk}, \code{alpha = mk}, and \code{beta = nk - mk}.
#'
#' @return Matrix containing population size estimates with confidence intervals for each method.
#'
#' @references Dettloff, K. (2023). Assessment of bias and precision among simple closed
#' population mark-recapture estimators. \emph{Fisheries Research}, 265, 106756.
#' \doi{10.1016/j.fishres.2023.106756}
#'
#' @importFrom stats rbeta rbinom quantile median
#' @export
#' @examples
#' M2 <- 2
#' n <- c(232, 524, 152, 98, 353)
#' m <- c(0, 5, 8, 6, 13)
#' set.seed(123)
#' closedCI(M2, n, m, ndraws = 10000)

closedCI <- function(marked, caught, recaptured, newmarks = NULL, alpha = 0.05, ndraws = 1e5) {

  # validate inputs
  stopifnot(
    length(marked) == 1,
    length(caught) == length(recaptured),
    alpha > 0, alpha < 1,
    ndraws >= 1e3,
    ndraws %% 1 == 0,
    all(c(marked, caught, recaptured) %% 1 == 0),
    marked > 0,
    all(caught > 0),
    all(recaptured >= 0),
    marked >= recaptured[1],
    all(recaptured <= caught),
    any(recaptured > 0)
  )

  # take newmarks to be caught minus recaptured by default
  if (!is.null(newmarks)) {
    stopifnot(
      all(newmarks %% 1 == 0),
      length(newmarks) == length(caught) - 1,
      all(newmarks >= 0),
      all(newmarks <= caught[-length(caught)] - recaptured[-length(recaptured)])
    )
  } else {
    if (length(caught) > 1) {
      newmarks <- caught[-length(caught)] - recaptured[-length(recaptured)]
    } else {
      newmarks <- numeric(0)
    }
  }

  # calculate marked individuals at large for sample
  M <- cumsum(c(marked, newmarks))
  # ensure number of recaptures is not greater than marked individuals at large
  stopifnot(all(recaptured <= M))

  # force data types to numeric instead of integer
  marked <- as.numeric(marked)
  caught <- as.numeric(caught)
  recaptured <- as.numeric(recaptured)

  ## create matrices to store bootstrap replicates
  n_p <- length(caught)

  # draw number of recaptures from beta-binomial distribution
  mk_vec <- rbinom(ndraws * n_p, rep(caught, each = ndraws),
                   rbeta(ndraws * n_p, rep(recaptured, each = ndraws), rep(caught - recaptured, each = ndraws)))
  mk_mat <- matrix(mk_vec, nrow = ndraws, ncol = n_p)
  # ensure number of recaptures is not greater than marked individuals at large in sample
  mk_mat[] <- pmin(mk_vec, rep(M, each = ndraws))

  # calculate marked individuals at large for random draw
  if (n_p > 1) {
    diffs <- rep(caught[-n_p], each = ndraws) - mk_mat[, -n_p, drop = FALSE]
    cum_diffs <- matrix(0, nrow = ndraws, ncol = n_p - 1)
    cum_diffs[, 1] <- diffs[, 1]
    if (n_p > 2) {
      for (j in 2:(n_p - 1)) {
        cum_diffs[, j] <- cum_diffs[, j - 1] + diffs[, j]
      }
    }
    Mk_mat <- cbind(marked, marked + cum_diffs)
  } else {
    Mk_mat <- matrix(marked, nrow = ndraws, ncol = 1)
  }

  # calculate adjusted Schnabel statistic
  S1_reps <- pmax(ceiling(as.vector(Mk_mat %*% caught) / (rowSums(mk_mat) + 1)), M[length(M)])

  # calculate adjusted Schumacher-Eschmeyer statistic
  penalty <- prod(c(M[1], caught))^-(1 / (n_p + 1)) # decimal correction for adjusted Schumacher-Eschmeyer estimator
  S2_reps <- pmax(as.vector(((Mk_mat + 1)^2) %*% (caught + 1)) / rowSums(Mk_mat * (mk_mat + 1)) - 2 - penalty, M[length(M)])

  # calculate adjusted Schnabel population size estimate
  estS1 <- max(ceiling(sum(M * caught) / (sum(recaptured) + 1)), M[length(M)])
  # calculate adjusted Schumacher-Eschmeyer population size estimate
  estS2 <- max(round(sum((M + 1)^2 * (caught + 1)) / sum(M * (recaptured + 1)) - 2 - penalty), M[length(M)])

  # store bootstrap replicates
  reps <- list(S1_reps, S2_reps)

  # compute bias adjusted percentile confidence intervals
  ci <- mapply(function(x, y) {
    # calculate quantiles
    q_raw <- quantile(x, probs = c(alpha / 2, 1 - alpha / 2))
    # adjust bias
    raw_ci <- round(q_raw - (median(x) - y))
    # clamp values
    lower <- pmax(pmin(min(raw_ci), y), M[length(M)])
    upper <- pmax(max(raw_ci), y)
    # combine and add names
    res <- c(lower, upper)
    names(res) <- names(q_raw)
    return(res)
  }, x = reps, y = c(estS1, estS2))

  # output matrix
  out <- cbind(Estimate = c(estS1, estS2), t(ci))
  rownames(out) <- c("Schnabel", "S-E")
  return(out)
}
