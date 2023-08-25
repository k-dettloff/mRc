#' Multi-visit closed population mark-recapture estimates
#'
#' @description Calculate adjusted Schnabel and Schumacher-Eschmeyer estimates with confidence intervals.
#'
#' @param marked number of animals marked on first visit (M2)
#' @param caught vector of catch on subsequent visits (nk)
#' @param recaptured vector of recaptures on subsequent visits (mk)
#' @param newmarks vector of newly marked animals on subsequent visits (default: nk-mk)
#' @param alpha type I error rate for confidence intervals (default: 0.05)
#' @param ndraws number of bootstrap draws (default: 10,000)
#'
#' @details Bias adjusted estimators are based on Dettloff (2023).
#' Bootstrap confidence intervals are computed using a beta-binomial distribution with n = nk, alpha = mk, beta = nk-mk.
#' @return Matrix containing population size estimates with confidence intervals for each method
#' @references Dettloff, K. (2023). Assessment of bias and precision among simple closed population mark-recapture estimators.
#' Fisheries Research 265, 106756. doi: <https://doi.org/10.1016/j.fishres.2023.106756>
#' @export
#' @examples
#' M2 = 2
#' n = c(232, 524, 152, 98, 353)
#' m = c(0, 5, 8, 6, 13)
#' closedCI(M2, n, m, ndraws = 1000)

closedCI = function(marked, caught, recaptured, newmarks = NULL, alpha = 0.05, ndraws = 1e5) {

  # validate inputs
  stopifnot(length(marked) == 1, length(caught) == length(recaptured),
            marked > 0, all(caught > 0), all(recaptured >= 0), any(recaptured > 0), all(recaptured <= caught),
            alpha > 0, alpha < 1,
            ndraws >= 1e3)

  # take newmarks to be caught minus recaptured by default
  if(!is.null(newmarks)) {
    stopifnot(length(newmarks) == length(caught) - 1,
              all(newmarks >= 0), all(newmarks <= caught[-length(caught)] - recaptured[-length(recaptured)]))
  } else newmarks = caught[-length(caught)] - recaptured[-length(recaptured)]

  # calculate marked individuals at large for sample
  M = cumsum(c(marked, newmarks))
  # ensure number of recaptures is not greater than marked individuals at large
  stopifnot(all(recaptured <= M))

  ### function to generate bootstrap replicates
  bootCI = function(M2, n, m) {

    # draw number of recaptures from beta-binomial distribution,
    # ensure number of recaptures is not greater than marked individuals at large in sample
    mk = pmin(mapply(function(n, m) stats::rbinom(1, n, stats::rbeta(1, m, n - m)), n, m), M)
    # calculate marked individuals at large for random draw
    Mk = cumsum(c(M2, n[-length(n)] - mk[-length(mk)]))
    # calculate adjusted Schnabel statistic
    S1 = ceiling(sum(Mk * n) / (sum(mk) + 1))
    # calculate adjusted S-E statistic
    S2 = sum((Mk + 1)^2 * (n + 1)) / sum(Mk * (mk + 1)) - 2
    # output vector
    c("Schnabel" = S1, "S.E." = S2)

  }

  # calculate adjusted Schnabel population size estimate
  estS1 = ceiling(sum(M * caught) / (sum(recaptured) + 1))
  # calculate adjusted S-E population size estimate
  estS2 = sum((M + 1)^2 * (caught + 1)) / sum(M * (recaptured + 1)) - 2
  # generate bootstrap replicates
  reps = replicate(ndraws, bootCI(marked, caught, recaptured))
  # compute bias adjusted percentile confidence intervals
  ci = mapply(function(x, y) pmax(round(stats::quantile(x, probs = c(alpha/2, 1 - alpha/2)) - (stats::median(x) - y)), M[length(M)]),
              x = asplit(reps, 1), y = c(estS1, estS2))
  # output matrix
  cbind(Estimate = c(estS1, round(estS2)), t(ci))

}
