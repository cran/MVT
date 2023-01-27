## ID: kurtosis.R, last updated 2023-01-25, F.Osorio

kurtosis.student <- function(x)
{ ## Mardia's multivariate kurtosis coefficient
  b2 <- kurtosis(x)
  attr(b2, 'excess') <- NULL
  # estimation of the kurtosis parameter for the t-distribution
  p <- ncol(x)
  k <- b2 / (p * (p + 2)) - 1
  if (k < 0)
    k <- NaN
  eta <- .5 * k / (1 + 2 * k)
  list(kurtosis = b2, kappa = k, eta = eta)
}
