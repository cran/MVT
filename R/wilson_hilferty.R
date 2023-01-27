## ID: wilson_hilferty.R, last updated 2021-03-05, F.Osorio

WH.student <- function(x, center, cov, eta = 0)
{ # Wilson-Hilferty transformation
  if (missing(center) && missing(cov)) {
    if (is.null(x$distances) && is.null(x$dims))
      stop("x is not a valid object.")
    distances <- x$distances
    n <- x$dims[1]
    p <- x$dims[2]
  } else {
    distances <- Mahalanobis(x, center, cov, inverted = FALSE)
    n <- nrow(x)
    p <- ncol(x)
    if (p != nrow(cov))
      stop("covariance matrix is not compatible.")
    if (p != length(center))
      stop("center vector is not compatible.")
  }

  if (eta == 0) {
    z <- .C("Wilson_Hilferty_chisq",
            distances = as.double(distances),
            n = as.integer(n),
            p = as.integer(p),
            z = double(n))$z
  } else {
    z <- .C("Wilson_Hilferty_F",
            distances = as.double(distances),
            n = as.integer(n),
            p = as.integer(p),
            eta = as.double(eta),
            z = double(n))$z
  }
  z
}
