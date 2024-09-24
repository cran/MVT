## ID: distn.R, last updated 2024-09-23, F.Osorio

dmt <-
function(x, mean = rep(0, nrow(Sigma)), Sigma = diag(length(mean)), eta = 0.25, log = FALSE)
{ # pdf for the multivariate Student-t distribution
  give.log <- log
  if (is.data.frame(x))
    x <- as.matrix(x)
  else if (!is.matrix(x))
    stop("'x' must be a matrix or a data frame")
  if (!all(is.finite(x)))
    stop("'x' must contain finite values only")
  n <- nrow(x)
  p <- ncol(x)

  if (is.null(mean))
    stop("'mean' must be provided")
  if (isFALSE(mean))
    mean <- double(p) # center is zeroed
  if (!is.vector(mean))
    stop("'mean' must be a vector")
  if (length(mean) != p)
    stop("'mean' has incorrect length")

  if (is.null(Sigma))
    stop("'Sigma' matrix must be provided")
  if (!is.matrix(Sigma))
    stop("'Sigma' must be a matrix")
  if (all(dim(Sigma) != c(p,p)))
    stop("'Sigma' has incorrect dimensions")
  if (!isSymmetric(Sigma))
    Sigma <- asSymmetric(Sigma)

  if ((eta < 0 ) || (eta >= 1/2))
    stop("eta must be in [0,1/2)")

  storage.mode(x) <- "double"
  storage.mode(Sigma) <- "double"

  # calling C code
  if (eta == 0) { # gaussian log-density
    y <- .C("pdf_mnorm",
            y = double(n),
            x = x,
            n = as.integer(n),
            p = as.integer(p),
            mean = as.double(mean),
            Sigma = Sigma)$y
  } else { # multivariate Student-t log-density
    y <- .C("pdf_mstudent",
            y = double(n),
            x = x,
            n = as.integer(n),
            p = as.integer(p),
            mean = as.double(mean),
            Sigma = Sigma,
            eta = as.double(eta))$y
  }

  if (!give.log) y <- exp(y)
  y
}
