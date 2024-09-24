## ID: RNG.R, last updated 2024-09-23, F.Osorio

rmt <-
function(n = 1, mean = rep(0, nrow(Sigma)), Sigma = diag(length(mean)), eta = .25)
{ ## multivariate Student-t random deviates
  if (n <= 0)
    stop("n must be a positive integer")
  if (length(mean) != nrow(Sigma))
    stop("mean and Sigma have non-conforming sizes")
  if ((eta < 0 ) || (eta >= 1/2))
    stop("eta must be in [0,1/2)")
  p <- nrow(Sigma)

  y <- matrix(0, nrow = n, ncol = p)
  dy <- dim(y)
  if (eta == 0) { # from fastmatrix
    y <- rmnorm(n, mean, Sigma)
  } else { # call C code
    y <- .C("RNG_mstudent",
            y = as.double(y),
            dims = as.integer(dy),
            mean = as.double(mean),
            Sigma = as.double(Sigma),
            eta = as.double(eta))$y
    y <- matrix(y, nrow = n, byrow = TRUE)
  }
  y
}
