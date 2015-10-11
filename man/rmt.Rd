\name{rmt}
\alias{rmt}

\title{Multivariate-t random deviates}
\description{
  Random number generation from the multivariate-t distribution.
}
\usage{
rmt(n = 1, mean = rep(0, nrow(Sigma)), Sigma = diag(length(mean)), eta = .25)
}

\arguments{
  \item{n}{ the number of samples requested}
  \item{mean}{ a vector giving the means of each variable}
  \item{Sigma}{ a positive-definite covariance matrix}
  \item{eta}{ shape parameter (must be in \code{[0,1/2)}). Default value is \code{0.25}}
}

\value{
  If \code{n = 1} a vector of the same length as \code{mean}, otherwise a
  matrix of \code{n} rows of random vectors.
}

\details{
  The function \code{rmt} is an interface to C routines, which make calls to
  subroutines from LAPACK. The matrix decomposition is internally done using
  the Cholesky decomposition. If \code{Sigma} is not non-negative definite then
  there will be a warning message.
  
  This parameterization of the multivariate-t includes the normal distribution as a particular
  case when \code{eta = 0}.
}

\references{
  Devroye, L. (1986).
  \emph{Non-Uniform Random Variate Generation}.
  Springer-Verlag, New York.
}

\seealso{
  \code{\link{rt}}
}

\examples{
# covariance matrix
Sigma <- matrix(c(10,3,3,2), ncol = 2)
Sigma

# generate the sample
y <- rmt(n = 1000, Sigma = Sigma)

# scatterplot of a random bivariate t sample with mean vector
# zero and covariance matrix 'Sigma'
par(pty = "s")
plot(y, xlab = "", ylab = "")
title("bivariate t sample (eta = 0.25)", font.main = 1)
}
\keyword{distribution}
\keyword{multivariate}
