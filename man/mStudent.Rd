\name{mStudent}
\alias{mStudent}
\alias{dmt}
\alias{rmt}

\title{Multivariate t distribution}
\description{
  These functions provide the density and random number generation from the multivariate 
  Student-t distribution.
}
\usage{
dmt(x, mean = rep(0, nrow(Sigma)), Sigma = diag(length(mean)), eta = 0.25, log = FALSE)
rmt(n = 1, mean = rep(0, nrow(Sigma)), Sigma = diag(length(mean)), eta = 0.25)
}

\arguments{
  \item{x}{ vector or matrix of data. }
  \item{n}{ the number of samples requested.}
  \item{mean}{ a vector giving the means of each variable}
  \item{Sigma}{ a positive-definite covariance matrix}
  \item{eta}{ shape parameter (must be in \code{[0,1/2)}). Default value is \code{0.25}}
  \item{log}{ logical; if TRUE, the logarithm of the density function is returned.}
}

\value{
  If \code{x} is a matrix with \eqn{n} rows, then \code{dmt} returns a \eqn{n\times 1}{n x 1} 
  vector considering each row of \code{x} as a copy from the multivariate \emph{t} distribution.

  If \code{n = 1}, then \code{rmt} returns a vector of the same length as \code{mean}, otherwise 
  a matrix of \code{n} rows of random vectors.
}

\details{
  A random vector \eqn{\bold{X} = (X_1,\dots,X_p)^T} has a multivariate \emph{t} distribution, 
  with a \eqn{\bold{\mu}} mean vector, covariance matrix \eqn{\bold{\Sigma}}, and \eqn{0 \leq \eta 
  < 1/2} shape parameter, if its density function is given by
  \deqn{
    f(\bold{x}) = K_p(\eta)|\bold{\Sigma}|^{-1/2}\left\{1 + c(\eta)(\bold{x} - \bold{\mu})^T 
    \bold{\Sigma}^{-1} (\bold{x} - \bold{\mu})\right\}^{-\frac{1}{2\eta}(1 + \eta p)}.
  }
  where
  \deqn{
    K_p(\eta) = \left(\frac{c(\eta)}{\pi}\right)^{p/2}\frac{\Gamma(\frac{1}{2\eta}(1 + \eta p))}
    {\Gamma(\frac{1}{2\eta})},
  }
  with \eqn{c(\eta)=\eta/(1 - 2\eta)}. This parameterization of the multivariate \emph{t} distribution
  is introduced mainly because \eqn{\bold{\mu}} and \eqn{\bold{\Sigma}} correspond to the mean vector 
  and covariance matrix, respectively.

  The function \code{rmt} is an interface to C routines, which make calls to subroutines from LAPACK. 
  The matrix decomposition is internally done using the Cholesky decomposition. If \code{Sigma} is not 
  non-negative definite then there will be a warning message.
  
  This parameterization of the multivariate-\emph{t} includes the normal distribution as a particular 
  case when \code{eta = 0}.
}

\references{
  Fang, K.T., Kotz, S., Ng, K.W. (1990).
  \emph{Symmetric Multivariate and Related Distributions}.
  Chapman & Hall, London.

  Gomez, E., Gomez-Villegas, M.A., Marin, J.M. (1998).
  A multivariate generalization of the power exponential family of distributions.
  \emph{Communications in Statistics - Theory and Methods} \bold{27}, 589-600.
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
