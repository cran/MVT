\name{kurtosis}
\alias{kurtosis}
\title{Mardia's multivariate kurtosis coefficient}
\description{
This function computes the kurtosis of a multivariate distribution and estimates the kurtosis parameter for the t-distribution using the method of moments.
}
\usage{
kurtosis(x, center, cov)
}
\arguments{
 \item{x}{ vector or matrix of data with, say, p columns.}
 \item{center}{ mean vector of the distribution or second data vector of length p.}
 \item{cov}{ covariance matrix (p x p) of the distribution.}
}
\value{
 A list with the following components :

 \item{kurtosis}{ returns the value of Mardia's multivariate kurtosis.}
 \item{kappa}{ returns the excess kurtosis related to a multivariate t-distribution.}
 \item{eta}{ estimated shape (kurtosis) parameter using the methods of moments, only valid if \eqn{0 \le \eta < 1/4}.}
}
\references{
  Mardia, K.V. (1970).
  Measures of multivariate skewness and kurtosis with applications.
  \emph{Biometrika} \bold{57}, 519-530.

  Osorio, F., and Galea, M. (2015).
  Statistical inference in multivariate analysis using the t-distribution.
  Unpublished manuscript.
}
\examples{
data(companies)
S <- cov(companies)
kurtosis(companies, colMeans(companies), S)
}
\keyword{multivariate}
