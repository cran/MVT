\name{kurtosis.student}
\alias{kurtosis.student}
\title{Mardia's multivariate kurtosis coefficient}
\description{
This function computes the kurtosis of a multivariate distribution and estimates the kurtosis parameter for the t-distribution using the method of moments.
}
\usage{
kurtosis.student(x)
}
\arguments{
 \item{x}{ vector or matrix of data with, say, p columns.}
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

  Osorio, F., Galea, M., Henriquez, C., Arellano-Valle, R. (2023). 
  Addressing non-normality in multivariate analysis using the t-distribution. 
  \emph{AStA Advances in Statistical Analysis} \bold{107}, 785-813.
}
\examples{
data(companies)
kurtosis.student(companies)
}
\keyword{multivariate}
