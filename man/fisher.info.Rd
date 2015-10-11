\name{fisher.info}
\alias{fisher.info}
\title{Fisher information matrix}
\description{
Calculate the Fisher information matrix for the multivariate t-distribution.
}
\usage{
fisher.info(object)
}
\arguments{
\item{object}{ an object of class \code{'studentFit'} representing the fitted model.}
}
\value{
  A square matrix of order \eqn{p(p + 3)/2 + 1} containing the Fisher information of the fitted model.
}
\references{
  Osorio, F., and Galea, M. (2015).
  Statistical inference in multivariate analysis using the t-distribution.
  Unpublished manuscript.
}
\seealso{
\code{\link{studentFit}}
}
\examples{
data(PFM)
fit <- studentFit(~ cuprum.D + habitat.D + planvital.D + provida.D, data = PFM,
  family = Student(eta = 0.25))
info <- fisher.info(fit)
}
\keyword{multivariate}
