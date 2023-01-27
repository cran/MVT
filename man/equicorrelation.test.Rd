\name{equicorrelation.test}
\alias{equicorrelation.test}
\title{Equicorrelation test}
\description{
Performs several test for testing that the covariance matrix follows an equicorrelation (or compound symmetry) structure.
Likelihood ratio test, score, Wald and gradient can be used as a test statistic.
}
\usage{
equicorrelation.test(object, test = "LRT")
}
\arguments{
  \item{object}{ object of class \code{'studentFit'} representing the fitted model.}
  \item{test}{ test statistic to be used. One of "LRT" (default), "Wald", "score" or "gradient".}
}
\value{
  A list of class 'equicorrelation.test' with the following elements:

  \item{statistic}{ value of the statistic, i.e. the value of either Likelihood ratio test, Wald, score or gradient test.}
  \item{parameter}{ the degrees of freedom for the test statistic, which is chi-square distributed.}
  \item{p.value}{ the p-value for the test.}
  \item{estimate}{ the estimated covariance matrix.}
  \item{null.value}{ the hypothesized value for the covariance matrix.}
  \item{method}{ a character string indicating what type of test was performed.}
  \item{null.fit}{ a list representing the fitted model under the null hypothesis.}
  \item{data}{ name of the data used in the test.}
}
\references{
  Sutradhar, B.C. (1993).
  Score test for the covariance matrix of the elliptical t-distribution.
  \emph{Journal of Multivariate Analysis} \bold{46}, 1-12.
}
\examples{
data(examScor)
fit <- studentFit(examScor, family = Student(eta = .25))
fit

z <- equicorrelation.test(fit, test = "LRT")
z
}
\keyword{htest}
