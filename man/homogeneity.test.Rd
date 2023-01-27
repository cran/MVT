\name{homogeneity.test}
\alias{homogeneity.test}
\title{Test of variance homogeneity of correlated variances}
\description{
Performs several test for testing equality of \eqn{p \ge 2} correlated variables. Likelihood ratio test,
score, Wald and gradient can be used as a test statistic.
}
\usage{
homogeneity.test(object, test = "LRT")
}
\arguments{
  \item{object}{ object of class \code{'studentFit'} representing the fitted model.}
  \item{test}{ test statistic to be used. One of "LRT" (default), "Wald", "score" or "gradient".}
}
\value{
  A list of class 'homogeneity.test' with the following elements:

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
  Harris, P. (1985).
  Testing the variance homogeneity of correlated variables.
  \emph{Biometrika} \bold{72}, 103-107.
  
  Modarres, R. (1993).
  Testing the equality of dependent variables.
  \emph{Biometrical Journal} \bold{7}, 785-790.
}
\examples{
data(examScor)
fit <- studentFit(examScor, family = Student(eta = .25))
fit

z <- homogeneity.test(fit, test = "LRT")
z
}
\keyword{htest}
