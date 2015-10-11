\name{center.test}
\alias{center.test}
\title{One-sample location test}
\description{
Performs several test for testing the mean of a multivariate-t population.
Generalized Hotelling's T-squared, likelihood ratio test, score, Wald and gradient
can be used as a test statistic.
}
\usage{
center.test(object, center, test = "LRT")
}
\arguments{
  \item{object}{ object of class \code{'studentFit'} representing the fitted model.}
  \item{center}{ a vector indicating the hypothesized value of the mean.}
  \item{test}{ test statistic to be used. One of "hotelling", "LRT" (default), "Wald", "score" or "gradient".}
}
\value{
  A list of class 'center.test' with the following elements:

  \item{statistic}{ value of the statistic, i.e. the value of either Hotelling T-squared, likelihood ratio test, Wald, score or gradient test.}
  \item{parameter}{ the degrees of freedom for the test statistic, which is chi-square distributed.}
  \item{p.value}{ the p-value for the test.}
  \item{estimate}{ the estimated mean vector.}
  \item{null.value}{ the hypothesized value for the mean vector.}
  \item{method}{ a character string indicating what type of test was performed.}
  \item{null.fit}{ a list representing the fitted model under the null hypothesis.}
  \item{data}{ name of the data used in the test.}
}
\references{
  Anderson, T.W. (2003).
  \emph{An Introduction to Multivariate Statistical Analysis}.
  Wiley, New York.

  Hotelling, H. (1931).
  The generalization of Student's ratio.
  \emph{Annals of Mathematical Statistics} \bold{2}, 360-378.
   
  Osorio, F., and Galea, M. (2015).
  Statistical inference in multivariate analysis using the t-distribution.
  Unpublished manuscript.
}
\examples{
data(cork)
x <- sweep(cork[,2:4], 1, cork[,1], "-")
colnames(x) <- c("E_N", "S_N", "W_N")
pairs(x)

fit <- studentFit(x, family = Student(eta = .25))
z <- center.test(fit, center = c(0,0,0), test = "score")
z
}
\keyword{htest}
