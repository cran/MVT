\name{WH.student}
\alias{WH.student}
\title{Wilson-Hilferty transformation}
\description{
  Returns the Wilson-Hilferty transformation of random variables with \eqn{F} distribution.
}
\usage{
WH.student(x, center, cov, eta = 0)
}
\arguments{
  \item{x}{ object of class \code{'studentFit'} from which is extracted the estimated Mahalanobis distances of the fitted model.
  Also \code{x} can be a vector or matrix of data with, say, \eqn{p} columns.}
  \item{center}{ mean vector of the distribution or second data vector of length \eqn{p}. Not required if \code{x} have class \code{'studentFit'}.}
  \item{cov}{ covariance matrix (\eqn{p} by \eqn{p}) of the distribution. Not required if \code{x} have class \code{'studentFit'}.}
  \item{eta}{ shape parameter of the multivariate t-distribution. By default the multivariate normal (\code{eta = 0}) is considered.}
}
\details{
  Let \eqn{F} the following random variable:
  \deqn{F = \frac{D^2/p}{1-2\eta}}{F = D^2/(p(1-2\eta))}
  where \eqn{D^2} denotes the squared Mahalanobis distance defined as
  \deqn{D^2 = (x - \mu)^T \Sigma^{-1} (x - \mu)}{D^2 = (x - \mu)^T \Sigma^-1 (x - \mu)}

  Thus the Wilson-Hilferty transformation is given by
  \deqn{z = \frac{(1 - \frac{2\eta}{9})F^{1/3} - (1 - \frac{2}{9p})}{(\frac{2\eta}{9}F^{2/3} + \frac{2}{9p})^{1/2}}%
  }{z = ((1 - 2\eta/9)F^1/3 - (1 - 2/(9p))) / (2\eta/9 F^2/3 + 2/(9p))^1/2}
  and \eqn{z} is approximately distributed as a standard normal distribution. This is useful, for instance, in the construction of
  QQ-plots.

  For \code{eta = 0}, we obtain
  \deqn{z = \frac{F^{1/3} - (1 - \frac{2}{9p})}{(\frac{2}{9p})^{1/2}}%
  }{z = (F^1/3 - (1 - 2/(9p))) / (2/(9p))^1/2}
  which is the Wilson-Hilferty transformation for chi-square variables.
}
\references{
  Osorio, F., Galea, M., Henriquez, C., Arellano-Valle, R. (2023).
  Addressing non-normality in multivariate analysis using the \emph{t}-distribution. 
  \emph{AStA Advances in Statistical Analysis} \bold{107}, 785-813.

  Wilson, E.B., and Hilferty, M.M. (1931).
  The distribution of chi-square.
  \emph{Proceedings of the National Academy of Sciences of the United States of America} \bold{17}, 684-688.
}
\seealso{\code{\link{cov}}, \code{\link{mahalanobis}}, \code{\link{envelope.student}}}
\examples{
data(companies)
x <- companies
z <- WH.student(x, center = colMeans(x), cov = cov(x))
par(pty = "s")
qqnorm(z, main = "Transformed distances Q-Q plot")
abline(c(0,1), col = "red", lwd = 2)
}
\keyword{multivariate}
