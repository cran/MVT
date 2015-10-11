\name{Weights}
\alias{Weights}
\alias{dweights}
\alias{pweights}
\alias{qweights}
\title{Distribution of the weights from a multivariate t-distribution}
\description{
  Density, distribution function and quantile function for the weights distribution arising from the multivariate t-distribution with dimension \code{dim} and
  shape parameter \code{eta}.
}
\usage{
dweights(x, eta = .25, dim, log = FALSE, scaled = TRUE)
pweights(q, eta = .25, dim, lower.tail = TRUE, log.p = FALSE, scaled = TRUE)
qweights(p, eta = .25, dim, lower.tail = TRUE, log.p = FALSE, scaled = TRUE)
}
\arguments{
  \item{x, q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{eta}{shape parameter of the multivariate t-distribution, must be in the interval \eqn{[0,1/2)}. Default value is \code{eta = 0.25}.}
  \item{dim}{dimension of the multivariate t-distribution.}
  \item{log, log.p}{logical; if TRUE, probabilities p are given as \eqn{\log(p)}.}
  \item{lower.tail}{logical; if TRUE (default), probabilities are \eqn{P(X \le x)}, otherwise, \eqn{P(X > x)}.}
  \item{scaled}{logical; if TRUE, the weights are scaled to belong the interval \code{(0,1)}.}
}
\details{
  The weights' distribution with parameters \code{eta} and \code{dim} \eqn{= p} has density
  \deqn{f(x) = \frac{m^{1-(1/\eta + p)/2}}{B(1/(2\eta),p/2)}{x}^{1/(2\eta)-1}{(m-x)}^{p/2-1}%
  }{f(x) = m^(1-(1/\eta + p)/2)/B(1/(2\eta),p/2)x^(1/(2\eta)-1)(m-x)^(p/2-1)}
  for \eqn{0 \le \eta < 1/2}, \eqn{p > 0} and \eqn{0 < x < m}, where \eqn{m = (1 + p\eta)/(1 - 2\eta)}.
  
  The mean is \eqn{E(X) = 1/(1-2\eta)} and the variance is
  \deqn{Var(X) = \frac{2p\eta^2}{(1+(p+2)\eta){(1-2\eta)}^2}}{2p\eta^2/((1+(p+2)\eta)(1-2\eta)^2)}
  
  The scaled version of the weights distribution has a Beta distribution with parameters \eqn{1/(2\eta)} and \eqn{p/2}.
}
\value{
  \code{dweights} gives the density, \code{pweights} the distribution function,
  and \code{qweights} the quantile function.

  Invalid arguments will result in return value \code{NaN}, with a warning.
}
\references{
  Abramowitz, M., and Stegun, I.A. (1972). 
  \emph{Handbook of Mathematical Functions}. Dover, New York.
  Chapter 6: Gamma and related functions.

  Johnson, N.L., Kotz, S., and Balakrishnan, N. (1995). 
  \emph{Continuous Univariate Distributions}, volume 2, 2nd Ed. Wiley, New York.
  Chapter 25: Beta distributions.
  
  Osorio, F., and Galea, M. (2015).
  Statistical inference in multivariate analysis using the t-distribution.
  Unpublished manuscript.
}
\seealso{
  \link{Distributions} for other standard distributions.

  \code{\link{beta}} for the Beta function.
}
\examples{
data(companies)
fit <- studentFit(companies, family = Student(eta = .25))

# compute the 5% quantile from the estimated distribution of the weights
p <- fit$dims[2]
eta <- fit$eta
wts <- fit$weights
cutoff <- qweights(.05, eta = eta, dim = p, scaled = FALSE)

# identify observations with 'small' weights
n <- fit$dims[1]
which <- seq_len(n)[wts < cutoff]
which
}
\keyword{distribution}
