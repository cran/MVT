\name{studentFit}
\alias{studentFit}

\title{Estimation of mean and covariance using the multivariate t-distribution}
\description{
  Estimates the mean vector and covariance matrix assuming the data came from a multivariate
  t-distribution: this provides some degree of robustness to outlier without giving a high
  breakdown point.
}
\usage{
studentFit(x, data, family = Student(eta = .25), subset, na.action, control)
}
\arguments{
  \item{x}{  a formula or a numeric matrix or an object that can be coerced to a numeric matrix.}
  \item{data}{  an optional data frame (or similar: see \code{\link{model.frame}}), used only if
  \code{x} is a formula. By default the variables are taken from \code{environment(formula)}.}
  \item{family}{ a description of the error distribution to be used in the model. 
  By default the multivariate t-distribution with 0.25 as shape parameter is considered
  (using \code{eta = 0} allows to tackle the multivariate normal distribution).}
  \item{subset}{ an optional expression indicating the subset of the rows of 
  data that should be used in the fitting process.}
  \item{na.action}{ a function that indicates what should happen when the data contain NAs.}
  \item{control}{ a list of control values for the estimation algorithm to replace 
  the default values returned by the function \code{\link{mvt.control}}.}
}
\value{
  A list with class \code{'studentFit'} containing the following components:
  \item{call}{  a list containing an image of the \code{studentFit} call that produced the object.}
  \item{family}{  the \code{\link{Student}} object used, with the estimated shape parameters (if requested).}
  \item{center}{  final estimate of the location vector.}
  \item{Scatter}{  final estimate of the scale matrix.}
  \item{logLik}{  the log-likelihood at convergence.}
  \item{numIter}{  the number of iterations used in the iterative algorithm.}
  \item{weights}{  estimated weights corresponding to the assumed heavy-tailed distribution.}
  \item{distances}{  estimated squared Mahalanobis distances.}
  \item{eta}{  final estimate of the shape parameter, if requested.}
  Generic function \code{print} show the results of the fit.
}
\references{
  Kent, J.T., Tyler, D.E., and Vardi, Y. (1994).
  A curious likelihood identity for the multivariate t-distribution.
  \emph{Communications in Statistics: Simulation and Computation}
  \bold{23}, 441-453.

  Lange, K., Little, R.J.A., and Taylor, J.M.G. (1989).
  Robust statistical modeling using the t distribution.
  \emph{Journal of the American Statistical Association}
  \bold{84}, 881-896.
  
  Osorio, F., and Galea, M. (2015).
  Statistical inference in multivariate analysis using the t-distribution.
  Unpublished manuscript.
}

\seealso{
\code{\link{cov}}, \code{\link{cov.rob}}, \code{\link{cov.trob}}
}

\examples{
data(PFM)
fit <- studentFit(~ cuprum.D + habitat.D + planvital.D + provida.D, data = PFM,
  family = Student(eta = 0.25))
fit
}
\keyword{models}
\keyword{multivariate}
