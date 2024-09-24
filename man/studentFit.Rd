\name{studentFit}
\alias{studentFit}

\title{Estimation of mean and covariance using the multivariate t-distribution}
\description{
  Estimates the mean vector and covariance matrix assuming the data came from a multivariate
  t-distribution: this provides some degree of robustness to outlier without giving a high
  breakdown point.
}
\usage{
studentFit(x, data, family = Student(eta = .25), covStruct = "UN", subset, na.action, 
control)
}
\arguments{
  \item{x}{  a formula or a numeric matrix or an object that can be coerced to a numeric matrix.}
  \item{data}{  an optional data frame (or similar: see \code{\link{model.frame}}), used only if
  \code{x} is a formula. By default the variables are taken from \code{environment(formula)}.}
  \item{family}{ a description of the error distribution to be used in the model.
    By default the multivariate t-distribution with 0.25 as shape parameter is considered
    (using \code{eta = 0} allows to tackle the multivariate normal distribution).}
  \item{covStruct}{ a character string specifying the type of covariance structure. The options 
    available are: \code{"UN"} (unstructured) general covariance matrix with no additional structure
    (default), \code{"CS"} (compound symmetry) corresponding to a constant correlation or equicorrelation, 
    \code{"DIAG"} (diagonal) representing a diagonal positive-definite matrix, \code{"HOMO"} (homogeneous) 
    meaning a covariance matrix with homogeneous variances.}
  \item{subset}{ an optional expression indicating the subset of the rows of
    data that should be used in the fitting process.}
  \item{na.action}{ a function that indicates what should happen when the data contain NAs.}
  \item{control}{ a list of control values for the estimation algorithm to replace
    the default values returned by the function \code{\link{MVT.control}}.}
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
  Kent, J.T., Tyler, D.E., Vardi, Y. (1994).
  A curious likelihood identity for the multivariate t-distribution.
  \emph{Communications in Statistics: Simulation and Computation}
  \bold{23}, 441-453.

  Lange, K., Little, R.J.A., Taylor, J.M.G. (1989).
  Robust statistical modeling using the t distribution.
  \emph{Journal of the American Statistical Association}
  \bold{84}, 881-896.

  Osorio, F., Galea, M., Henriquez, C., Arellano-Valle, R. (2023).
  Addressing non-normality in multivariate analysis using the \emph{t}-distribution. 
  \emph{AStA Advances in Statistical Analysis} \bold{107}, 785-813.
}

\seealso{
\code{\link{cov}}, \code{\link[MASS:cov.rob]{cov.rob}} and \code{\link[MASS:cov.trob]{cov.trob}} 
in package \CRANpkg{MASS}.
}

\examples{
data(PSG)
fit <- studentFit(~ manual + automated, data = PSG, family = Student(eta = 0.25))
fit
}
\keyword{models}
\keyword{multivariate}
