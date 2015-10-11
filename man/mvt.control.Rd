\name{mvt.control}
\alias{mvt.control}

\title{Set control parameters}
\description{
  Allows users to set control parameters for the estimation routine available in \code{MVT}.
}
\usage{
mvt.control(maxIter = 2000, tolerance = 1e-6, fix.shape = FALSE)
}

\arguments{
  \item{maxIter}{  maximum number of iterations. The default is 2000.}
  \item{tolerance}{  the relative tolerance in the iterative algorithm.}
  \item{fix.shape}{  whether the shape parameter should be kept fixed in
  the fitting processes. The default is \code{fix.shape = FALSE}.}
}

\value{
  A list of control arguments to be used in a call to \code{\link{studentFit}}.

  A call to \code{mvt.control} can be used directly in the \code{control} argument 
  of the call to \code{\link{studentFit}}.
}

\examples{
ctrl <- mvt.control(maxIter = 500, tol = 1e-04, fix.shape = TRUE)
data(PFM)
studentFit(~ cuprum.D + habitat.D + planvital.D + provida.D, data = PFM,
  family = Student(eta = 0.25), control = ctrl)
}
\keyword{models}
