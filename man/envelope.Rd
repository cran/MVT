\name{envelope}
\alias{envelope}
\title{QQ-plot with simulated envelopes}
\description{
Constructs a normal QQ-plot using a Wilson-Hilferty transformation for the estimated Mahalanobis distances obtained from the fitting procedure.
}
\usage{
envelope(object, reps = 50, conf = 0.95, plot.it = TRUE)
}
\arguments{
\item{object}{ an object of class \code{'studentFit'} representing the fitted model.}
\item{reps}{ number of simulated point patterns to be generated when computing the envelopes. The default number is 50.}
\item{conf}{ the confidence level of the envelopes required. The default is to find 95\% confidence envelopes.}
\item{plot.it}{ if TRUE it will draw the corresponding plot, if FALSE it will only return the computed values. }
}
\value{
  A list with the following components :
  
  \item{transformed}{ a vector with the \code{z}-scores obtained from the Wilson-Hilferty transformation.}
  \item{envelope}{ a matrix with two columns corresponding to the values of the lower and upper pointwise confidence
  envelope.}
}
\references{
  Atkinson, A.C. (1985).
  \emph{Plots, Transformations and Regression}.
  Oxford University Press, Oxford.

  Osorio, F., and Galea, M. (2015).
  Statistical inference in multivariate analysis using the t-distribution.
  Unpublished manuscript.
}
\seealso{\code{\link{wilson.hilferty}}}
\examples{
data(PFM)
fit <- studentFit(~ cuprum.D + habitat.D + planvital.D + provida.D, data = PFM,
  family = Student(eta = 0.25))
envelope(fit, reps = 500, conf = 0.95)
}
\keyword{dplot}
