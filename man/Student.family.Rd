\name{Student}
\alias{Student}
\alias{Student.family}

\title{Family object for the multivariate t-distribution}
\description{
  Provide a convenient way to specify the details of the model used by function \code{\link{studentFit}}.
}
\usage{
Student(eta = .25)
}

\arguments{
  \item{eta}{ shape parameter for the multivariate t-distribution, must be confined to \code{[0,1/2)}.}
}
\details{
  \code{Student} is a generic function to create info about the t-distribution which
  is passed to the estimation algorithm.
}

\examples{
MyFmly <- Student(eta = .4)
MyFmly
}
\keyword{models}
