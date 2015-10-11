\name{commutation}
\alias{commutation}

\title{Commutation matrix for square matrices}
\description{
  This function returns a commutation matrix of order \code{n} which transforms,
  for an \code{n} by \code{n} matrix \code{x}, \code{vec(x)} into \code{vec(t(x))}.
}
\usage{
commutation(n = 2)
}
\arguments{
  \item{n}{ a non-negative integer.}
}
\value{
  Returns an square commutation matrix of order \eqn{n^2} .
}
\references{
  Magnus, J.R., and Neudecker, H. (1979).
  The commutation matrix: some properties and applications.
  \emph{The Annals of Statistics} \bold{7}, 381-394.

  Magnus, J.R., and Neudecker, H. (1999).
  \emph{Matrix Differential Calculus with Applications in Statistics and Econometrics}, 2nd Edition.
  Wiley, New York.
}
\examples{
a <- matrix(1:9, ncol = 3)
kmat <- commutation(nrow(a))
vec <- as.vector(a)
## vectorization of t(a):
kmat \%*\% vec
}
\keyword{array}
\keyword{algebra}
