\name{companies}
\docType{data}
\alias{companies}

\title{Financial data}
\description{
  Data extracted from Standard & Poor's Compustat PC Plus. This dataset has been used to illustrate some influence diagnostic techniques.
}
\usage{data(companies)}
\format{
  A data frame with 26 observations on the following 3 variables.
  \describe{
    \item{book}{ book value in dollars per share at the end of 1992.}
    \item{net}{ net sales in millions of dollars in 1992.}
    \item{ratio}{ sales to assets ratio in 1992.}
  }
}
\source{
  Hadi, A.S., and Nyquist, H. (1999).
  Frechet distance as a tool for diagnosing multivariate data.
  \emph{Linear Algebra and Its Applications} \bold{289}, 183-201.
  
  Hadi, A.S., and Son, M.S. (1997).
  Detection of unusual observations in regression and multivariate data.
  In: A. Ullah, D.E.A. Giles (Eds.) \emph{Handbook of Applied Economic Statistics}.
  Marcel Dekker, New York. pp. 441-463.
}
\keyword{datasets}