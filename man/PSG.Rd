\name{PSG}
\docType{data}
\alias{PSG}

\title{Transient sleep disorder}
\description{
  Clinical study designed to compare the automated and semi-automated scoring of 
  Polysomnographic (PSG) recordings used to diagnose transient sleep disorders. 
  The study considered 82 patients who were given a sleep-inducing drug (Zolpidem 
  10 mg). Measurements of latency to persistent sleep (LPS: lights out to the 
  beginning of 10 consecutive minutes of uninterrupted sleep) were obtained using 
  six different methods.
}
\usage{data(PSG)}
\format{
  A data frame with 82 observations on the following 3 variables.
  \describe{
    \item{manual}{ fully manual scoring.}
    \item{automated}{ automated scoring by the Morpheus software.}
    \item{partial}{ Morpheus automated scoring with manual review.}
  }
}
\source{
  Svetnik, V., Ma, J., Soper, K.A., Doran, S., Renger, J.J., Deacon, S., Koblan, K.S. (2007).
  Evaluation of automated and semi-automated scoring of polysomnographic recordings from a clinical trial using zolpidem in the treatment of insomnia.
  \emph{SLEEP} \bold{30}, 1562-1574.
}
\keyword{datasets}