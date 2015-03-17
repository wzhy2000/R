\name{summary.BLS.ret}
\alias{summary.BLS.ret}

\title{Summary of detection result}

\description{ 
 Summarize a result object exported from Bls method.
} 

\usage{ summary(r.bls) }

\arguments{
\item{r.bls}{an object return by \code{\link{bls.simple}} or \code{\link{bls.plink}}or \code{\link{bls.snpmat}}}
}

\details{

The summary command gives some tables of the result object r.bls which contains multiple values, and the command will export all of this values to the to the console.
}

\author{Zhong Wang, \email{zzw2@psu.edu} }

\examples{
summary(r.bls);
}

\keyword{summary}