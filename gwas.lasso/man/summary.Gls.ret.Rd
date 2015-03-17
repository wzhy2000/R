\name{summary.GLS.ret}
\alias{summary.GLS.ret}

\title{Summary of detection result}

\description{ 
 Summarize a result object exported from Gls method.
} 

\usage{ summary(r.gls) }

\arguments{
\item{r.gls}{an object return by \code{\link{gls.simple}} or \code{\link{gls.plink}}or \code{\link{gls.snpmat}}}
}

\details{

The summary command gives some tables of the result object r.gls which contains multiple values, and the command will export all of this values to the to the console.
}

\author{Zhong Wang, \email{zzw2@psu.edu} }

\examples{
summary(r.gls);
}

\keyword{summary}