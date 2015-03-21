\name{summary.BLS.ret}
\alias{summary.BLS.ret}

\title{Summary of detection result}

\description{ 
 Summarize a result object exported from Bls method.
} 

\usage{ 
\method{summary}{BLS.ret}(object, ...) 
}

\arguments{
\item{object}{an object return by \code{\link{bls.simple}} or \code{\link{bls.plink}}or \code{\link{bls.snpmat}}}
\item{...}{Other parameters}
}

\details{
The summary command gives some tables of the result object r.bls which contains multiple values, and the command will export all of this values to the to the console.
}

\author{Zhong Wang, \email{zzw2@psu.edu} }

\seealso{See \code{\link{bls.simple}} or \code{\link{bls.plink}}or \code{\link{bls.snpmat}} to generate an object }

\examples{

##r.bls is a result object obtained from bls.simpe, bls.plink, bls.snpmat
##e.g. r.bls <- bls.simple(...);
#
# summary(r.bls);

}

\keyword{summary}