\name{summary.GLS.ret}
\alias{summary.GLS.ret}

\title{Summary of detection result}

\description{ 
 Summarize a result object exported from Gls method.
} 

\usage{ 
\method{summary}{GLS.ret}(object, ...) 
}

\arguments{
\item{object}{an object return by \code{\link{gls.simple}} or \code{\link{gls.plink}}or \code{\link{gls.snpmat}}}
\item{...}{Other parameters}
}

\details{
The summary command gives some tables of the result object r.gls which contains multiple values, and the command will export all of this values to the to the console.
}

\author{Zhong Wang, \email{zzw2@psu.edu} }

\examples{

##r.gls is a result object obtained from gls.simpe, gls.plink, gls.snpmat
##e.g. r.gls <- gls.simple(...);
#
# summary(r.gls);

}

\keyword{summary}