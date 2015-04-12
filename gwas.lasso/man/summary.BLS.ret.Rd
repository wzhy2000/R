\name{summary.BLS.ret}
\alias{summary.BLS.ret}

\title{Summary of detection result obtained by the BLS model }

\description{ 
 Summarize a result object obtained by the BLS method.
} 

\usage{ 
\method{summary}{BLS.ret}(object, ...) 
}

\arguments{
\item{object}{an object obtained by \code{\link{bls.simple}} or \code{\link{bls.plink}}or \code{\link{bls.snpmat}}}
\item{...}{Other parameters}
}

\details{
The summary command gives some tables of the result object r.bls which contains multiple values, and the command will export all of this values to the to the console.
}

\author{Zhong Wang and Nating Wang }

\references{
 Li, J., Das, K., Fu, G., Li, R., & Wu, R. (2011). The Bayesian lasso for genome-wide association studies. Bioinformatics, 27(4), 516-523.
}

\seealso{See \code{\link{bls.simple}} or \code{\link{bls.plink}}or \code{\link{bls.snpmat}} to generate an object }

\examples{

##r.bls is a result object obtained from bls.simpe, bls.plink, bls.snpmat
##e.g. r.bls <- bls.simple(...);
#
# summary(r.bls);

}

\keyword{summary}
\keyword{BLS model}
