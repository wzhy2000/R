\name{summary.GLS.ret}
\alias{summary.GLS.ret}

\title{Summary of detection result obtained by the GLS model }

\description{ 
 Summarize a result object obtained by the GLS method.
} 


\usage{ 
\method{summary}{GLS.ret}(object, ...) 
}

\arguments{
\item{object}{an object obtained by \code{\link{gls.simple}} or \code{\link{gls.plink}}or \code{\link{gls.snpmat}}}
\item{...}{Other parameters}
}

\details{
The summary command gives some tables of the result object r.gls which contains multiple values, and the command will export all of this values to the to the console.
}

\author{Zhong Wang and Nating Wang }

\references{
 Li. J., Wang, Z., Li, R., Wu, R.(2015). Bayesian Group Lasso for nonparametric varying-coefficient models with application to functional genome-wide association studies. The Annals of Applied Statistics. 9(1).
}

\examples{
##r.gls is a result object obtained from gls.simpe, gls.plink, gls.snpmat
##e.g. r.gls <- gls.simple(...);
#
# summary(r.gls);
}

\keyword{summary}
\keyword{GLS model}
