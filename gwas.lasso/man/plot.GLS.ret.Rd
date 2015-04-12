\name{plot.GLS.ret}
\alias{plot.GLS.ret}

\title{Plotting the detection result obtained by the GLS model}

\description{ 
  Plot a result object obtained by the GLS method.
} 

\usage{ 
\method{plot}{GLS.ret}(x, y=NULL, ..., fig.prefix=NULL) 
}

\arguments{
\item{x}{An object obtained by \code{\link{gls.simple}} or \code{\link{gls.plink}}or \code{\link{gls.snpmat}}}
\item{y}{An null parameter}
\item{fig.prefix}{String, the prefix of output file name}
\item{...}{Other parameters}
}

\details{
The function plot can output three types of PDF figure, including: \cr

1) The Manhattan figure exported by the fGWAS method (*.fgwas.pdf) \cr

2) The genetic effects of all SNPs estimated by the variable select procedure (*. varsel.pdf.) \cr

3) The genetic effects of significant SNPs estimated by the refit procedure (*.refit.pdf). \cr

The Manhattan figure gives -log10 (p-values) for each SNP, from which the SNP with -log10 p-values greater than 
the threshold value specified in the control parameters will be selected to variable selection. \cr

In the GLS model, the figures of genetic effects will output the time-varying additive and dominant effects for each significant SNP.
}

\author{Zhong Wang and Nating Wang }

\references{
 Li. J., Wang, Z., Li, R., Wu, R.(2015). Bayesian Group Lasso for nonparametric varying-coefficient models with application to functional genome-wide association studies. The Annals of Applied Statistics. 9(1).
}

\seealso{See \code{\link{gls.simple}} or \code{\link{gls.plink}}or \code{\link{gls.snpmat}} to generate an object }

\examples{
##r.gls is a result object obtained from gls.simpe, gls.plink, gls.snpmat
##e.g. r.gls <- gls.simple(...);
#
# plot(r.gls, fig.prefix="r.gls.pdf");
}

\keyword{plot}
\keyword{GLS model}

