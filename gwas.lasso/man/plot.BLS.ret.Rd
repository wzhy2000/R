\name{plot.BLS.ret}
\alias{plot.BLS.ret}

\title{Plotting the detection result obtained by the BLS model}

\description{ 
  Plot a result object obtained by the BLS method.
} 

\usage{ 
\method{plot}{BLS.ret}(x, y=NULL, ..., fig.prefix=NULL) 
}

\arguments{
\item{x}{An object obtained by \code{\link{bls.simple}} or \code{\link{bls.plink}}or \code{\link{bls.snpmat}}.}
\item{y}{An null parameter.}
\item{fig.prefix}{String, the prefix of output file name.}
\item{...}{Other parameters.}
}

\details{
The function plot can output three types of PDF figure, including: \cr

1) The Manhattan figure exported by the fGWAS method (*.fgwas.pdf) \cr

2) The genetic effects of all SNPs estimated by the variable select procedure (*. varsel.pdf.) \cr

3) The genetic effects of significant SNPs estimated by the refit procedure (*.refit.pdf). \cr

The Manhattan figure gives -log10 (p-values) for each SNP, from which the SNP with -log10 p-values 
greater than the threshold value specified in the control parameters will be selected to variable selection.\cr

In the BLS model, the figures of genetic effects will output heritability sub-graph.\cr
}

\author{Zhong Wang and Nating Wang }

\references{
 Li, J., Das, K., Fu, G., Li, R., & Wu, R. (2011). The Bayesian lasso for genome-wide association studies. Bioinformatics, 27(4), 516-523.
}

\seealso{See \code{\link{bls.simple}} or \code{\link{bls.plink}}or \code{\link{bls.snpmat}} to generate an object. }

\examples{

##r.bls is a result object obtained from bls.simpe, bls.plink, bls.snpmat
##e.g. r.bls <- bls.simple(...);
#
# plot(r.bls, fig.prefix="r.bls.pdf");

}

\keyword{plot}
\keyword{BLS model}

