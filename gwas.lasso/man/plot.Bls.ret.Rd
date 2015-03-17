\name{plot.Bls.ret}
\alias{plot.Bls.ret}

\title{Plotting the detection result}

\description{ 
  Plot a result object exported from Bls method.
} 

\usage{ plot(r.bls, Brefit=T) }

\arguments{
\item{r.gls}{an object return by \code{\link{bls.simple}} or \code{\link{bls.plink}}or \code{\link{bls.snpmat}}}
\item{refit}{Logical, if TRUE, the refit result will be exported.}
}

\value{
The result will be saved as PDF format files, which contain three types: *.fgwas.pdf, *.refit.pdf and *. varsel.pdf.\cr
*.fgwas.pdf holds the figure which gives ¨Clog10 p-values for each SNP, from which the SNP with ¨Clog10 p-values greater than the threshold value set before will be selected.\cr
*. varsel.pdf holds the figures including addictive effects and dominant effects of selected SNPs and their heritability. Their time-varying additive effects and dominant effects are plotted in the figures using the Legendre coefficient estimates. It is possible to calculate the proportion of the phenotypic variance explained by a particular SNP heritability. Since heritability estimates are unitless, they could guide variable selection and identify SNPs that have relatively large effects on the phenotype. \cr
*.refit.pdf holds the figures including addictive effects and dominant effects of selected significant SNPs and their heritability.\cr
}

\details{
The plot command can export figures based on the result object, which contains multiple values.  

}

\author{Zhong Wang, \email{zzw2@psu.edu} }

\examples{
plot(r.bls);
}

\keyword{plot}