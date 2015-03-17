\name{plot.Gls.ret}
\alias{plot.Gls.ret}

\title{Plotting the detection result}

\description{ 
  Plot a result object exported from Gls method.
} 

\usage{ plot(r.gls, refit=T) }

\arguments{
\item{r.gls}{an object return by \code{\link{gls.simple}} or \code{\link{gls.plink}}or \code{\link{gls.snpmat}}}
\item{refit}{Logical, if TRUE, the refit result will be exported.}
}

\value{
The result will be saved as PDF format files, which contain three types: *.fgwas.pdf, *.refit.pdf and *. varsel.pdf.\cr
*.fgwas.pdf holds the figure which gives ¨Clog10 p-values for each SNP, from which the SNP with ¨Clog10 p-values greater than the threshold value set before will be selected.\cr
*. varsel.pdf holds the figures including addictive effects and dominant effects of selected SNPs and their heritability. Using the Legendre coefficient estimates, we plot their time-varying additive effects and dominant effects in the figures. Utilizing the poster samples of a, d, and the observed genotypes, we can calculate the proportion of the phenotypic variance explained by a particular SNP heritability. Since heritability estimates are unitless, they could guide variable selection and identify SNPs that have relatively large effects on the phenotype. \cr
*.refit.pdf holds the figures including addictive effects and dominant effects of selected significant SNPs and their heritability. Formulating Bayesian Group Lasso, we select a subset of significant SNPs.\cr
}

\details{
The plot command can export figures based on the result object, which contains multiple values.  

}

\author{Zhong Wang, \email{zzw2@psu.edu} }

\examples{
plot(r.gls);
}

\keyword{plot}