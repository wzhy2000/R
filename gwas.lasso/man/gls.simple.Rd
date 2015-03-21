\name{gls.simple}
\alias{gls.simple}

\title{ The Glasso model }

\description{ 
  A function used for identifying significant SNPs that control phenotypic variation and estimating additive and dominant genetic effects of SNPs.
} 

\usage{ gls.simple(file.phe, file.snp, Y.prefix, Z.prefix, covar.names, refit=TRUE, 
add.used=T, dom.used=T, fgwas.filter=F, options=NULL )  }

\arguments{
\item{file.phe}{Input phenotypic data file representing longitudinal phenotypic data and covariate data.}
\item{file.snp}{Input genotypic data file representing chromosome, position and genotype of SNP.}
\item{Y.prefix}{String, the prefix of the column names for response variables}
\item{Z.prefix}{String, the prefix of the measurement names}
\item{covar.names}{String or vector with strings representing the names of covariates }
\item{refit}{ Logical, if TRUE, refit should be processed. Default: T.}
\item{add.used}{Logical, if TRUE, the program will estimate additive effects of SNPs. Default: T.}
\item{dom.used}{Logical, if TRUE, the program will estimate dominant effects of SNPs. Default: T.}
\item{fgwas.filter}{Logical, If TRUE, the filtering function will be performed. Default: F.}
\item{options}{A list that contains nParallel.cpu, nPiecewise.ratio, nMcmcIter, fBurnInRound, fRhoTuning, fQval.add, fQval.dom, fgwas.cutoff, and debug.\cr
 Defaults to "0", "2", "2000", "0.3", "0.095", "0.05", "0.05", "0.05", "F", respectively}
}

\value{
Returns a list with multiple values: varsel_add, varsel_dom, varsel_cov, refit_add, refit_dom, refit_cov and fgwas.filter. varsel_add and varsel_dom represent additive and dominant effects of SNPs selected by the variable selection method, and varsel_cov represents covariate effects in the variable selection procedure. refit _add and refit _dom represent additive and dominant effects of SNPs selected by the refit method, and refit _cov represents covariate effects in the refit procedure.
}

\details{

A two-stage procedure for multi-SNP modeling and analysis in GWASs is proposed, by first filtering significant SNPs if the number of SNPs is greatly large, and then formulating Bayesian Group lasso to select a subset of significant SNPs by selecting Legendre coefficients that are not identically zero. The unknown parameters could be estimated by using the MCMC algorithm.\cr
Then, to ameliorate the bias of the parameter estimates introduced by group lasso penalties, it is necessary to refit the fGWAS model after variable selection, where only selected SNPs are included in the final model. The parameter refit indicates whether to carry out the refit function, and the default is TRUE. If the refit is TRUE, the program will refit the result, otherwise, it will not.\cr
For the function gls.simple, the input data are in a simple format. The phenotype file should be filled with numerical values and missing values recorded as NA. The package requires individual ids in the phenotype file and the genotype file to be consistent. In this file one row represents one subject. The file format represents the phenotypic data of gls model where X_* fields mean covariates data, Z_* fields mean measurements and Y_* fields mean longitudinal phenotypic data.\cr
The parameter named "option" provides a list. In the list, nParallel.cpu means the number of the CPU used to compute, nPiecewise.ratio means the ratio of grouping the SNPs, and the number of SNPs of each group is less than or equal to the number of the individuals times nPiecewise.ratio means the number of the CPU used to compute. nMcmcIter means the number of times of the MCMC algorithm execution. fBurnInRound means a ratio to cut off unstable calculations in order to make the result precise. fRhoTuning means a ratio used for narrowing the scope of parameter pho to reduce parameter oscillation. fQval.add and fQval.dom means the Q-value of additive and dominant genetic effects respectively. fgwas.cutoff means a threshold value, and debug means a logical value, indicating whether or not debug should be outputted.\cr

}

\examples{
library(gwas.lasso)
phe.out <- "gls.test.simple.phe";
snp.out <- "gls.test.simple.snp";
gls.simulate(phe.out, snp.out)
if(0) ret <- gls.simple(phe.out, snp.out, 
	Y.prefix="Y", Z.prefix="Z", covar.names=c("X_1","X_2"), fgwas.filter = TRUE );	
}

\keyword{functions}

