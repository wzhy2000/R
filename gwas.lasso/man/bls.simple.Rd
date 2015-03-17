\name{bls.simple}
\alias{bls.simple}

\title{ The Blasso model }

\description{ 
   A function used for identifying significant SNPs that control phenotypic variation and estimating additive and dominant genetic effects of SNPs.
} 

\usage{ bls.simple (file.phe, file.snp, Y.name, covar.names, refit=TRUE, add.used=T, dom.used=T, fgwas.filter=F, options=NULL) }

\arguments{
\item{file.phe}{Input phenotypic data file representing phenotypic data with single point in time and covariate data.}
\item{file.snp}{Input genotypic data file representing chromosome, position and genotype of SNPs.}
\item{Y.name}{String, the column name for the response variables}
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

A two-stage procedure for multi-SNP modeling and analysis in GWASs is proposed, by first filtering the less significant SNPs if the number of SNPs is greatly large, and then formulating Bayesian lasso to select a subset of significant SNPs by selecting functional coefficients that are not identically zero. The unknown parameters are estimated by using the MCMC algorithm. \cr
The parameter fgwas.filter indicates whether to filter the original SNPs information, and the default is TRUE. If the fgwas.filter is TRUE, the filtering function will be performed, otherwise, it will not be. \cr
Then, to ameliorate the bias of the parameter estimates introduced by group lasso penalties, it is necessary to refit the fGWAS model after variable selection, where only selected SNPs are included in the final model. The parameter refit indicates whether to carry out the refit function, and the default is TRUE. If the refit is TRUE, the program will refit the result, otherwise, it will not be.\cr
Input data are in a simple format. The phenotype file should be filled with numerical values and missing values recorded as NA. The package requires individual ids in the phenotype file and the genotype file to be consistent. In this file one row represents one subject. The file format represents the phenotypic data of Bls model where X_* fields mean covariates data and Y field means phenotypic data from a single time point.\cr
option provides a list. In the list, nParallel.cpu means the number of the CPU used to compute, nPiecewise.ratio means the ratio of grouping the SNPs, and the number of SNPs of each group is less than or equal to the number of the individuals times nPiecewise.ratio means the number of the CPU used to compute. nMcmcIter means the number of times of the MCMC algorithm execution. fBurnInRound means a ratio to cut off unstable calculations in order to make the result precise. fRhoTuning means a ratio used for narrowing the scope of parameter ¦Ñ to reduce parameter oscillation. fQval.add and fQval.dom means the Q-value of additive and dominant genetic effects respectively. fgwas.cutoff means a threshold value, and debug means a logical value, indicating whether or not debug should be outputted.\cr

}

\author{Zhong Wang, \email{zzw2@psu.edu} }

\examples{
library(glasso)
bls.simulate(phe.out="bls.simple.phe", snp.out="bls.simple.snp");
ret <- bls.simple(phe.out="bls.simple.phe", snp.out="bls.simple.snp", Y.name="Y", covar.names=c("X_1","X_2"), 
refit=T, add.used=T, dom.used=F, fgwas.filter=T );
}

\keyword{functions}