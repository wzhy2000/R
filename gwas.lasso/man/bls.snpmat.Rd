\name{bls.snpmat}
\alias{bls.snpmat}

\title{ The Blasso model }

\description{ 
  A function used for identifying significant SNPs that control phenotypic variation and estimating additive and dominant genetic effects of SNPs.
} 

\usage{ bls.snpmat (phe.mat, snp.mat, Y.name, covar.names, refit=TRUE, 
add.used=T, dom.used=T, fgwas.filter=F, options=NULL)  }

\arguments{
\item{phe.mat}{Matrix, phenotypic data with single point in time and covariate data.}
\item{snp.mat}{A matrix representing chromosome, position and genotype of SNPs.}
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

Input data are in a matrix format.The phenotype file should be filled with numerical values and missing values recorded as NA. In this file one row represents one subject. The file format represents the phenotypic data of Bls model where X_* fields mean covariates data and Y field means phenotypic data from a single time point.\cr
Each line of the Genotype file describes a single SNP and must contain exactly 2 columns: chromosome and position of SNP. Three genotypes (aa=0,Aa=1,AA=2) and missing data (coded as -1) are valid SNP values.\cr
}

\examples{
library(gwas.lasso)
phe.out <- "bls.test.simple.phe"  
snp.out <- "bls.test.simple.snp"
bls.simulate(phe.out, snp.out)
tb.phe<-read.csv(phe.out);
tb.snp<-read.csv(snp.out);
if(0)
ret <- bls.snpmat (tb.phe, tb.snp, Y.name="Y", covar.names=c("X_1","X_2"), 
		   fgwas.filter = FALSE );
}

\keyword{functions}