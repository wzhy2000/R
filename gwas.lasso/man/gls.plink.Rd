\name{gls.plink}
\alias{gls.plink}

\title{ The Glasso model }

\description{ 
  A function used for identifying significant SNPs that control phenotypic variation and estimating additive and dominant genetic effects of SNPs.
} 

\usage{ gls.plink (file.phe, file.plink.bed, file.plink.bim, file.plink.fam, Y.prefix, Z.prefix,
covar.names, refit=TRUE, add.used=T, dom.used=T, fgwas.filter=F, options=NULL ) }

\arguments{
\item{file.phe}{Input phenotypic data file representing longitudinal phenotypic data and covariate data.}
\item{file.plink.bed}{Input file, a binary file containing genotype information}
\item{file.plink.bim}{Input file, an extended MAP file}
\item{file.plink.fam}{Input file stored phenotype information}
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

The phenotype file should be filled with numerical values and missing values recorded as NA. The package requires individual ids in the phenotype file and the genotype file to be consistent. In this file one row represents one subject. The file format represents the phenotypic data of gls model where X_* fields mean covariates data, Z_* fields mean measurements and Y_* fields mean longitudinal phenotypic data.\cr
The genotype files are in a plink format.The *.bed file is a compressed file and when you try to view it you will only see lots of strange characters on the screen. \cr
The phenotype information is stored in the *.fam file where one row represent one subject and the first six columns are mandatory: Family ID, Individual ID, Paternal ID, Maternal ID, Sex (1=male; 2=female; other=unknown) and Phenotype. \cr
The the *.bim file is an extended MAP file where each line of file describes a single individual and it must contain exactly 4 columns: chromosome, snp identifier, genetic distance and base-pair position (bp units). The following two extra cols are allele names. Genetic distance can be specified in centimorgans with the --cm flag. The .bim file can be viewed with a standard text editor as well.The *.fam file and *.bim file are plain text files that can be viewed with a standard text editor.\cr
}

\examples{
if(0)
{

	library(gwas.lasso)
	gls.plink( file.phe.long, file.plink.bed, file.plink.bim, file.plink.fam, 
	Y.prefix="Y", Z.prefix="Z", covar.names=c("X"), fgwas.filter = TRUE );
}
}

\keyword{functions}