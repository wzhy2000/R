\name{bls.plink.tped}
\alias{bls.plink.tped}


\title{ The Blasso model }

\description{ 
   A function used for identifying significant SNPs that control phenotypic variation and estimating additive and dominant genetic effects of SNPs.
} 

\usage{ bls.plink.tped( file.phe, file.plink.tped, file.plink.tfam, Y.name, 
covar.names, refit=TRUE, add.used=T, dom.used=T, options=NULL)}

\arguments{
\item{file.phe}{Input phenotypic data file representing phenotypic data with single point in time and covariate data.}
\item{file.plink.tped}{Input file, a file listing genotypes}
\item{file.plink.tfam}{Input file, a file containing genotype information}
\item{Y.name}{String, column name for the response variables}
\item{covar.names}{String or vector with strings representing the names of covariates }
\item{refit}{ Logical, if TRUE, refit should be processed. Default: T.}
\item{add.used}{Logical, if TRUE, the program will estimate additive effects of SNPs. Default: T.}
\item{dom.used}{Logical, if TRUE, the program will estimate dominant effects of SNPs. Default: T.}
\item{options}{A list that contains nParallel.cpu, nPiecewise.ratio, nMcmcIter, fBurnInRound, fRhoTuning, fQval.add, fQval.dom, fgwas.cutoff, and debug.\cr
 Defaults to "0", "2", "2000", "0.3", "0.095", "0.05", "0.05", "0.05", "F", respectively}
}

\value{
Returns a list with multiple values: varsel_add, varsel_dom, varsel_cov, refit_add, refit_dom, refit_cov and fgwas.filter. varsel_add and varsel_dom represent additive and dominant effects of SNPs selected by the variable selection method, and varsel_cov represents covariate effects in the variable selection procedure. refit _add and refit _dom represent additive and dominant effects of SNPs selected by the refit method, and refit _cov represents covariate effects in the refit procedure.
}

\details{

The phenotype file should be filled with numerical values and missing values recorded as NA. The package requires individual ids in the phenotype file and the genotype file to be consistent. In this file one row represents one subject. The file format represents the phenotypic data of Bls model where X_* fields mean covariates data and Y field means phenotypic data from a single time point.\cr
Input genotypic files contains two text files: one (TPED) containing SNP and genotype information where one row is a SNP; one (TFAM) containing individual and family information, where one row is an individual.\cr
The first 4 columns of a TPED file are chromosome, rs# or snp identifier, Genetic distance (morgans) and Base-pair position (bp units).Then all genotypes are listed for all individuals for each particular SNP on each line. The TFAM file is just the first six columns of a standard PED file. In otherwords, we have just taken the standard PED/MAP file format, but swapped all the genotype information between files, after rotating it 90 degrees.\cr
This kind of format can be convenient to work with when there are very many more SNPs than individuals.
}

\author{Zhong Wang, \email{zzw2@psu.edu} }


\keyword{functions}



