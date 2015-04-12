\name{bls.plink}
\alias{bls.plink}

\title{ scan PLINK dataset with single measument by the BLS model}

\description{ 
Scan PLINK dataset with single measument and detect the joint additive and dominant effects by the Bayesian Lasso Model.
} 

\usage{
bls.plink(file.phe, 
      file.plink.bed, file.plink.bim, file.plink.fam,
      Y.name, 
      covar.names, 
      refit=TRUE, 
      add.used=T, 
      dom.used=T, 
      fgwas.filter=F, 
      options=NULL)}

\arguments{
\item{file.phe}{Input CSV file containing one single measument and covariate data.}
\item{file.plink.bed}{Input PLINK data file, a binary file containing genotype information.}
\item{file.plink.bim}{Input PLINK data file, an extended map file.}
\item{file.plink.fam}{Input PLINK data file, family information.}
\item{Y.name}{String, column name of response variables.}
\item{covar.names}{String or vector of string indicating column names of covariates. }
\item{refit}{ Logical, if TRUE, the refit procedure should be processed. Default: TRUE.}
\item{add.used}{Logical, if TRUE, additive effects will be estimated jointly. \cr Default: TRUE.}
\item{dom.used}{Logical, if TRUE, dominant effects will be estimated jointly. \cr Default: TRUE.}
\item{fgwas.filter}{Logical, If TRUE, the SNPs will be filtered by fGWAS midel before the LASSO methods are applied. Default: FALSE.}
\item{options}{A list containing control parameters, including \code{nParallel.cpu}, \code{nMcmcIter}, \code{nPiecewise.ratio}, \code{fBurnInRound}, \code{fRhoTuning}, \code{fQval.add},\cr
\code{fQval.dom}, \code{fgwas.cutoff}, and \code{debug}.\cr
see details.
}
}

\value{
An object of class \code{BLS.ret} including the following elements. 
\item{varsel}{Matrix, additive and dominant effects estimated by the procedure of variable selection for each SNP. 
11 columns are available, 
\itemize{
\item{(1) chromosome group}{}
\item{(2) position}{}
\item{(3) significant sign of additive effect}{}
\item{(4) mean value of additive effect}{}
\item{(5) minimum value of addtive effect}{}
\item{(6) maximum value of addtive effect}{}
\item{(7) significant sign of dominant effect}{}
\item{(8) mean value of dominant effect}{}
\item{(9) minimum value of dominant effect}{}
\item{(10) maximum value of dominant effect}{}
\item{(11) Heritability rate}{}
}
}
\item{varsel_cov}{Matrix, covariate effects estimated by the variable selection procedure, 4 columns as follows:
\itemize{
\item{(1) significant sign of covariate effect}{}
\item{(2) overall mean value of covariate effect}{}
\item{(3) minimum value of covariate effect}{}
\item{(4) maximum value of covariate effect}{}
}
}
\item{refit}{Matrix, additive and dominant effects estimated by the refit procedure for the SNPs selected by variable selection, 11 columns are same as the matrix of varible selection(\code{varsel}).}
\item{refit_cov}{Matrix, covariate effects destimated by the refit procedure, 4 columns are same as the matrix of varible selection(\code{varsel_cov}).}
\item{fgwas.filter}{Matrix, the p-values estimated by fGWAS method for each SNP, 5 columns as follows:  
\itemize{
\item{(1) SNP ID}{}
\item{(2) chromosome group}{}
\item{(3) position}{}
\item{(4) likelihood ratio}{}
\item{(5) p-value}{}
}
}
\item{options}{ List, data file names and all parameters, including the control parameters defined in \code{options} }
}

\details{
(1) PLINK dataset\cr
The PLINK dataset include three files, including \emph{bed},\emph{bim} and \emph{fam} file.\cr
The \emph{bed} file is a compressed binary file containing genotype information.\cr
The \emph{bim} file is an extended map file where each line of this file defines a SNP information 
containing exactly 4 columns: chromosome, SNP identifier, genetic distance and base-pair 
position (bp units). \cr
The \emph{fam} file is a family file containing the individuals and phenotypic information including 
six columns: Family ID, Individual ID, Paternal ID, Maternal ID, Sex (1=male; 2=female; 
other=unknown) and phenotype.\cr
This package does \strong{not} use the phenotype information stored in the fam file. The Phenotype 
values and covariates are stored into a CSV file which should be filled with numerical values and 
missing values encoded as NA.\cr

(2) fGWAS model\cr
A SNP filter strategy named fGWAS is proposed in this package. This method is used to detect the
marginal genetic effects. If the parameter \code{fgwas.filter} is true, the fGWAS method will be
employed to remove the non-significant SNPs. Once SNPs are filtered by the fGWAS method, significant 
SNPs are remained to do LASSO process.\cr

(3) Varible selection and refit\cr
A two-stage procedure based on Bayesian Lasso or Group Lasso is proposed, including variable selection 
and refit process. It is necessary to refit the fGWAS model after variable selection, where only selected SNPs 
are included in the final process. The parameter \code{refit} indicates whether to carry out the refit function.\cr

(4) Control parameters\cr
The control parameters include the following items.\cr
\itemize{
\item{\code{nParallel.cpu}}    {Default:0, the number of the CPU used to compute.}
\item{\code{nPiecewise.ratio}} {Default:2, the ratio of grouping SNPs.}
\item{\code{nMcmcIter}}        {Default:2000, the number of times of the MCMC algorithm execution. }
\item{\code{fBurnInRound}}     {Default:0.3, a ratio to cut off unstable calculations in order to make the result precise.}
\item{\code{fRhoTuning}}       {Default:0.095, a ratio used for narrowing the scope of parameter rho to reduce parameter oscillation.}
\item{\code{fQval.add}}        {Default:0.05, }
\item{\code{fQval.dom}}        {Default:0.05,the Q-value of additive and dominant genetic effects respectively.}
\item{\code{fgwas.cutoff}}     {Default:FALSE, a threshold value, and \bold{debug} means a logical value, indicating whether or not debug should be outputted.}
}

}

\author{Zhong Wang and Nating Wang }

\references{
Li, J., Das, K., Fu, G., Li, R., & Wu, R. (2011). The Bayesian lasso for genome-wide association studies. Bioinformatics, 27(4), 516-523.
}

\examples{
if(0)
{
bls.plink( file.phe.long, file.plink.bed, file.plink.bim, file.plink.fam, 
	Y.name="Y", covar.names="sex", fgwas.filter = T, refit = F);
}	
}

\keyword{functions}
\keyword{BLS model}
\keyword{PLINK}

