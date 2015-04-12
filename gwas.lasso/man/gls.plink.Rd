\name{gls.plink}
\alias{gls.plink}

\title{ scan PLINK dataset with longitudinal measuments by the GLS model }

\description{ 
 Scan PLINK dataset with longitudinal measument and detect the joint additive and dominant effects by the Group Lasso Model.
} 

\usage{
gls.plink(file.phe,
      file.plink.bed,file.plink.bim,file.plink.fam,
      Y.prefix,
      Z.prefix,
      covar.names,
      refit=TRUE,
      add.used=T,
      dom.used=T,
      fgwas.filter=F, 
      options=NULL) 
}
\arguments{
\item{file.phe}{Input CSV file containing longitudinal measuments and covariate data.}
\item{file.plink.bed}{Input PLINK data file, a binary file containing genotype information.}
\item{file.plink.bim}{Input PLINK data file, an extended map file.}
\item{file.plink.fam}{Input PLINK data file, family information.}
\item{Y.prefix}{String, prefix of column names of response variables.}
\item{Z.prefix}{String, prefix of measurement names.}
\item{covar.names}{String or vector of string indicating column names of covariates.}
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
An object of class \code{GLS.ret} including the following elements. 
\item{varsel_add}{Matrix, additive effects estimated by the procedure of variable selection for each SNP. 
21 columns are available, 
\itemize{
\item{(1) chromosome group}{}
\item{(2) position}{}
\item{(3) the first Legendre coefficient of additive effcts}{}
\item{(4) the secnod Legendre coefficient of additive effcts}{}
\item{(5) the third Legendre coefficient of additive effcts}{}
\item{(6) the forth Legendre coefficient of additive effcts}{}
\item{(7) quadratic sum of mean values of additive effect}{}
\item{(8) the first mean value of additive effect}{}
\item{(9) the second mean value of additive effect}{}
\item{(10) the third mean value of additive effect}{}
\item{(11) the forth mean value of additive effect}{}
\item{(12) quadratic sum of minimum value of addtive effect}{}
\item{(13) the first minimum value of addtive effect}{}
\item{(14) the second minimum value of addtive effect}{}
\item{(15) the third minimum value of addtive effect}{}
\item{(16) the forth minimum value of addtive effect}{}
\item{(17) quadratic sum of maximum value of addtive effect}{}
\item{(18) the first maximum value of addtive effect}{}
\item{(19) the second maximum value of addtive effect}{}
\item{(20) the third maximum value of addtive effect}{}
\item{(21) the forth maximum value of addtive effect}{}
}
}
\item{varsel_dom}{Matrix, dominant effects estimated by the procedure of variable selection for each SNP. 
21 columns are available, 
\itemize{
\item{(1) chromosome group}{}
\item{(2) position}{}
\item{(3) the first Legendre coefficients of dominant effcts}{}
\item{(4) the secnod Legendre coefficients of dominant effcts}{}
\item{(5) the third Legendre coefficients of dominant effcts}{}
\item{(6) the forth Legendre coefficients of dominant effcts}{}
\item{(7) quadratic sum of mean values of dominant effect}{}
\item{(8) the first mean value of dominant effect}{}
\item{(9) the second mean value of dominant effect}{}
\item{(10) the third mean value of dominant effect}{}
\item{(11) the forth mean value of dominant effect}{}
\item{(12) quadratic sum of minimum value of dominant effect}{}
\item{(13) the first minimum value of dominant effect}{}
\item{(14) the secnod minimum value of dominant effect}{}
\item{(15) the third minimum value of dominant effect}{}
\item{(16) the forth minimum value of dominant effect}{}
\item{(17) quadratic sum of maximum value of dominant effect}{}
\item{(18) the first maximum value of dominant effect}{}
\item{(19) the secnod maximum value of dominant effect}{}
\item{(20) the third maximum value of dominant effect}{}
\item{(21) the forth maximum value of dominant effect}{}
}
}
\item{varsel_cov}{Matrix, covariate effects estimated by the variable selection procedure, 19 columns as follows:
\itemize{
\item{(1) the first Legendre coefficients of covariate effcts}{}
\item{(2) the secnod Legendre coefficients of covariate effcts}{}
\item{(3) the third Legendre coefficients of covariate effcts}{}
\item{(4) the forth Legendre coefficients of covariate effcts}{}
\item{(5) quadratic sum of mean values of dominant effect}{}
\item{(6) the first mean value of covariate effect}{}
\item{(7) the second mean value of covariate effect}{}
\item{(8) the third mean value of covariate effect}{}
\item{(9) the forth mean value of covariate effect}{}
\item{(10) quadratic sum of minimum value of dominant effect}{}
\item{(11) the first minimum value of covariate effect}{}
\item{(12) the secnod minimum value of covariate effect}{}
\item{(13) the third minimum value of covariate effect}{}
\item{(14) the forthminimum value of covariate effect}{}
\item{(15) quadratic sum of maximum value of dominant effect}{}
\item{(16) the first maximum value of covariate effect}{}
\item{(17) the secnod maximum value of covariate effect}{}
\item{(18) the third maximum value of covariate effect}{}
\item{(19) the forthmaximum value of covariate effect}{}
}
}
\item{refit_add}{Matrix, additive effects estimated by the refit procedure for the SNPs selected by variable selection, 21 columns are same as the matrix of varible selection \cr (\code{varsel_add}).}
\item{refit_add}{Matrix, dominant effects estimated by the refit procedure for the SNPs selected by variable selection, 21 columns are same as the matrix of varible selection \cr (\code{varsel_dom}).}
\item{refit_add}{Matrix, covariate effects estimated by the refit procedure for the SNPs selected by variable selection, 21 columns are same as the matrix of varible selection \cr (\code{refit_cov}).}
\item{fgwas.filter}{Matrix, the p-values estimated by fGWAS method for each SNP, 5 columns as follows:  
\itemize{
\item{(1) SNP ID}{}
\item{(2) chromosome group}{}
\item{(3) position}{}
\item{(4) likelihood ratio}{}
\item{(5) p-value}{}
}
}
\item{varsel_PSRF}{Matrix, the potential scale reduction factor in the variable selection procedure. }
\item{refit_PSRF}{Matrix, the potential scale reduction factor in the refit procedure.}
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
 Li. J., Wang, Z., Li, R., Wu, R.(2015). Bayesian Group Lasso for nonparametric varying-coefficient models with application to functional genome-wide association studies. The Annals of Applied Statistics. 9(1).
}

\examples{
if(0)
{
gls.plink( file.phe.long, file.plink.bed, file.plink.bim, file.plink.fam, 
Y.prefix="Y", Z.prefix="Z", covar.names=c("X"), fgwas.filter = TRUE );
}
}

\keyword{functions}
\keyword{GLS model}
\keyword{PLINK}