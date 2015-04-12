\name{gls.simple}
\alias{gls.simple}

\title{ scan Simple dataset with single measument by the GLS model}

\description{ 
Scan Simple dataset with longitudinal measuments and detect the joint additive and dominant effects by the Group Lasso Model.
} 

\usage{
gls.simple(file.phe, 
      file.snp, 
      Y.prefix, 
      Z.prefix, 
      covar.names,
      refit=TRUE, 
      add.used=T, 
      dom.used=T, 
      fgwas.filter=F, 
      options=NULL)}

\arguments{
\item{file.phe}{Input CSV file containing longitudinal measuments and covariate data.}
\item{file.snp}{Input CSV file containing chromosome, position and genotypes of SNPs.}
\item{Y.prefix}{String, prefix of column names of response variables.}
\item{Z.prefix}{String, prefix of measurement names.}
\item{covar.names}{String or vector of string indicating column names of covariates. }
\item{refit}{Logical, if TRUE, the refit procedure should be processed. Default: TRUE.}
\item{add.used}{Logical, if TRUE, additive effects will be estimated jointly. \cr Default: TRUE.}
\item{dom.used}{Logical, if TRUE, dominant effects will be estimated jointly. \cr Default: TRUE.}
\item{fgwas.filter}{Logical, If TRUE, the SNPs will be filtered by fGWAS midel before the LASSO methods are applied. Default: FALSE.}
\item{options}{A list containing control parameters, including \code{nParallel.cpu}, \code{nMcmcIter}, \code{nPiecewise.ratio}, \code{fBurnInRound}, \code{fRhoTuning}, \code{fQval.add},\cr
\code{fQval.dom}, \code{fgwas.cutoff}, and \code{debug}.\cr
see details.}
}

\value{
An object of class \code{GLS.ret}, see "Details" in the \code{\link{gls.plink}}.
}

\details{
(1) simple dataset\cr
The simple format is a user-defined format which is designed to store smallamount of SNPs 
for users who do not use PLINK.\cr
The genotypic data are stored in the CSV format, where each line describes a single SNP and must 
start with 2 columns of chromosome information (chromosomenumber and SNP position).Three genotypes 
(aa=0,Aa=1,AA=2) and missing data (coded as -1 or NA ) are valid SNP values.\cr

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
phe.out <- "gls.test.simple.phe"
snp.out <- "gls.test.simple.snp"
gls.simulate(phe.out, snp.out);
ret <- gls.simple(phe.out, snp.out, Y.prefix="Y", Z.prefix="Z", 
covar.names=c("X_1","X_2"), fgwas.filter = FALSE );	
}

\keyword{functions}
\keyword{GLS model}
\keyword{simple}

