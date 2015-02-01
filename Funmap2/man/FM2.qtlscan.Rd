\name{FM2.qtlscan}
\alias{FM2.qtlscan}

\title{Main function to perform all analyses}

\description{ 
  Perform all analyses for the real data as a main function.
} 

\usage{ FM2.qtlscan(phenp_file, geno_file, marker_file, model_type, cross_type, covar_type=COVAR_AR1, options=list()) }

\arguments{
\item{pheno_file}{a phenotype file. The format is described in the details section of \code{\link{FM2.load_data}}.}
\item{geno_file}{a genotype file. The format is described in the details section of \code{\link{FM2.load_data}}.}
\item{marker_file}{a marker file. The format is described in the details section of \code{\link{FM2.load_data}}.}
\item{model_type}{CURVE_LC, CURVE_BI, CURVE_PC, CURVE_NP, CURVE_EXP. see \code{\link{FM2.param}}.}
\item{cross_type}{CROSS_BC, CROSS_F2 or CROSS_RIL. see \code{\link{FM2.param}}.}
\item{options}{optional values for system configuration, see the following details.} 
}

\value{
A result object of hypothesis tests is returned. Although three objects besides the result object are generated, only then result object can be returned as a return value. Other objects is saved into a RDATA file which its name can be specified by \code{options} argument.

}

\details{


\strong{FM2.qtlscan} is main function of the Funmap package, it encapsulates a consecutive procedures, including:

1) Loading the phenotype, genotype and marker file.\cr
2) Performing the hypothesis tests.\cr
3) Storing the intermediate results into RDATA(.rdata) file before the permutation.\cr
4) Outputing a PDF report which includes the summary of the data and hypothesis test 10.\cr
5) Taking a long time to execute permutation parallelly or not.\cr
6) Storing all objects(raw data, hypothesis result, permutation result) into RDATA(.rdata) file.\cr
7) Summarize all objects.\cr

Three objects generated in this function can be retrieved by the following variables.

1) Raw data object describled in \code{FM2.load_data}: \strong{dat_os};\cr
2) Hypothesis result object describled in \code{FM2.qtlscan}: \strong{ret_os};\cr
3) Permutation result object describled in \code{FM2.permutation}: \strong{ret_perm_os};\cr

If RDATA(.rdata) filename is not specified in the argument of options, the naming follows the phenotype filename.
If the phenotype filename is "PHENO_ABC.csv", the .rdata file would be "PHENO_ABC.csv.rdata".
RDATA(.rdata) file can be loaded by the \strong{load} command in R, so the user can make further analysis based on these objects.

In the \code{\link{FM2.set_value}} page, all system variables are listed. But for \code{FM2.one_step} function, the following items are important.


1) \code{cluster_count}, default=1, indicates whether the clusters (parallel computation) are available.\cr
2) \code{permu_loop}, default=1000, specifies the permutation count.\cr
3) \code{file.summary}, default=NA, a filename where the summary information which displayed in the console will be saved.\cr
4) \code{file.rdata}, default=NA, a filename for RDATA format file where the data object, result object and permutation object will be saved.\cr
5) \code{file.report}, default=NA, a filename for PDF report file where the summaries and figures will be saved.\cr

}

\author{Zhong Wang, \email{zzw2@psu.edu} }

\examples{

library(Funmap);
ret <- FM2.qtlscan( "pheno.csv", "geno.csv", "marker.csv", CURVE_LC, CROSS_BC );

}

\keyword{functions}
