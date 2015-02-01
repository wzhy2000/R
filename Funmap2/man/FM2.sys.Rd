\name{System}
\alias{FM2.set_value}
\alias{FM2.get_value}
\title{System environment configuration}

\description{ 
	Configure system variables for Funmap Package.
} 

\usage{ 
FM2.set_value( key, value ); 
FM2.get_value( key, def=NA );
}

\arguments{
\item{key}{the count of clusters.}
\item{value}{the permutation count.}
\item{def}{the default value.}
}

\details{

The Funmap employs a system object to configure the system variables. This system object is created after \code{\link{FM2.param}} or \code{\link{FM2.load_data}}, \code{\link{FM2.simu_test}} or \code{\link{FM2.one_step}} is called. 
The user can use \code{FM2.set_value()} and \code{FM2.get_value()} to manipulate these system variables.

In the input arguments of \code{\link{FM2.one_step}} and \code{\link{FM2.qtlscan}}, the system variables can be adjusted via the argumant of \code{options}. The following section lists all system variables.

1) \strong{cluster_count}, default=1, the cluster count for parallel permutation.\cr
2) \strong{permu_loop}, default=1000, the count of permutation loop.\cr
3) \strong{file.summary}, default=NA, a filename where the summary information which displayed in the console will be saved.\cr
4) \strong{file.rdata}, default=NA, a filename for RDATA format file where the data object, result object and permutation object will be saved.\cr
5) \strong{file.report}, default=NA, a filename for PDF report file where the summaries and figures will be saved.\cr
6) \strong{scan_step}, default=2, an interval distance used to scan flanking marker, default is 2cm.\cr
7) \strong{peak_count}, default=5, a number shows how many significant QTL will be selected.\cr
8) \strong{plot_doctype}, default=pdf, the figure output type for summary command.\cr
9) \strong{np.order}, default=6, the order of Legendre polynomial for nonparametric method.\cr

For the items of \code{file.summary}, \code{file.rdata} and \code{file.report}, if they are undefined, the Funmap will give filenames based on the phenotype
filename.

}


\author{Zhong Wang, \email{zzw2@psu.edu} }

\examples{
FM2.set_value("cluster_count", 10 ); 
FM2.set_value("permu_loop", 2000 );
}

\keyword{System}


