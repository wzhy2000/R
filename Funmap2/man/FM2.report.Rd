\name{FM2.report}
\alias{FM2.report}

\title{PDF report for data and result object}

\description{ 
Output a PDF report including the summary of data and result object.
} 

\usage{ FM2.report( pdf.file, dat_obj, res_obj, options=list()  ) }

\arguments{
\item{ pdf.file}{PDF file name.}
\item{dat_obj}{a data object returned by \code{\link{FM2.simulate}} or \code{\link{FM2.load_data}} }
\item{res_obj}{a result object returned by \code{\link{FM2.qtlscan}} }
}

\value{
A PDF file name is returned if PDF report is generated successfully.
}

\details{

This function don't use the HaruPDF package anymore!!!
It outputs the summary information and figures into a PDF file for the data object and the hypothesis test 10 of the result object.

The following link is an example for \code{FM2.report}.

\url{http://statgen.psu.edu/software/funmap/report_demo.pdf}

}

\author{Zhong Wang, \email{zzw2@psu.edu} }

\examples{
par <- FM2.param(par_LC, CURVE_LC, CROSS_BC);
dat <- FM2.simulate(par);
res <- FM2.qtlscan(dat);
FM2.report("test.pdf", dat, ret);
}

\keyword{functions}
