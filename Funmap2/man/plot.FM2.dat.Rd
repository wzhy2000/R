\name{plot.FM2.dat}
\alias{plot.FM2.dat(dat_obj)}
\title{Figures of a data object}

\description{ 
  Draw figures for a data object.
} 

\usage{ plot( dat_obj, plot_type=NA) }

\arguments{
\item{dat_obj}{A data object returned by \code{\link{FM2.simulate}} or\code{\link{FM2.load_data}}   .}
\item{plot_type}{The plot type, 1 is for tiled curves and 2 is for overlapping curves.}
}

\details{


Two figures can be outputted to R console.

1) tiled curves for every individuals.

2) overlapping curves for every individuals

These figures can also be exported by the summary command. The summary command will generate this figure to a PDF file, however the plot command will output these figure to the R console.

An example of this command is available in the following URL.

\url{http://statgen.psu.edu/software/funmap/plot.data1.jpg}.
\url{http://statgen.psu.edu/software/funmap/plot.data2.jpg}.

}


\author{Zhong Wang, \email{zzw2@psu.edu} }

\examples{
par <- FM2.param(par_LC, CURVE_LC, CROSS_BC);
dat <- FM2.simulate(par);
plot(dat);
}

\keyword{functions}
