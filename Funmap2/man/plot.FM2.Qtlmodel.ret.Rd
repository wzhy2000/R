\name{plot.FM2.Qtlmodel.ret}
\alias{plot.FM2.Qtlmodel.ret( res_obj, dat_obj )}

\title{ Figures of qtlscan tests}

\description{ 
  Plot some figures based on the results of hypothesis test.
} 

\usage{plot(res_obj, dat_obj ) }

\arguments{
\item{res_obj}{a result object of hypothesis tests returned by \code{\link{FM2.qtlscan}}.}
\item{dat_obj}{a data object  returned by \code{\link{FM2.simulate}} or \code{\link{FM2.load_data}}.}
}

\details{

For logistic curve, the result of hypothesis test 10 and 15 can draw some figures. For other two curves, becuase only hypothesis test 10 is available, the parameter of test is not meaningful.

The Hypothesis test 10 can output three kinds of figure on the basis of the result, including:

1) The LR profile for all chromosomes.\cr
2) The LR profile for QTL postion. \cr
3) The curve for QTL postion. \cr

The hypothesis test 15 can output a correlation figure for time*time matrix on the basis of the result.

The examples can be viewed in the following url.

\url{http://statgen.psu.edu/software/funmap/plot.t10-1.jpg}.\cr
\url{http://statgen.psu.edu/software/funmap/plot.t10-2.jpg}.\cr
\url{http://statgen.psu.edu/software/funmap/plot.t10-3.jpg}.\cr
\url{http://statgen.psu.edu/software/funmap/plot.t15.jpg}.\cr

}


\author{Zhong Wang, \email{zzw2@psu.edu} }

\examples{
par <- FM2.param( par_LC, CURVE_LC, CROSS_BC);
dat <- FM2.simulate(par);
ret <- FM2.qtlscan(dat);
plot(ret, dat);
}

\keyword{functions}
