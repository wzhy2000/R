\name{FM2.qtlmodel}
\alias{FM2.qtlmodel}

\title{ Hypothesis tests}

\description{ 
  Perform the hypothesis tests for a data object.
} 

\usage{ FM2.qtlmodel( dat_obj) }

\arguments{
\item{dat_obj}{a data object returned by \code{\link{FM2.simulate}} or \code{\link{FM2.load_data}}}
\item{test}{indicates which null hypothesis will be performed. \cr 
For logistic curve(CURVE_LC), valid values are 10 to 15. \cr
Except logistic curve(CURVE_LC), valid value only is 10.}
}

\value{
A result object which can be summarized by the summary(\code{\link{summary.FM2.ret.hp}}) command.
}

\details{

\bold{Hypothesis Test 10:}\cr
\strong{For different genotypes, all parameters are identical.}

The hypothesis test 10 scans every marker by the specified step (2cm). It maybe take a long time, so the Funmap package displays its progress after each chromosome (linkage group) has been calculated.

After the marker scanning, the package identifies the five most significant QTLs. In one choromosome (group) at most one significant QTL is selected. The most significant QTL is strongly displayed at the head of report. 

}


\author{Zhong Wang, \email{zzw2@psu.edu} }

\examples{
par <- FM2.param(par_LC, CURVE_LC, CROSS_BC);
dat <- FM2.simulate(par);
ret <- FM2.qtlmodel(dat);
summary(ret);
}

\keyword{functions}
