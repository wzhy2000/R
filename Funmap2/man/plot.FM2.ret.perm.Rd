\name{plot.FM2.ret.perm}
\alias{plot.FM2.ret.perm(perm_obj) }

\title{Figure of permutation result}

\description{ 
  Draw cutoff profile on the basis of the permutation result.
} 

\usage{ plot(perm_obj) }

\arguments{                                                  
\item{perm_obj}{an object of permutation result returned by \code{\link{FM2.permutation}}.}
}

\details{

This figure can also be exported by the summary command, see (\code{\link{summary.FM2.ret.perm}}) for details. 
There is a little difference between the summary command and the plot command. 
The summary command can generate this figure to a PDF file, 
but the plot command can output this figure in the R console.

An example of this command is available in the following URL.

\url{http://statgen.psu.edu/software/funmap/plot.perm.jpg}.

}


\author{Zhong Wang, \email{zzw2@psu.edu} }

\examples{
par <- FM2.param(par_LC, CURVE_LC, CROSS_BC);
dat <- FM2.simulate(par);
ret_perm <- FM2.permutation(dat, options=list(permu_loop=100));
plot(ret_perm);
}

\keyword{functions}
