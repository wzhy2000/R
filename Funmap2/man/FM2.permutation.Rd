\name{FM2.permutation}
\alias{FM2.permutation}

\title{Permutation}

\description{ 
  Execute permutation to get the cutoff value for significance levels p=0.05 and 0.01.
} 

\usage{ FM2.permutation(dat, options=list()) }

\arguments{
\item{dat}{A data object returned by \code{\link{FM2.simulate}} or \code{\link{FM2.load_data}} }
\item{options}{.Optional value for permutation control.\cr
cluster_count : default=1, the cluster count for parallel permutation.\cr
permu_loop: default=1000, the count of permutation loop.\cr
}
}

\value{
The p-value table will be returned. It has two columns, the first column is for significance level 
and the second column is for cutoff value.
}

\details{
If permutation count is greater than 100, more precise cutoff will be obtained.
For example, 10000 times permutation can give the significance table which
looks like the following table.

\tabular{rrrrrrrrr}{
0.9\tab0.8\tab0.7\tab0.6\tab0.5\tab0.4\tab0.3\tab0.2\tab0.1\cr
0.09\tab0.08\tab0.07\tab0.06\tab0.05\tab0.04\tab0.03\tab0.02\tab0.01\cr
0.009\tab0.008\tab0.007\tab0.006\tab0.005\tab0.004\tab0.003\tab0.002\tab0.001\cr
0.0009\tab0.0008\tab0.0007\tab0.0006\tab0.0005\tab0.0004\tab0.0003\tab0.0002\tab0.0001\cr
}

If clusters are available, the permutaion can use \strong{SNOW} package to do parallel computation.
In order to do that, the following is necessary.

1) \strong{SNOW} is installed.

2) The cluster count should be specified( \bold{FM2.set_value( "cluster_count", 2 )} ).

}

\author{Zhong Wang, \email{zzw2@psu.edu} }

\examples{
par <- FM2.param( par_LC, CURVE_LC, CROSS_BC);
dat <- FM2.simulate(par);
ret <- FM2.permutation(dat);
summary( ret );
plot(ret);
}

\keyword{functions}
