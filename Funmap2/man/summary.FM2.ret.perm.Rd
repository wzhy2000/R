\name{summary.FM2.ret.perm}
\alias{summary.FM2.ret.perm}

\title{Summary of permutation result}

\description{ 
  Summarize the permutation result.
} 

\usage{ summary( perm_obj, file=NA, append=TRUE ) }

\arguments{
\item{perm_obj}{an object of permutation result returned by \code{\link{FM2.permutation}}.}
\item{file}{ an file name to save the summarized information.}
}

\details{

The summary command gives a table of cutoff values 
which starts at 90% and stops at 1% or even smaller.
If the p-value of x% want to be estimated, the 
permutation count should be greater than 100/x.

This summary can also export a figure based on the cutoff table.
The following gives an example of this summary command.

\preformatted{
Permutation result:   
------------------------------------
          Model: Logistic Curve
          Cross: BC
           Loop: 100
------------------------------------

p-value    Cutoff
0.90000    0.76583
0.80000    1.45845
0.70000    2.21657
0.60000    3.09488
0.50000    4.16930
0.40000    5.52038
0.30000    7.33172
0.20000    9.97031
0.10000   14.86003
0.09000   15.63818
0.08000   16.50743
0.07000   17.51136
0.06000   18.68964
0.05000   20.09660
0.04000   21.77278
0.03000   23.96467
0.02000   26.98845
0.01000   32.80592
 
*1:The figure 1 for the permutation result is saved to pheno.csv.1.pdf.
}

}

\author{Zhong Wang, \email{zzw2@psu.edu} }

\examples{
par <- FM2.param( par_LC, CURVE_LC, CROSS_BC);
dat <- FM2.simulate(par);
ret_perm <- FM2.permutation( dat, options=list(permu_loop=100) );
summary( ret_perm );
}

\keyword{summary}
