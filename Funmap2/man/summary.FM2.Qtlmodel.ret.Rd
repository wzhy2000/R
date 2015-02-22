\name{summary.FM2.Qtlmodel.ret}
\alias{ summary.FM2.Qtlmodel.ret }

\title{Summary of qtlscan tests.}

\description{ 
  Summarize information for the hypothesis test results.
} 

\usage{ summary(res_object, dat_obj ) }

\arguments{
\item{res_obj}{a result object returned by \code{\link{FM2.qtlscan}} which stores the results of hypothesis tests.}
\item{dat_obj}{a data object returned by \code{\link{FM2.simulate}} or \code{\link{FM2.load_data}}.}
}

\details{
For logistic curve, there are at most 6 hypothesis tests and for other two curves 
only 1 hypothesis test. When summarization is performed,
the report is exported to the R console or log file respectively.

The following sections demostrates the context of summary report.

\bold{Hypothesis Test 10}

Three parts are included. 

1) Estimated parameters. \cr
2) The significant QTL postions. \cr
3) Three kinds of figures are drawn on the basis of the result, including:\cr

a) The LR profile for all chromosomes.\cr
b) The LR profile for QTL postion. \cr
c) The curve for QTL postion. \cr

\preformatted{
Hypothesis test 10: 
    a2=a1 and b2=b1 and r2=r1
------------------------------------
          Model: Logistic Curve
          Cross: Backcross
       QTL pos.: 50.1  (Group:8)
         QTL LR: 66.516   
    QTL p-value: 0.000    
Grwoth para(Qq): a2=  30.615, b2=  10.776, r2=   0.538
Grwoth para(qq): a1=  23.707, b1=   9.449, r1=   0.615
            rho: 0.953   
         sigma2: 8.637   
------------------------------------
 No.  Grp      Pos.       LR        a1        b1        r1        a0        b0        r0
   1    8   50.100    66.516    30.615    10.776     0.538    23.707     9.449     0.615
   2   12  113.100    55.190    29.865     9.736     0.528    25.207     8.725     0.586
   3   13   12.000    50.963    29.518     9.723     0.526    24.926     8.906     0.602
   4   18   10.000    25.684    29.236     9.550     0.536    25.523     8.963     0.584
   5    1  151.300    24.162    25.998     8.520     0.575    28.801     9.672     0.536

*1:The figure 1 for the hypothesis test(10) is saved to pheno.csv.1.pdf.
*2:The figure 2 for the hypothesis test(10) is saved to pheno.csv.2.pdf.
*2:The figure 2 for the hypothesis test(10) is saved to pheno.csv.3.pdf.
*2:The figure 2 for the hypothesis test(10) is saved to pheno.csv.4.pdf.
*2:The figure 2 for the hypothesis test(10) is saved to pheno.csv.5.pdf.
*2:The figure 2 for the hypothesis test(10) is saved to pheno.csv.6.pdf.
*3:The figure 3 for the hypothesis test(10) is saved to pheno.csv.7.pdf.
}


}


\author{Zhong Wang, \email{zzw2@psu.edu} }

\examples{
par <- FM2.param(par_LC, CURVE_LC, CROSS_BC);
dat <- FM2.simulate(par);
ret <- FM2.qtlscan(dat);
summary( ret, dat );
}

\keyword{summary}
