\name{summary.FM2.dat}
\alias{summary.FM2.dat}

\title{Summary of raw data object}

\description{ 
  Summarize information for the raw data object.
} 

\usage{ summary( dat_obj, file=NA , append = TRUE ) }

\arguments{
\item{dat_obj}{a data object return by \code{\link{FM2.simulate}} or \code{\link{FM2.load_data}} }
\item{file}{a file name to save summarized information. }
}

\details{

The data format is described in \code{\link{FM2.simulate}}.

The following example demonstrates summary command for a data object.

\preformatted{
Data set for FunMap model:
------------------------------------
          Date: 2010-03-19 03:49:36
          Model: Logistic Curve
          Cross: F2
    Pheno. file: simu.pheno.LC.F2
     Geno. file: simu.geno.LC.F2
     Maker file: simu.marker.LC.F2
    Sample size: 100       
   Sample times: 7         
   Marker count: 6         
          LC  a: 19.83678  
              b: 8.96370   
              r: 0.47202   
            rho: 0.75430   
         sigma2: 0.58849   
------------------------------------
*1:The figure 1 for all individuals is saved to simu.pheno.LC.F2.3.pdf.
*2:The figure 2 for all individuals is saved to simu.pheno.LC.F2.4.pdf.
}

While the summary is outputed on the R console, two figures are drawn and saved to two files.

}


\author{Zhong Wang, \email{zzw2@psu.edu} }

\examples{
par <- FM2.param(par_LC, CURVE_LC, CROSS_BC);
dat <- FM2.simulate(par);
summary( dat );
}

\keyword{functions}
