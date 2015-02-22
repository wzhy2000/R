\name{summary.FM2.par}
\alias{ summary.FM2.par }

\title{Summary of parameter object}

\description{ 
  Summarized a parameter object by summary command.
} 

\usage{ summary( par_obj ) }

\arguments{
\item{par_obj}{a parameter object returned by \code{\link{FM2.param}} }
}

\details{

Below is an output example for a parameter object returned by \code{\link{FM2.param}}(par_LC, CURVE_LC, CROSS_F2).

\preformatted{
------------------------------------
           Date: 2010-03-19 03:45:13
          Model: Logistic Curve
          Cross: F2
    Sample size: 100       
   Sample times: 7         
    Marker pos.: 0,20,40,60,80,100
       QTL pos.: 50        
       LC(QQ) a: 21.98240  
              b: 9.77680   
              r: 0.46990   
       LC(Qq) a: 19.98240  
              b: 8.77680   
              r: 0.46990   
       LC(qq) a: 15.95070  
              b: 7.57370   
              r: 0.48360   
            rho: 0.75430   
         sigma2: 0.22600   
------------------------------------
}

}

\author{Zhong Wang, \email{zzw2@psu.edu} }

\examples{
par <- FM2.param(par_LC, CURVE_LC, CROSS_F2);
summary( par )
}

\keyword{summary}
