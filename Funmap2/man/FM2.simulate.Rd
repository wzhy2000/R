\name{FM2.simulate}
\alias{FM2.simulate}

\title{Data simulation}

\description{ 
  Create a simulation data object for model demonstration.
} 

\usage{ FM2.simulate(par_obj) }

\arguments{
\item{par_obj}{a parameter object. see \code{\link{FM2.param}} for details. }
}

\value{
A data object can be returned, its structure is same as the real data.
}

\details{

The structure of simulation data is identical to real data object. It has the following structure.

\preformatted{
$ name        : String "LC.dat" (it is not important for the user)
$ sample_N    : sample size
$ sample_times: sample times
$ pheno_file  : string "simu.pheno.LC.F2"
$ geno_file   : string "simu.geno.LC.F2"
$ marker_file : string "simu.marker.LC.F2"
$ phenos_table: a matrix of dim=c(sample_N, sample_times) 
$ genos_table : a matirx of dim=c(sample_N, makers count) 
$ marker_table: data.frame of dim=c(makers count, 4)
  ..$ Marker : marker name 
  ..$ Dist   : marker distance
  ..$ grp_idx: group index
  ..$ Group  : group name
}

1) The items of pheno_file, geno_file and marker_file are 
made up by the Funmap package and will be used to assign 
the output filename as a filename prefix.

2) In the phenos_table, the genotype is coded by 1=Qq 2=QQ for backcross 
, 0=qq, 1=Qq, 2=QQ for F2 intercross and 0=qq, 2=QQ for RILs intercross..

}

\author{Zhong Wang, \email{zzw2@psu.edu} }

\examples{
par <- FM2.param( par_LC, CURVE_LC, CROSS_BC);
dat <- FM2.simulate(par);
#summarize the data information.
summary( dat );
plot(dat);
str(dat);
}

\keyword{functions}
