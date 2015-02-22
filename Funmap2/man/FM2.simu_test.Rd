\name{FM2.simu_test}
\alias{FM2.simu_test}

\title{Simulation test}

\description{ 
  Demostrate the simulation test for five curves as a main function.
} 

\usage{ FM2.simu_test( par_value, curve_type, cross_type, covar_type=COVAR_AR1 ) }

\arguments{
\item{par_value}{5 parameter value objects can be selected, see \code{\link{FM2.param}} for details. }
\item{curve_type}{5 curves can be specified, see \code{\link{FM2.param}} for details. }
\item{cross_type}{3 crosses can be specified, CROSS_BC, CROSS_F2 or CROSS_RIL. see \code{\link{FM2.param}} for details. }
\item{covar_type}{2 covariance type can be specified, COVAR_AR1 or COVAR_SAD2. see \code{\link{FM2.param}} for details. }
}

\value{
A list including raw data object(\strong{dat_st}), qtl scan result(\strong{ret_st}) and permutation result(\strong{ret_permu_st}).
}

\details{

\strong{FM2.simu_test} just demostrates how to use the Funmap to do a simulation test, 
it includes the following steps:

1) Create a parameter object according to the variables of model_type and cross_type.\cr
2) Simulate a raw data object on the basis of the parameter object.\cr
3) Perform the hypothesis tests.\cr
4) Execute permutation to get a cutoff for significant estimation.\cr
5) Store all objects(raw data, hypothesis result, permutation result).\cr
6) Summarize all objects and plot all figures.\cr


Three objects can be retrieved by the following variables.

1) Raw data object describled in \code{\link{FM2.simulate}}: \strong{dat_st};\cr
2) Hypothesis result object describled in \code{\link{FM2.qtlscan}}: \strong{ret_st};\cr
3) Permutation result object describled in \code{\link{FM2.permutation}}: \strong{ret_perm_st};\cr

It is very simple to name .rdata file where all objects are stored .
The .rdata filename starts with "simu_test_", 
and concacts the model type(LC, BC or PC) and cross type("BC" for backcross, "F2", "RIL") ,
and ends with ".rdata". For examples, a simulation test for Logistic curve and F2 intercross will 
have a rdata file named by "simu_test_LC_F2.rdara".


The source codes of \code{FM2.simu_test} are shown in the following.

\preformatted{

FM2.simu_test<-function( model_type, cross_type )
{
	if (is.null( FM_sys ) )
		FM2.start();

	par_st<- FM2.param( par_obj, curve_type, cross_type, covar_type );
	summary(par_st);
	
	dat_st<- FM2.simulate( par_st);
	summary(dat_st);

	ret_st<- FM2.qtlmodel( dat_st);
	summary(ret_st, dat_st);

	##Step 1
	rdata.name <- paste( "simu_test_",FM2.curve$name,"_",FM2.cross$name,".rdata", sep="");
	save(dat_st, ret_st, file = rdata.name );

	##Step 2(Report) 
	pdf.name <- paste( "simu_test_",FM2.curve$name,"_",FM2.cross$name,".pdf", sep="");
	FM2.report(pdf.name, dat_st, ret_st);

	##Step 3(Ppermutation) 
	FM2.set_value("debug", TRUE)
	ret_perm_st<- FM2.permutation( dat_st, options=list(permu_loop=100) );
	summary( ret_perm_st );
	try(plot(ret_perm_st));
	
	##Step 4(Save All Data)
	save(dat_st, ret_st, ret_perm_st, file = rdata.name );
	cat("\nNotice: The results(dat_st, ret_st and ret_perm_st) are stored into ", rdata.name, "\n");

	return( list(dat_st, ret_st, ret_perm_st) );
}
}


}

\author{Zhong Wang, \email{zzw2@psu.edu} }

\examples{

lst <- FM2.simu_test(par_LC, CURVE_LC, CROSS_BC);
#plot the data information.
plot(lst$dat_st);

}

\keyword{functions}
