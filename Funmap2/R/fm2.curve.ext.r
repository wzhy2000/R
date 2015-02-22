#################################################################
#
# New Systems Mapping Application(Funmap2)
#
# Implemented Curve:
#
# CURVE_LC:  Logistic Curve
#
# History:
# 02/15/2012 Version 1.0
#
##################################################################

#-----------------------------------------------------------------
# CURVE_LC:  Logistic Curve
#-----------------------------------------------------------------

LC.get_mu <- function(par, times, options=list())
{
	return ( par[[1]]/(1+par[[2]]*exp(-par[[3]]*times) ) );
}

class_LC<-list(
	type 	 	= 1,
	name 	 	= "LC",
	desc 	 	= "Logistic Curve",
	par_num    	= 3, 
	par_name 	= c( "a", "b", "r"),
	trait_num  	= 1,
	get_mu      	= LC.get_mu );


par_LC<-list(
	#sample size
	simu_N      = 100,
	#trait measured in time 1 to 7
	simu_times  = c(1:7),
	#marker distance
	simu_mrkdist= c(0, 20, 20, 20, 20, 20),
	#qtl position
	simu_qtlpos = 50,			

	simu_cross = list(),
	
	simu_covar = list(),
	
	#AR1 model
	covar_ar1 = list(
		#rho
		rho   = 0.7543,
		s2    = 1 ),
	
	# 	QQ2         
	QQ2 = list(
		simu_a   = 21.9824,
		simu_b   = 9.7768,
		simu_r   = 0.4699),

	# 	Qq1         
	Qq1 = list(
		simu_a   = 19.9810,
		simu_b   = 8.776,
		simu_r   = 0.47),
	
	#	qq0
	qq0 = list(
		simu_a   = 15.95,
		simu_b   = 7.58,
		simu_r   = 0.48)
);

#-----------------------------------------------------------------
# CURVE_BI:  Bi-exponential Curve
#-----------------------------------------------------------------
BI.get_mu <- function(par, times, options=list())
{
	return(  par[[1]]*exp(-par[[2]]*times) + par[[3]]*exp(-par[[4]]*times) );
}

class_BI<-list(
	type 	 	= 1,
	name 	 	= "BI",
	desc 	 	= "Bi-exponential Curve",
	par_num    	= 4, 
	par_name 	= c( "p1", "d1", "p2", "d2"),
	trait_num  	= 1,
	get_mu      	= BI.get_mu );


par_BI<-list(
	#sample size
	simu_N      = 100,
	#trait measured in time 1 to 7
	simu_times  = c(1:7),
	#marker distance
	simu_mrkdist= c(0, 20, 20, 20, 20, 20),
	#qtl position
	simu_qtlpos = 50,			

	simu_cross = list(),
	
	simu_covar = list(),
	
	#AR1 model
	covar_ar1 = list(
		#rho
		rho   = 0.7543,
		s2    = 1 ),
	
	# 	QQ2         
	QQ2 = list(
		simu_p1  = 19.9824,
		simu_d1  = 0.4699,
		simu_p2  = 8.7768,
		simu_d2  = 1.4699),

	# 	Qq1         
	Qq1 = list(
		simu_p1  = 17.9824,
		simu_d1  = 0.0699,
		simu_p2  = 9.7768,
		simu_d2  = 1.0699),
	
	#	qq0
	qq0 = list(
		simu_p1  = 15.9507,
		simu_d1  = 0.1836,
		simu_p2  = 10.5737,
		simu_d2  = 1.8836)
);

#-----------------------------------------------------------------
# CURVE_PC:  Logistic Curve
#-----------------------------------------------------------------
PC.get_mu <- function(par, times, options=list())
{
	return(  par[[1]] + par[[3]]*times/(par[[2]] + times) );
}

class_PC<-list(
	type 	 	= 1,
	name 	 	= "PC",
	desc 	 	= "Pharmacology Curve",
	par_num    	= 3, 
	par_name 	= c( "E0", "E50", "Emax" ),
	trait_num  	= 1,
	get_mu      	= PC.get_mu );


par_PC<-list(
	#sample size
	simu_N      = 100,
	#trait measured in time 1 to 7
	simu_times  = c(1:7),
	#marker distance
	simu_mrkdist= c(0, 20, 20, 20, 20, 20),
	#qtl position
	simu_qtlpos = 50,			

	simu_cross = list(),
	
	simu_covar = list(),
	
	#AR1 model
	covar_ar1 = list(
		#rho
		rho   = 0.7543,
		s2    = 1 ),
	
	# 	QQ2         
	QQ2 = list(
		simu_E0  = 10.9824,
		simu_E50 = 15.909,
		simu_Emax= 20.7768),

	# 	Qq1         
	Qq1 = list(
		simu_E0  = 8.9824,
		simu_E50 = 16.098,
		simu_Emax= 20.7768),
	
	#	qq0
	qq0 = list(
		simu_E0  = 6.9507,
		simu_E50 = 12.090,
		simu_Emax= 18.5737)
);


#-----------------------------------------------------------------
# CURVE_EXP:  Exponentiation Curve
#-----------------------------------------------------------------
EXP.get_mu <- function(par, times, options=list())
{
	return(  par[[1]]*( times^par[[2]]) );
}

class_EXP <- list(
	type 	 	= 1,
	name 	 	= "EXP",
	desc 	 	= "Exponentiation Curve",
	par_num    	= 2, 
	par_name 	= c( "a", "b" ),
	trait_num  	= 1,
	get_mu      	= EXP.get_mu );


par_EXP<-list(
	#sample size
	simu_N      = 100,
	#trait measured in time 1 to 7
	simu_times  = c(1:7),
	#marker distance
	simu_mrkdist= c(0, 20, 20, 20, 20, 20),
	#qtl position
	simu_qtlpos = 50,			

	simu_cross = list(),
	
	simu_covar = list(),
	
	#AR1 model
	covar_ar1 = list(
		#rho
		rho   = 0.7543,
		s2    = 1 ),
	
	# 	QQ2         
	QQ2 = list( simu_a = 11.049,
		    simu_b = 1.151),

	# 	Qq1         
	Qq1 = list( simu_a = 9.049,
		    simu_b = 1.251),
	
	#	qq0
	qq0 = list( simu_a = 7.148,
		    simu_b = 1.359)
);


#-----------------------------------------------------------------
# CURVE_NP:  Nonparametric Curve
#-----------------------------------------------------------------


# private: Legendre.model
Legendre.model<-function( t, mu, tmin=NULL, tmax=NULL )
{
	u <- -1; 
	v <- 1;
	if (is.null(tmin)) tmin<-min(t);
	if (is.null(tmax)) tmax<-max(t);
	
	ti    <- u + ((v-u)*(t-tmin))/(tmax - tmin);
	np.order <- length(mu)-1;

	L <- mu[1] + ti*mu[2];
	if (np.order>=2)
		L <- L + 0.5*(3*ti*ti-1)* mu[3] ;
	if (np.order>=3)
		L <- L + 0.5*(5*ti^3-3*ti)*mu[4] ;
	if (np.order>=4)
		L <- L + 0.125*(35*ti^4-30*ti^2+3)* mu[5];
	if (np.order>=5)
		L <- L + 0.125*(63*ti^5-70*ti^3+15*ti)*mu[6];
	if (np.order>=6)
		L <- L + (1/16)*(231*ti^6-315*ti^4+105*ti^2-5)* mu[7];
	if (np.order>=7)
		L <- L + (1/16)*(429*ti^7-693*ti^5+315*ti^3-35*ti)* mu[8];
	if (np.order>=8)
		L <- L + (1/128)*(6435*ti^8-12012*ti^6+6930*ti^4-1260*ti^2+35)* mu[9];
	if (np.order>=9)
		L <- L + (1/128)*(12155*ti^9-25740*ti^7+18018*ti^5-4620*ti^3+315*ti)* mu[10];
	if (np.order>=10)
		L <- L + (1/256)*(46189*ti^10-109395*ti^8+90090*ti^6-30030*ti^4+3465*ti^2-63)* mu[11];
	if (np.order>=11)
	{
		for(r in 11:(np.order))
		{
			kk <- ifelse(r%%2==0, r/2, (r-1)/2);
			for (k in c(0:kk) )
			{
				L <- L + (-1)^k*factorial(2*r-2*k)/factorial(k)/factorial(r-k)/factorial(r-2*k)/(2^r)*ti^(r-2*k)*mu[r+1];
			}
		}
	}
	
	
	return(L);
}


NP.get_mu <- function(par, times, options=list())
{
	return( Legendre.model(times, unlist(par))  );
}

class_NP <- list(
	type 	 	= 1,
	name 	 	= "NP",
	desc 	 	= "Nonparametric Curve",
	par_num    	= 5, 
	par_name 	= c( "u0", "u1", "u2", "u3", "u4"),
	trait_num  	= 1,
	get_mu      	= NP.get_mu );

par_NP<-list(
	#sample size
	simu_N      = 100,
	#trait measured in time 1 to 7
	simu_times  = c(1:7),
	#marker distance
	simu_mrkdist= c(0, 20, 20, 20, 20, 20),
	#qtl position
	simu_qtlpos = 50,			

	simu_cross = list(),
	
	simu_covar = list(),
	
	#AR1 model
	covar_ar1 = list(
		#rho
		rho   = 0.7543,
		s2    = 1 ),
	
	# 	QQ2         
	QQ2 = list(simu_u0 = 11.049,
		   simu_u1 = 1.551,
		   simu_u2 = -8.019,
		   simu_u3 = 3.151,
		   simu_u4 = 0.652,
		   simu_u5 = -0.597,
		   simu_u6 = 0.8211),

	# 	Qq1         
	Qq1 = list( simu_u0 = 9.049,
		    simu_u1 = 1.151,
		    simu_u2 = -6.019,
		    simu_u3 = 2.651,
		    simu_u4 = 0.652,
		    simu_u5 = -0.797,
		    simu_u6 = 0.621),
	
	#	qq0
	qq0 = list( simu_u0 = 7.148,
		    simu_u1 = 1.379,
		    simu_u2 = -4.489,
		    simu_u3 = 2.004,
		    simu_u4 = 0.662,
		    simu_u5 = -0.836,
		    simu_u6 = 0.432)
);


ZZZ.regcurve<-function()
{
	CURVE_LC <<- FM2.reg_curve( class_LC );
	CURVE_BI <<- FM2.reg_curve( class_BI );
	CURVE_PC <<- FM2.reg_curve( class_PC );
	CURVE_EXP <<- FM2.reg_curve( class_EXP );
	CURVE_NP <<- FM2.reg_curve( class_NP );
}