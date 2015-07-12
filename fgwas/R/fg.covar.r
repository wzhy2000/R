## Create a genral covariance model 
## 
## Two methods fortfbs.db
## tfbs.group: get the statistical summary( SQL/sum ) 
##             according to a field
## CisBP.create : get the subset by the specified parameter.
##

setClass("fgwas.covar", 
  representation(
    name        = "character",
	par_num    	= "integer", 
	par_name 	= "character",
	par_simu 	= "matrix",
	par_lower   = "numeric",
	par_upper   = "numeric"	)
  )

setGeneric("fgwas.covar.getMatrix", 
    def=function(covar.obj, pars, times) {
	  standardGeneric("fgwas.covar.getMatrix")
	})


## Create a CisBP database class (extending tfbs.db)

setClass("AR1.Covar", 
  representation(
    desc = "character"   
    ), contains = "fgwas.covar"
  )


setMethod("fgwas.covar.getMatrix", signature(covar.obj="AR1.Covar"),
    function( covar.obj, pars, times )
    {
		n <- length(times);
		rho<- pars[1];
		s2 <- pars[2];

		sigma <- array(0, dim=c(n,n));
		for(i in 1:n)
		{
			sigma[i,i:n] <- rho^abs( i - c(i:n) );
			sigma[i:n,i] <- sigma[i,i:n];				
		}
	
		sigma <- sigma*abs(s2);	

		return(sigma);
})



setClass("SAD.Covar", 
  representation(
    desc = "character"   
    ), contains = "fgwas.covar"
  )


setMethod("fgwas.covar.getMatrix", signature(covar.obj="SAD.Covar"),
    function( covar.obj, pars, times )
    {
		n <- length(times);
		phi<- pars[1];
		v2 <- pars[2];

		tmp <- (1-phi^2);
		sigma <- array(1, dim=c(n,n));
		for(i in 1:n)
		{
			#for(j in i:n)
			#{
			#	sigma[i,j] <- phi^(j-i) * (1-phi^(2*i))/(1-phi^2)				
			#	sigma[j,i] <- sigma[i,j];				
			#}
			sigma[i,i:n] <- phi^( c(i:n) - i ) * (1-phi^(2*i))/tmp;				
			sigma[i:n,i] <- sigma[i,i:n];				
		}

		sigma <- sigma*abs(v2);	
		return(sigma);
})


#------------------------------------------------------------------------------------------
# covariance
#------------------------------------------------------------------------------------------
fg.getCovar<- function( covar.name )
{
	for(i in 1:length(fgwas.sys$covars))
	{
		if ( fgwas.sys$covars[[i]]== covar.name ) 
			return( fgwas.sys$covars[[i]] );
	}

	return(NULL);	
}

ZZZ.regcovar<-function()
{
	fgwas.sys$covars <<- list();

	fgwas.sys$covars[[1]] <<- new("AR1.Covar", name="AR1", para_num=2, para_simu=NA, para_lower=c(0,0), para_upper=c(1,Inf)) ;
	cat(" ", fgwas.sys$covars[[1]]@name, " covariance is registered.\n");

	fgwas.sys$covariances[[2]] <<- new("SAD.Covar", name="SAD", para_num=2, para_simu=NA, para_lower=c(0,0), para_upper=c(2,Inf));
	cat(" ", fgwas.sys$covars[[2]]@name, " covariance is registered.\n");
}
