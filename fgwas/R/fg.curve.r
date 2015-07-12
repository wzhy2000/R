## Create a genral curve model 
## 
## Two methods fortfbs.db
## tfbs.group: get the statistical summary( SQL/sum ) 
##             according to a field
## CisBP.create : get the subset by the specified parameter.
##

setClass("fgwas.curve", 
  representation(
    name        = "character",
	par_num    	= "integer", 
	par_name 	= "character",
	par_simu 	= "matrix" )
  )

setGeneric("fgwas.curve.getMu", 
    def=function(curve.obj, param, times) {
	  standardGeneric("fgwas.curve.getMu")
	})


## Create a CisBP database class (extending tfbs.db)
##  
## Download Link: http://cisbp.ccbr.utoronto.ca/bulk_archive.php
## paper: Determination and inference of eukaryotic transcription 
##        factor sequence specificity, Cell, 2014

setClass("Logistic.Curve", 
  representation(
    desc = "character"   
    ), contains = "fgwas.curve"
  )


setMethod("fgwas.curve.getMu", signature(curve.obj="Logistic.Curve"),
    function( curve.obj, param, times )
    {
	  y <- param[1]/(1+param[2]*exp(-param[3]*times) ) 
      return(y);	    
})


setClass("Logistic2.Curve", 
  representation(
    desc = "character"   
    ), contains = "fgwas.curve"
  )

setMethod("fgwas.curve.getMu", signature(curve.obj="Logistic2.Curve"),
    function( curve.obj, param, times )
    {
	  y <- par[1]/(1+par[2]*exp(-par[3]*times) ) + par[4]/(1+par[5]*exp(-par[6]*times) ) 
      return(y);	    
})


setClass("ABRK.Curve", 
  representation(
    desc = "character"   
    ), contains = "fgwas.curve"
  )


setMethod("fgwas.curve.getMu", signature(curve.obj="ABRK.Curve"),
    function( curve.obj, param, times )
    {
  	  y<- par[1]*( 1 + par[2]*exp(-par[3]*times) )^(1/(1-par[4]) );
      return(y);	    
})


fg.getCurve<- function( curve.name )
{
	for(i in 1:length(fgwas.sys$curves))
	{
		if ( fgwas.sys$curves[[i]]== curve.name ) 
			return( fgwas.sys$curves[[i]] );
	}

	return(NULL);	
}

ZZZ.regcurve<-function()
{
	fgwas.sys$curves <<- list();

	fgwas.sys$curves[[1]] <<- new("Logistic.Curve", name="Logistic");
	cat(" ", fgwas.sys$curves[[1]]@name, " curve is registered.\n");

	fgwas.sys$curves[[2]] <<- new("Logistic2.Curve", name="Logistic2");
	cat(" ", fgwas.sys$curves[[2]]@name, " curve is registered.\n");

	fgwas.sys$curves[[3]] <<- new("ABRK.Curve", name="ABRK");
	cat(" ", fgwas.sys$curves[[3]]@name, " curve is registered.\n");
}

if(0)
{
    Mu_logistic2 <<- list( type="CURVE",
    			 name="Double Logistic", 
    			 func=get_log2_mu, 
    			 pnum=6, 
    			 lower=c(-Inf,-Inf,-Inf,-Inf,-Inf,-Inf), 
    			 upper=c(Inf,Inf,Inf,Inf,Inf,Inf) );

    Mu_logistic  <<- list( 
    		        type="CURVE",
    			name="Logistic Curve", 
    			func=get_log_mu, 
    			pnum=3, 
    			par=c(), 
    			simu_par=list(QQ=c(14.5, 0.8, 2.7), Qq=c(15.5, 0.75, 2.6), qq=c( 16.5, 0.77, 2.5)), 
    			lower=c(-Inf,-Inf,-Inf), 
    			upper=c(Inf,Inf,Inf) );
    			
    Mu_abrk      <<- list( 
    			type="CURVE",
    			name="Growth Curve(ABRK)", 
    			func=get_abrk_mu, 
    			pnum=4, 
    			par=c(), 
    			lower=c(-Inf,-Inf,-Inf), 
    			upper=c(Inf,Inf,Inf) );
}