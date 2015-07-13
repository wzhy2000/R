summary.fGWAS.scan<-function(object, ...)
{
 	r.fgwas <- object;
 	#fgwas
 	#filter
 	#options
 	#params
 	#curve
 	#covariance
 	#est.values
 	
	r.sum.ret <- list();

	if(!is.null(r.gls$fgwas))
	{
		re7 <- r.gls$fgwas;
		fgwas.sig <- which( re7[,7] <= r.gls$options$fgwas.cutoff );
		if(length(fgwas.sig)>0)
		{
			fgwas_sigs <- re7[ fgwas.sig, , drop=F];
			fgwas.sig.inc <- order(fgwas_sigs[,7]);
			r.sum.ret$fgwas_sig <- fgwas_sigs[fgwas.sig.inc,];
		}
		
		if(!is.null(r.sum.ret$varsel))
			r.sum.ret$varsel <- cbind(r.sum.ret$varsel, fgwas.pvalue=find_fgwas_pvalue( r.gls$fgwas, rownames(r.sum.ret$varsel) ) ) ;

		if(!is.null(r.sum.ret$refit))
			r.sum.ret$refit <- cbind(r.sum.ret$refit, fgwas.pvalue=find_fgwas_pvalue( r.gls$fgwas, rownames(r.sum.ret$refit) ) ) ;
		
	}

	class(r.sum.ret) <- "sum.fGWAS.scan";
	
	r.sum.ret
}

print.sum.fGWAS.scan<-function(x, ...)
{
 	r.sum.ret <- x;

	if(!is.null(r.sum.ret$fgwas_sig))
	{
		cat("--- Significant SNPs Estimate by fGWAS method:", NROW(r.sum.ret$fgwas_sig), "SNPs\n");
		if( NROW(r.sum.ret$fgwas_sig)>25 )
		{
			cat("Top 25 SNPs:\n");
			show(r.sum.ret$fgwas_sig[1:25,,drop=F]);
		}
		else	
			show(r.sum.ret$fgwas_sig);
	}
}


plot.fGWAS.scan<-function( x, y=NULL, ... , fig.prefix=NULL )
{
	r.gls <- x;

	if( missing(fig.prefix)) fig.prefix <- "gls.plot";

	if(!is.null(r.gls$fgwas))
	{
		filter.man <- r.gls$fgwas[, c(1,2,7), drop=F]
		draw_man_fgwas( filter.man, fig.prefix, "fgwas" );
	}
	else
		cat("! No fGWAS filter results.\n");		
		
	if( !is.null(r.gls$varsel_add) || !is.null(r.gls$varsel_dom))
	{
		if ( !is.null(r.gls$varsel_add) )  varsel <- r.gls$varsel_add[, c(1,2), drop=F]
		if ( !is.null(r.gls$varsel_dom) )  varsel <- r.gls$varsel_dom[, c(1,2), drop=F]

		if ( !is.null(r.gls$varsel_add) ) varsel<- cbind( varsel, r.gls$varsel_add[,7] );
		if ( !is.null(r.gls$varsel_dom) ) varsel<- cbind( varsel, r.gls$varsel_dom[,7] );

		draw_man_adh2( varsel, fig.prefix, "varsel" );
	}
	else
		cat("! No varible selection results.\n");		

	if( !is.null(r.gls$refit_add) || !is.null(r.gls$refit_dom) )
	{
		refit<- merge_add_dom( r.gls$refit_add, r.gls$refit_dom);

		draw_refit_curve( refit, fig.prefix, "curve" );
	}
	else
		cat("! No refit results.\n");		
}


print.fGWAS.scan<-function(x, ...)
{
}

summary.fGWAS.dat<-function( x,..., fig.prefix=NULL )
{

}

summary.fGWAS.perm<-function( x,..., fig.prefix=NULL )
{

}

print.fGWAS.dat<-function( x,..., fig.prefix=NULL )
{

}

print.fGWAS.perm<-function( x,..., fig.prefix=NULL )
{

}

plot.fGWAS.dat<-function( x,..., fig.prefix=NULL )
{

}

plot.fGWAS.perm<-function( x,..., fig.prefix=NULL )
{

}