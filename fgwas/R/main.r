fgwas.simple<-function( file.phe, file.snp, Y.prefix, Z.prefix, covariate.names, curve=NA, covariance=NA, fgwas.filter=FALSE, options=NULL )  
{
	cat( "[ fGWAS SIMPLE ] Procedure.\n");
	cat( "Checking the parameters ......\n");

	if ( missing(file.phe) || missing(file.snp) || missing(Y.prefix) || missing(Z.prefix) || missing(covariate.names) )
         stop("! file.phe, file.snp, Y.prefix, Z.prefix and covariate.names must be assigned with the valid values.");

	if ( !(is.character(Y.prefix) && length(Y.prefix)==1 ) )
		stop("! The parameter of Y.prefix should be assigned with a prefix of outcome column in the phenotypic data.");
	if ( !(is.character(Z.prefix) && length(Z.prefix)==1 ) )
		stop("! The parameter of Z.prefix should be assigned with a prefix of time column in the phenotypic data.");
	if ( !missing("covariate.names") && length(covariate.names)>0 && !is.character(covariate.names) )
		stop("! The parameter of covariate.names should be assigned with covariate names in the phenotypic data.");
	if ( !(is.logical(fgwas.filter) && length(fgwas.filter)==1 ) )
		stop("! The parameter of fgwas.filter should be a logical value(TRUE or FALSE).");
	
	cat("* Phenotypic Data File = ",  file.phe, "\n");
	cat("* Simpe SNP File = ",  file.snp, "\n");

	show_fgwas_parameters( curve, covariance, Y.prefix, Z.prefix, covariate.names, fgwas.filter ) ;

	if (missing(options)) 
		options <- get_default_options()
	else	
	{
		options0 <- get_default_options();
        options0[names(options)] <- options;
        options <- options0;
    }
	
	cat( "Checking the optional items......\n");
	show_options( options );
	
	r.simple <- list();
	r.simple$options <- options;
	r.simple$params <- list( 
				file.phe       = file.phe, 
				file.snp       = file.snp, 
				Y.prefix       = Y.prefix, 
				Z.prefix       = Z.prefix, 
				covariate.names= covariate.names, 
				fgwas.filter   = fgwas.filter);

	r.phe <- read_simple_phenotype( file.phe, Y.prefix, Z.prefix, covariate.names );	
	if(r.phe$error)
		stop(r.phe$err.info);
	
	r.simpe$phe.mat <- r.phe$phe.mat;


	r.est  <- fg.estimate( r.simple$phe.mat, curve, covariance )
	if( r.est$error )
		stop(r.est$err.info)
	else
	{
		if(is.na(curve)) curve <- r.est$curve;	
		if(is.na(covariance)) covariance <- r.est$covariance;

		r.simple$curve <- curve;
		r.simple$covariance <- covariance;
		r.simple$est.values <- r.est$est.values;
	}

	r.snp <- read_simple_genotype( file.snp, r.simpe$phe.mat );	
	if(r.snp$error)
		stop(r.snp$err.info);

	r.simpe$snp.mat <- r.snp$snp.mat;

	if(fgwas.filter)
	{
		r.filter <- snpmat_lm_filter( r.simple$phe.mat, r.simple$snp.mat, Y.prefix, Z.prefix, covariate.names, options$nParallel.cpu, options$fgwas.cutoff, "LONG");
		if( r.filter$error )
			stop(r.filter$err.info);
			
		
		r.simple$filter <- r.filter$r.fgwas;
		r.simple$org.snp.mat <- r.simple$snp.mat;
		r.simple$snp.mat <- r.filter$snp.mat,
	}

	subset_op <- function(snpmat, sub.idx)
	{
		return( snpmat[sub.idx,,drop=F] );
	}
	
	r.fgwas <- fgscan_parallel( 
				NROW( r.simple$snp.mat ),
				subset_op,
				r.simple$phe.mat,
				Y.prefix, 
				Z.prefix,
				covariate.names,
				options$debug,
				options$nParallel.cpu );
	
	#reduce size
	r.simple$org.snp.mat <- NULL;
	r.simple$snp.mat <- NULL;
	
	if(!is.null(r.fgwas) && !is.na(r.fgwas) )
	{
		r.simple$fgwas <- r.fgwas;
		return(r.simple);		   
	}
	else
	{
		cat("! No results\n");

		r.simple$fgwas <- "try-error";
		return(r.simple);		   
	}
}

fgwas.plink<-function( file.phe, file.plink.bed, file.plink.bim, file.plink.fam, Y.prefix, Z.prefix, covariate.names, curve=NA, covariance=NA, fgwas.filter=FALSE, options=NULL, force.split=FALSE, plink.command=NULL )        
{
	cat( "[ fGWAS PLINK ] Procedure.\n");
	cat( "Checking the parameters ......\n");

	if ( missing(file.phe) || missing(file.plink.bed) || missing(file.plink.bim) || missing(file.plink.fam) || 
         missing(Y.prefix) || missing(Z.prefix) || missing(covariate.names) )
        stop("! file.phe, file.plink.bed, file.plink.bim, file.plink.fam, Y.prefix, Z.prefix and covariate.names must be assigned with the valid values.");

	if ( !(is.character(Y.prefix) && length(Y.prefix)==1 ) )
		stop("! The parameter of Y.prefix should be assigned with a prefix of outcome column in the phenotypic data.");
	if ( !(is.character(Z.prefix) && length(Z.prefix)==1 ) )
		stop("! The parameter of Z.prefix should be assigned with a prefix of time column in the phenotypic data.");
	if ( !missing("covariate.names") && length(covariate.names)>0 && !is.character(covariate.names) )
		stop("! The parameter of covariate.names should be assigned with covariate names in the phenotypic data.");
	if ( !(is.logical(fgwas.filter) && length(fgwas.filter)==1 ) )
		stop("! The parameter of fgwas.filter should be a logical value(TRUE or FALSE).");
	if ( !(is.logical(force.split) && length(force.split)==1 ) )
		stop("! The parameter of force.split should be a logical value(TRUE or FALSE).");

	cat("* Phenotypic Data File = ",  file.phe, "\n");
	cat("* PLINK BED File = ",  file.plink.bed, "\n");
	cat("* PLINK BIM File = ",  file.plink.bim, "\n");
	cat("* PLINK FAM File = ",  file.plink.fam, "\n")
	cat("* PLINK Command = ",   plink.command, "\n")
	cat("* Force Split by PLINK Command = ", force.split, "\n")
	
	show_fgwas_parameters( curve, covariance, Y.prefix, Z.prefix, covariate.names, fgwas.filter ) ;

	if (missing(options)) 
		options <- get_default_options()
	else	
	{
		options0 <- get_default_options();
        options0[names(options)] <- options;
        options <- options0;
    }
	
	cat( "Checking the optional items......\n");
	show_options( options);
	
	r.plink <- list();
	r.plink$options <- options;
	r.plink$params <- list( 
				file.phe       = file.phe, 
				file.plink.bed = file.plink.bed, 
				file.plink.bim = file.plink.bim, 
				file.plink.fam = file.plink.fam,				
				Y.prefix       = Y.prefix, 
				Z.prefix       = Z.prefix, 
				covariate.names= covariate.names, 
				fgwas.filter   = fgwas.filter);


	r.phe <- read_simple_phenotype( file.phe, Y.prefix, Z.prefix, covariate.names );	
	if(r.phe$error)
		stop(r.phe$err.info);
	
	r.plink$phe.mat <- r.phe$phe.mat;

	r.est  <- fg.estimate( r.plink$phe.mat, Y.prefix, Z.prefix, covariate.names, curve, covariance )
	if( r.est$error )
		stop(r.est$err.info);
	else
	{
		if(is.na(curve)) curve <- r.est$curve;	
		if(is.na(covariance)) covariance <- r.est$covariance;

		r.plink$curve <- curve;
		r.plink$covariance <- covariance;
		r.plink$est.values <- r.est$est.values;
	}
	
	pd <- list();
	if( force.split || !try_load_plink( file.plink.bed,  file.plink.bim, file.plink.fam ) )
	{
		# It is bigdata which need to split it into chromosome unit
	    # The following will split the data and force to do fGWAS filter.

	    r.filter <- plink_fgwas_bigdata ( file.plink.bed,  file.plink.bim, file.plink.fam, file.phe, plink.command, 
	    						          Y.prefix, Z.prefix, covariate.names, options$nParallel.cpu, options$fgwas.cutoff, "GLS");
		if( r.filter$error )
			stop(r.filter$err.info);

	    fgwas.filter <- TRUE;
		r.plink$filter <- r.filter$r.fgwas;
	    pd <- r.filter$snp.mat;
	}	
	else
	{
		pd <- load_plink_binary( file.plink.bed,  file.plink.bim, file.plink.fam, file.phe );
		if( is.null(pd) )
			stop("Failed to load PLINK dataset!");

		if(fgwas.filter)
		{
			# call FGWAS.R to do FILTER and the gls__snpmat
			r.filter <- plink_fgwas_filter( pd, Y.prefix, Z.prefix, covariate.names, options$nParallel.cpu, options$fgwas.cutoff, "GLS");
			if( r.filter$error )
				stop(r.filter$err.info);
			
			pd <- r.filter$snp.mat;

			r.plink$filter <- r.filter$r.fgwas;
		}				
	}
	
	if( fgwas.filter)
	{
		subset_op <- function(snpmat, sub.idx)
		{
			return( snpmat[sub.idx,,drop=F] );
		}

		r.fgwas <- fgscan_parallel( 
				NROW( pd ),
				subset_op,
				pd,
				r.plink$phe.mat,
				Y.prefix, 
				Z.prefix,
				covariate.names,
				options$debug,
				options$nParallel.cpu );

	}
	else
	{
		subset_op <- function(snpmat, sub.idx)
		{
			snp.sub <- get_plink_subsnp(snpmat, sub.idx );
			snp.mat <- cbind( snp.sub$info[,c(2,3)], snp.sub$snp )
			return( snp.mat );
		}
		
		r.fgwas <- fgscan_parallel( 
			NCOL( pd$snp.mat$genotypes ),
			subset_op,
			pd$snp.mat,
			pd$phe.mat,
  		   	Y.prefix, 
  		   	Z.prefix,
  		  	covariate.names,
			options$debug,
			options$nParallel.cpu );
	}
	
	if(!is.null(r.fgwas) && !is.na(r.fgwas) )
	{
		r.plink$fgwas <- r.fgwas;
		return(r.plink);		   
	}
	else
	{
		cat("! No results\n");

		r.plink$fgwas <- "try-error";
		return(r.plink);		   
	}

}

fgwas.snpmat<-function( phe.mat, snp.mat, Y.prefix, Z.prefix, covariate.names, curve=NA, covariance=NA, fgwas.filter=FALSE, options=NULL)     
{
	cat( "[ fGWAS SNPMAT ] Procedure.\n");
	cat( "Checking the parameters ......\n");

	if ( missing(phe.mat) || missing(snp.mat) || missing(Y.prefix) || missing(Z.prefix) || missing(covariate.names) )
             stop("! phe.mat, snp.mat, file.plink.tfam, Y.prefix, Z.prefix and covariate.names must be assigned with the valid values.");

	if ( !(is.character(Y.prefix) && length(Y.prefix)==1 ) )
		stop("! The parameter of Y.prefix should be assigned with a prefix of outcome column in the phenotypic data.");
	if ( !(is.character(Z.prefix) && length(Z.prefix)==1 ) )
		stop("! The parameter of Z.prefix should be assigned with a prefix of time column in the phenotypic data.");
	if ( !missing("covariate.names") && length(covariate.names)>0 && !is.character(covariate.names) )
		stop("! The parameter of covariate.names should be assigned with covariate names in the phenotypic data.");
	if ( !(is.logical(fgwas.filter) && length(fgwas.filter)==1 ) )
		stop("! The parameter of fgwas.filter should be a logical value(TRUE or FALSE).");

	cat("* Phenotypic Matrix = ",  dim(phe.mat), "\n");
	cat("* SNP Matrix = ",  dim(snp.mat), "\n");

	show_fgwas_parameters( curve, covariance, Y.prefix, Z.prefix, covariate.names, fgwas.filter ) ;

	if (missing(options)) 
		options <- get_default_options()
	else	
	{
		options0 <- get_default_options();
        options0[names(options)] <- options;
        options <- options0;
    }
	
	cat( "Checking the optional items......\n");
	show_options( options);

	if( class(phe.mat)=="data.frame" )
 	{
 		cat("Phenotypic data frame is converted to the matrix class.\n");  
 		phe.colnames <- colnames(phe.mat); 
 		phe.rownames <- rownames(phe.mat); 
 		phe.mat <- matrix(as.numeric(as.matrix(phe.mat, rownames.force=NA)), ncol=NCOL(phe.mat))
 		colnames(phe.mat) <- phe.colnames;
 		rownames(phe.mat) <- phe.rownames;
	}


	r.snpmat <- list();
	r.snpmat$options <- options;
	r.snpmat$params <- list( 
				file.phe       = file.phe, 
				file.snp       = file.snp, 
				Y.prefix       = Y.prefix, 
				Z.prefix       = Z.prefix, 
				covariate.names= covariate.names, 
				fgwas.filter   = fgwas.filter);

	r.snpmat$phe.mat <- phe.mat;

	r.est  <- fg.estimate( r.snpmat$phe.mat, curve, covariance )
	if( r.est$error )
		stop(r.est$err.info)
	else
	{
		if(is.na(curve)) curve <- r.est$curve;	
		if(is.na(covariance)) covariance <- r.est$covariance;

		r.snpmat$curve <- curve;
		r.snpmat$covariance <- covariance;
		r.snpmat$est.values <- r.est$est.values;
	}


	if(fgwas.filter)
	{
		r.filter <- snpmat_fgwas_filter( phe.mat, snp.mat, Y.prefix, Z.prefix, covariate.names, options$nParallel.cpu, options$fgwas.cutoff, "GLS");
		if( r.filter$error )
			stop(r.filter$err.info);
		
		r.snpmat$filter <- r.filter$r.fgwas;
		snp.mat <- r.filter$snp.mat;
	}
	

	subset_op <- function(snpmat, sub.idx)
	{
		return( snpmat[sub.idx,,drop=F] );
	}

	r.fgwas <- snpmat_parallel(
			NROW(r.snpmat$snp.mat),
			subset_op,
			r.snpmat$snp.mat,
			phe.mat,
			Y.prefix, 
			Z.prefix,
			covariate.names,
			options$debug,
			options$nParallel.cpu);

	if(!is.null(r.fgwas) && !is.na(r.fgwas) )
	{
		r.snpmat$fgwas <- r.fgwas;
		return(r.snpmat);		   
	}
	else
	{
		cat("! No results\n");

		r.snpmat$fgwas <- "try-error";
		return(r.snpmat);		   
	}
}

summary.fGWAS.ret<-function(object, ...)
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

	class(r.sum.ret) <- "sum.fGWAS.ret";
	
	r.sum.ret
}

print.sum.fGWAS.ret<-function(x, ...)
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


show_fgwas_parameters<-function( curve, covariance, Y.prefix, Z.prefix, covariate.names, fgwas.filter ) 
{
	cat( "* Response Variable =",   Y.prefix, "\n");
	cat( "* Time Variable =",   Z.prefix, "\n");
	cat( "* Covariate Columns =",  covariate.names, "\n");
	cat( "* fGWAS Filter Used =",  ifelse( fgwas.filter, "Yes", "No"), "\n");
}

get_sig_fgwas_snp <- function( r.gls )
{
	if( is.null(r.gls$varsel_add) && is.null(r.gls$varsel_dom) ) return( NULL );

	idx.sig.add <- c();
	if( !is.null(r.gls$varsel_add) )
		idx.sig.add <- which( rowSums(r.gls$varsel_add[,c(3:6),drop=F]) > 0 );

	idx.sig.dom <- c();
	if( !is.null(r.gls$varsel_dom) )
		idx.sig.dom <- which( rowSums(r.gls$varsel_dom[,c(3:6),drop=F]) > 0 );

	idx.sig <- unique(c(idx.sig.dom, idx.sig.add));
	if (length( idx.sig )==0) return(NULL);
	
	return(idx.sig);
}


plot.fGWAS.ret<-function( x, y=NULL, ... , fig.prefix=NULL )
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

read_simple_fgwas_data <- function( file.phe, file.snp, bImputed=TRUE )
{
	tb.phe <- read.csv(file.phe, header=T);
	rownames(tb.phe) <- tb.phe[,1];
	tb.phe <- tb.phe[,-1, drop=F];
	
	tb.snp <- read.csv(file.snp, header=T);
	
	cat("Checking data files......\n");
	cat("* Individuals:", NCOL(tb.snp)-2, "\n");
	cat("* SNPs:", NROW(tb.snp), "\n");
	
	if(bImputed) tb.snp <- impute_simple_snp(tb.snp);

	return(list(phe.mat=tb.phe, snp.mat=tb.snp));
}
