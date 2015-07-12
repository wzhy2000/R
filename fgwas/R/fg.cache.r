proc_cache_reset<-function( phe.csv, cache.reset=T )
{
	phe.name <- remove_extname(phe.csv);
	rdata.cache <- paste(phe.name, ".cache.rdata", sep="");

	fg.reset <- cache.reset;
	if( !fg.reset)
	{
		x<-try( load(rdata.cache), TRUE );
		if( class(x)=="try-error" ) fg.reset <- TRUE;
	}

	if (fg.reset)
	{
		fg.cache <<- list( 
			mu     = c(), 
			umiss  = c(), 
			rdata  = rdata.cache, 
			perm.sample  = NULL,
			perm.unsample  = NULL,
			found  = 0 );
		save( fg.cache, file=fg.cache$rdata, envir = .GlobalEnv );
	}
	
	return;
}
proc_cache_find<-function( umiss )
{
	if (length(which(umiss!=9))>0 )	umiss[which(umiss!=9)]<-1;
	if (length(which(umiss==9))>0 )	umiss[which(umiss==9)]<-0;
	if (!is.null(fg.cache$perm.sample))
		umiss <- umiss[fg.cache$perm.unsample];
	
	str <- paste(umiss, sep="", collapse="");
	find <- which(fg.cache$umiss==str);

	if (length(find)==0)
		return(NULL)
	else
	{
		if (fg.sys$log>=LOG_DEBUG) 
			cat("CacheHit[", fg.cache$found, "]", fg.cache$umiss[ find[1] ], "\n");		

		fg.cache$found <<- fg.cache$found+1;
		n.len <- dim(fg.cache$mu)[2];
		return( fg.cache$mu[ find[1], 2:n.len] )
	}
}

proc_cache_save<-function( umiss, par.mu, value )
{
	if (length(which(umiss!=9))>0 )	umiss[which(umiss!=9)]<-1;
	if (length(which(umiss==9))>0 )	umiss[which(umiss==9)]<-0;
	if (!is.null( fg.cache$perm.sample ))
		umiss <- umiss[ fg.cache$perm.unsample ];

	str <- paste(umiss, sep="", collapse="")

	find <- which(fg.cache$umiss==str);

	if (length(find)==0)
	{
		fg.cache$umiss <<- c( fg.cache$umiss, str);
		fg.cache$mu <<- rbind( fg.cache$mu, c(value, par.mu) );
	}
	else
	{
		if (fg.cache$mu[ find[1],1] >value)
			fg.cache$mu[ find[1],] <<- c(value, par.mu)
	}

	#cat("CAHCE ITEMS", length(fg.cache$umiss), "\n");	
}

proc_cache_merge<-function()
{
	t.cache.umiss <- fg.cache$umiss;
	t.cache.mu    <- fg.cache$mu;

	x<-try( load(fg.cache$rdata), TRUE );
	if (class(x)=="try-error" || is.null(fg.cache$umiss) )
	{
		fg.cache$umiss <<- t.cache.umiss;
		fg.cache$mu    <<- t.cache.mu;
	}
	else
	{
		if( dim(fg.cache$mu)[1] != dim( t.cache.mu ) [1] )
		{
			warning("The dimensions of cache file are different with the current data, the mergeoperation is aborted.")
			return;		
		}
		
		c.mu    <- fg.cache$mu;
		c.umiss <- fg.cache$umiss;

		for(i in 1:length(t.cache.umiss))
		{
			find <- which(fg.cache$umiss == t.cache.umiss[i] );
			if (length(find)==0)
			{
				c.umiss <- c( c.umiss, t.cache.umiss[i]);
				c.mu    <- rbind( c.mu, t.cache.mu[i,] );
			}
			else
			{
				if ( c.mu[ find[1],1] > t.cache.mu[i,1] )
					c.mu[ find[1],] <- t.cache.mu[i,];
			}
		}

		fg.cache$umiss <<- c.umiss;
		fg.cache$mu    <<- c.mu;


	}

	if( fg.sys$log >= LOG_DEBUG ) 
		cat("SAVE CAHCE ITEMS", length(fg.cache$umiss), "\n");	
	
	save( fg.cache, file=fg.cache$rdata, envir = .GlobalEnv );
	assign("fg.cache", fg.cache, envir = .GlobalEnv)
}

.onAttach <- function(library, pkg) 
{
	fg.cache <<- list( 
		mu     = c(), 
		umiss  = c(), 
		rdata  = "fg.cache.rdata", 
		perm.sample  = NULL,
		perm.unsample  = NULL,
		found  = 0 );
}
