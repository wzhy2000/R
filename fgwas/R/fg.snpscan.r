fg.snpscan<-function( fg.dat, fgwas.filter = F, options = NULL )
{
	check_FG_DAT( fg.dat ); 
	
	options0 <- get_default_options();
	if (! missing(options)) options0[names(options)] <- options;
       	options <- options0;

	r.filter <- list();
	r <- list();

	subset_simple_op <- function(snpmat, sub.idx)
	{
		return( snpmat[sub.idx,,drop=F] );
	}

	subset_plink_op <- function(snpmat, sub.idx)
	{
		return( snpmat[sub.idx,,drop=F] );
	}

	if( fgwas.filter)
	{
		r.filter <- snpmat_fgwas_filter( 
				fg.dat$phe.mat, 
				fg.dat$snp.mat, 
				"Y", 
				"T", 
				colnames(fg.dat$pheX),
				options$nParallel.cpu, 
				options$fgwas.cutoff );
				
		if( r.filter$error )
			stop(r.filter$err.info);
	
		r <- mle_parallel(
			NROW( r.filter$snp.mat ),
			subset_simple_op,
			r.filter$snp.mat,
			fg.dat$pheY,
			fg.dat$pheX,
			fg.dat$pheT,
			fg.dat$phe.est,
			options$nParallel.cpu,
			options$debug );

	}
	else
	{
		r <- mle_parallel(
			NROW( fg.dat$snp.mat ),
			subset_simple_op,
			fg.dat$snp.mat,
			fg.dat$pheY,
			fg.dat$pheX,
			fg.dat$pheT,
			fg.dat$phe.est,
			options$nParallel.cpu,
			options$debug );
	}
	
	r.scan <- list( 
			filter = r.filter,
			snp.range = fg.dat$snp.range, 
			phe.est = fg.dat$phe.est, 
			results = r, 
			top = proc_top_select( r, 10 ) ,
			options = options );
			
	class(r.scan) <- "FGWAS.SCAN";		

	return(r.scan);
}

fg.detect.sigsnp<-function( fg.dat, fg.scan, fg.perm=NULL )
{
	check_FG_DAT( fg.dat );
	check_FG_SCAN( fg.scan );
	
	if ( !missing(fg.perm) ) 
	{
		check_FG_PERM( fg.perm );
	
		sig.05 <- NA;
		sig.01 <- NA;
		if ( !is.null(fg.perm$pcut.05) || !is.null(fg.perm$pcut.01) )
		{
			sig <- proc_sig_select(fg.scan$full, fg.perm);
			sig.05 <- sig$sig.05;
			sig.01 <- sig$sig.01;
		}
	}
	else
	{
		#todo...
	}
	
	return(list(error=F, sig.05=sig.05, sig.01=sig.01));
}

mle_parallel<-function( n.snp,
			f_subset_op,
			snp.mat,
			pheY,
			pheX,
			pheT,
			phe.est,
			op.cpu,
			op.debug )
{
	cat( "Genetic Effect Analysis by fGWAS method......\n");

	n.inv <- NROW(phe.mat);
	n.snpratio <- 10;
	real.task <- ifelse( op.cpu>0, op.cpu*2, 1);
	snp.sect0 <- seq( 1, n.snp, n.inv * n.snpratio );
	snp.sect1 <- c( snp.sect0[-1]-1, n.snp );
	sect.list <- seq( 1, length(snp.sect0), real.task );
	
	varsel.snp <- c();				
	for(i in sect.list )
	{
		snpmat.list <- list();
		for(k in 1:real.task )
		{
			if( i + k - 1 <= length(snp.sect0) )
			{
				idx.snp <- c(snp.sect0[ i + k - 1]:snp.sect1[ i + k - 1 ]);
				if ( length(idx.snp) < n.inv * n.snpratio )
				{
					n.add <- n.inv * n.snpratio - length(idx.snp);
					idx.snp <- c( sample(1:snp.sect0[ i + k - 1])[1:n.add], idx.snp );
					if( length( which( is.na( idx.snp ) ) ) >0 ) 
						idx.snp <- idx.snp[ -which(is.na(idx.snp)) ];
				}
				
				snpmat.list[[k]] <- f_subset_op(snp.mat, idx.snp );
			}
		}
	
		r.cluster.init <- mle_parallel_list( 
					snpmat.list, 
					pheY,
					pheX,
					pheT,
					phe.est,
					op.cpu, 
					op.debug );

		varsel.snp <- rbind( varsel.snp, sig$snp.mat);		
	}
	
	cat("*", NROW(varsel.snp), "SNPs are selected in the first run.\n");
	
	return(r);	
}

mle_parallel_list<-function( snpmat.list,
				pheY,
				pheX,
				pheT,
				phe.est,
				op.ncpu,
				op.debug )
				
{
	cpu.fun<-function( sect )
	{
		library(glasso);
		
		r.xls.i <- mle_call(
				as.matrix( snpmat.list[[ sect ]]),
				pheY,
				pheX,
				pheT,
				phe.est,
				op.debug );

		return(r.xls.i);
	}

	r.cluster <- list();
	
	if( op.ncpu>1 && require(snowfall) )
	{
		cat("Starting parallel computing, snowfall/snow......\n"); 
		sfInit(parallel = TRUE, cpus = op.ncpu, type = "SOCK")
		
		sfExport("phe.mat", "snpmat.list", "est.H0", "op.debug" );

		r.cluster <- sfClusterApplyLB( 1:length(snpmat.list), cpu.fun );

		sfStop();

		cat("Stopping parallel computing......\n");
	}		
	else
	{
		cat("Starting piecewise analysis......\n");
		for(i in 1:length(snpmat.list))
			r.cluster[[i]] <- cpu.fun( i );
	}
	
	return( r.cluster );
}

mle_call<-function(snp.mat,
		  pheY,
		  pheX,
		  pheT,
		  phe.est,
		  op.debug )
{			
	#proc_cache_reset( fg.dat$file.phe.csv, cache.reset );
	
	rlist   <- list();
	for ( k in 1:NOR(snp.mat) )
	{
		e <- try( r0 <- proc_mle( k, snp.mat[k,], pheY, pheX, pheT, phe.est, op.debug ), FALSE );
		if( class(e) == "try-error" )
			return(list(error=T, err.info="ERR-MSG"));
		
		rlist[[k]] <- r0;
	}

	r.scan   <- lapply( rlist, rbind );

	#proc_cache_merge();

	return(r.scan);
}	

summary.FGWAS.SCAN <- function(fg.scan)
{
#todo
}