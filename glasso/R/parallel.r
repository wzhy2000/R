snpmat_parallel<-function( n.snp,
			f_subset_op,
			snp.mat,
			phe.mat,
  		   	Y.name, 
  		   	Z.name,
  		   	covar.names, 
  		   	refit,
  		   	add.used,
  		   	dom.used,
  		   	op.piecewise.ratio,
  		   	op.nMcmcIter,
		   	op.fBurnInRound,
		   	op.fRhoTuning,
	           	op.fQval.add,
	           	op.fQval.dom,
			op.debug,
			op.cpu,
			lasso.method)
{
	cat( "Genetic Effect Analysis by BLASSO/GLASSO method......\n");

	n.inv <- NROW(phe.mat);
	
	if( op.piecewise.ratio == 0 || n.snp < n.inv * op.piecewise.ratio )
	{
		
		cat("* Final LASSO calling.\n");
		r.xls <- snpmat_call( snp.mat,
				phe.mat,
				Y.name, 
				Z.name,
				covar.names,
				refit,
				add.used,
				dom.used,
				op.nMcmcIter,
				op.fBurnInRound,
				op.fRhoTuning,
				op.fQval.add,
				op.fQval.dom,
				op.debug,
				lasso.method);
		return(r.xls);	
	}
	
	real.task <- ifelse( op.cpu>0, op.cpu*2, 1);
	snp.sect0 <- seq(1, n.snp, n.inv * op.piecewise.ratio);
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
				if ( length(idx.snp) < n.inv * op.piecewise.ratio )
				{
					n.add <- n.inv * op.piecewise.ratio - length(idx.snp);
					idx.snp <- c( sample(1:snp.sect0[ i + k - 1])[1:n.add], idx.snp );
					if( length( which( is.na( idx.snp ) ) ) >0 ) 
						idx.snp <- idx.snp[ -which(is.na(idx.snp)) ];
				}
				
				snpmat.list[[k]] <- f_subset_op(snp.mat, idx.snp );
			}
		}
	
		r.cluster.init <- snpmat_parallel_list( phe.mat, 
					snpmat.list, 
					Y.name, 
					Z.name,
					covar.names, 
					refit=F,
					add.used,
					dom.used,
					op.nMcmcIter,
					op.fBurnInRound,
					op.fRhoTuning,
					op.fQval.add,
					op.fQval.dom,
					op.debug,
					op.cpu, 
					lasso.method);

		for(k in 1:length(r.cluster.init))
		{
			if(lasso.method=="BLS")
				sig <- get_sig_bls_snp( r.cluster.init[[k]], snpmat.list[[k]] )	
			else
				sig <- get_sig_gls_snp( r.cluster.init[[k]], snpmat.list[[k]] );		
			
			if(!is.null(sig)) varsel.snp <- rbind( varsel.snp, sig$snp.mat);		
		}
	}
	
	cat("*", NROW(varsel.snp), "SNPs are selected in the first run.\n");
	
	R <- 1;
	while( NROW(varsel.snp) > n.inv * op.piecewise.ratio )
	{
		n.snp0 <- NROW(varsel.snp);
		snp.sect0 <- seq(1, n.snp0, n.inv * op.piecewise.ratio);
		snp.sect1 <- c( snp.sect0[-1]-1, n.snp0 );
		sect.list <- seq(1, length(snp.sect0), real.task );

		varsel.snp0 <- c();				

		for(i in sect.list )
		{
			snpmat.list <- list();
			for (k in 1:length(snp.sect0))
			{
				if( i + k - 1 > length(snp.sect0) )
					next;
				
				idx.snp <- snp.sect0[ i + k - 1]:snp.sect1[i + k -1];
				if ( length(idx.snp) < n.inv * op.piecewise.ratio )
				{
					n.add <- n.inv * op.piecewise.ratio - length(idx.snp);
					idx.snp <- c( sample(1:snp.sect0[ i + k - 1])[1:n.add], idx.snp );
					if( length( which( is.na( idx.snp ) ) ) >0 ) 
						idx.snp <- idx.snp[ -which(is.na(idx.snp)) ];
				}
				
				snpmat.list[[k]] <- varsel.snp[ idx.snp, ]
			}
			
			r.cluster <- snpmat_parallel_list( phe.mat, 
						snpmat.list, 
						Y.name, 
						Z.name,
						covar.names, 
						refit=F,
						add.used,
						dom.used,
						op.nMcmcIter,
						op.fBurnInRound,
						op.fRhoTuning,
						op.fQval.add,
						op.fQval.dom,
						op.debug,
						op.cpu,
						lasso.method);

			for(k in 1:length(r.cluster))
			{
				if(lasso.method=="BLS")
					sig <- get_sig_bls_snp( r.cluster[[k]], snpmat.list[[k]] )	
				else
					sig <- get_sig_gls_snp( r.cluster[[k]], snpmat.list[[k]] );		

				if(!is.null(sig)) varsel.snp <- rbind( varsel.snp, sig$snp.mat);		
			}

		}

		varsel.snp <- varsel.snp0; 
		
		R <- R + 1;
		cat("*", NROW(varsel.snp), "SNPs are selected in ", R, "(th) run.\n")
	}
	
	cat("* Final LASSO calling.\n");
	r.xls <- snpmat_call( 
			varsel.snp,
			phe.mat,
			Y.name, 
			Z.name,
			covar.names,
			refit=T,
			add.used,
			dom.used,
			op.nMcmcIter,
			op.fBurnInRound,
			op.fRhoTuning,
			op.fQval.add,
			op.fQval.dom,
			op.debug, 
			lasso.method);
	
	r <- wrap_xls_parallel(r.xls, r.cluster.init);
	
	return(r);	
}

snpmat_parallel_list<-function( phe.mat,
				snpmat.list,
				Y.name, 
				Z.name,
				covar.names,
				refit=F,
				add.used,
				dom.used,
				op.nMcmcIter,
				op.fBurnInRound,
				op.fRhoTuning,
				op.fQval.add,
				op.fQval.dom,
				op.debug,
				op.ncpu,
				lasso.method)
				
{
	cpu.fun<-function( sect )
	{
		library(glasso);
		
		r.xls.i <- snpmat_call(
				as.matrix( snpmat.list[[ sect ]]),
				phe.mat,
				Y.name, 
				Z.name,
				covar.names,
				refit=F,
				add.used,
				dom.used,
				op.nMcmcIter,
				op.fBurnInRound,
				op.fRhoTuning,
				op.fQval.add,
				op.fQval.dom,
				op.debug,
				lasso.method);

		return(r.xls.i);
	}

	r.cluster <- list();
	
	if( op.ncpu>1 && require(snowfall) )
	{
		cat("Starting parallel computing, snowfall/snow......\n"); 
		sfInit(parallel = TRUE, cpus = op.ncpu, type = "SOCK")
		
		sfExport("phe.mat", "Y.name", "Z.name", "covar.names", 
				"snpmat.list" ,
				"refit",
				"add.used",
				"dom.used",
				"op.nMcmcIter",
				"op.fBurnInRound",
				"op.fRhoTuning",
				"op.fQval.add",
				"op.fQval.dom",
				"op.debug",
				"lasso.method");

		r.cluster <- sfClusterApplyLB( 1:length(snpmat.list), cpu.fun);

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

snpmat_call<-function(  snp.mat,
			phe.mat,
			Y.name, 
  		   	Z.name,
  		   	covar.names, 
  		   	refit,
  		   	add.used,
  		   	dom.used,
  		   	op.nMcmcIter,
		   	op.fBurnInRound,
		   	op.fRhoTuning,
	           	op.fQval.add,
	           	op.fQval.dom,
			op.debug,
			lasso.method)
{			
	if(lasso.method=="BLS")
	{
		r <- .Call("bls_snpmat", 
			as.matrix( phe.mat ),
			as.matrix( snp.mat*1.0 ),
			Y.name, 
			paste(covar.names, collapse=","), 
			refit,
			add.used,
			dom.used,
			op.nMcmcIter,
			op.fBurnInRound,
			op.fRhoTuning,
			op.fQval.add,
			op.fQval.dom,
			ifelse( op.debug, 3, 1) );
	}
	else
	{
		r <- .Call("gls_snpmat", 
			as.matrix( phe.mat ),
  		   	as.matrix( snp.mat*1.0  ),
  		   	Y.name, 
  		   	Z.name, 
	   		paste(covar.names, collapse=","), 
  		   	refit,
  		   	add.used,
  		   	dom.used,
  		   	op.nMcmcIter,
		   	op.fBurnInRound,
		   	op.fRhoTuning,
	           	op.fQval.add,
	           	op.fQval.dom,
			ifelse(op.debug, 3, 1));
	}
	
	return(r);
}	

wrap_xls_parallel<-function(r.xls, r.cluster)
{
	varsel <- c();
	varsel_add <- c();
	varsel_dom <- c();
	varsel_Qbest <- c();
	varsel_PSRF <- c();
	
	get_nodup<-function(r.org, r.new)
	{
		row.inter <- intersect( rownames(r.new), rownames(r.org) );
		if( length(row.inter)>0 )
		{
			dup.pos <- match( row.inter, rownames(r.new) );
			return( r.new[ -dup.pos,,drop=F] );
		}
		else
			return( r.new );
	}
	
	for(i in 1:length(r.cluster) )
	{
		if ( !is.null(r.cluster[[i]]$varsel) )
		{
			r.nodup <- get_nodup( varsel, r.cluster[[i]]$varsel);
			varsel <- rbind( varsel,  r.nodup );
		}
		if ( !is.null(r.cluster[[i]]$varsel_add) )
		{
			r.nodup <- get_nodup( varsel_add, r.cluster[[i]]$varsel_add);
			varsel_add <- rbind( varsel_add, r.nodup );
		}
		if ( !is.null(r.cluster[[i]]$varsel_dom) )
		{
			r.nodup <- get_nodup( varsel_dom, r.cluster[[i]]$varsel_dom);
			varsel_dom <- rbind( varsel_dom, r.nodup );
		}
		if ( !is.null(r.cluster[[i]]$varsel_Qbest) )
		{
			r.nodup <- get_nodup(varsel_Qbest, r.cluster[[i]]$varsel_Qbest);
			varsel_Qbest <- rbind( varsel_Qbest, r.nodup );
		}
		if ( !is.null(r.cluster[[i]]$varsel_PSRF) )
		{
			r.nodup <- get_nodup( varsel_PSRF, r.cluster[[i]]$varsel_PSRF);
			varsel_PSRF <- rbind( varsel_PSRF,r.nodup );
		}
	}

	if ( !is.null(r.xls$varsel) )       r.xls$varsel <- varsel;
	if ( !is.null(r.xls$varsel_add) )   r.xls$varsel_add <- varsel_add;
	if ( !is.null(r.xls$varsel_dom) )   r.xls$varsel_dom <- varsel_dom;
	if ( !is.null(r.xls$varsel_Qbest) ) r.xls$varsel_Qbest <- varsel_Qbest;
	if ( !is.null(r.xls$varsel_PSRF) )  r.xls$varsel_PSRF <- varsel_PSRF;

	return(r.xls);
}
