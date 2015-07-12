fg_permu_core<-function( fg.dat, n.loop=5 )
{
	check_FG_DAT( fg.dat );
	if ( !missing(n.loop) )  check_n_loop( n.loop );
	
	perm.sample   <- sample( 1 : fg.dat$n.obs);
	perm.unsample <- rep(0, fg.dat$n.obs);
	for(i in 1:length( perm.sample ) )
		perm.unsample[i] <- which( perm.sample==i )[1];

	proc_cache_reset( fg.dat$file.phe.csv );

	fg.sys$perm.sample    <<-  perm.sample;
	fg.sys$perm.unsample  <<-  perm.unsample;
	fg.cache$perm.sample  <<-  perm.sample;
	fg.cache$perm.unsample<<-  perm.unsample;

	fg.dat$phe.org <- fg.dat$phe;
	fg.dat$phe     <- fg.dat$phe[perm.sample,]
	
	r.perm <-c();
	maxLR2 <- -Inf;
	for (k in 1:fg.dat$n.snp )
	{
		r0 <- proc_mle( k, fg.dat, n.loop, b.permu=T );
		r.perm  <- rbind( r.perm,   r0 );

		if( maxLR2 < r0[4] ) maxLR2 <- r0[4];
		if( k %% 100 == 0 ) proc_chache_merge();
	}

	proc_cache_merge();
	return(list(error=F, maxLR2=maxLR2, perm.sample=perm.sample, perm.unsample=perm.unsample ) );
}

fg.permu.cutoff<-function( r.perm )
{
	pcut.05<- NA;
	if(dim(r.perm)[1]>=20)
	{
		n.05 <- round(dim(r.perm)[1]*0.05);
		if (n.05==0) n.05 <-1;
		pcut.05 <- sort(r.perm[,1], decreasing=T)[n.05]
	}
	
	pcut.01<- NA;
	if(dim(r.perm)[1]>=100)
	{
		n.01 <- round(dim(r.perm)[1]*0.01);
		if (n.01==0) n.01 <-1;
		pcut.01 <- sort(r.perm[,1], decreasing=T)[n.01]
	}

	p_cut <- c();
	if (dim(r.perm)[1]>=100)
		p_cut <- c( 0.9,  0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.09, 0.08, 0.07, 0.06, 0.05, 0.04, 0.03, 0.02, 0.01);
	if (dim(r.perm)[1]>=1000)
		p_cut <- c(p_cut, 0.009, 0.008, 0.007, 0.006, 0.005, 0.004, 0.003, 0.002, 0.001 );
	if (dim(r.perm)[1]>=10000)
		p_cut <- c(p_cut, 0.0009, 0.0008, 0.0007, 0.0006, 0.0005, 0.0004, 0.0003, 0.0002, 0.0001);
	if (dim(r.perm)[1]>=100000)
		p_cut <- c(p_cut, 0.00001);

	pv2<-NULL;
	if(length(p_cut)>0)
	{
		pv2<-array(0, dim=c(length(p_cut), 2) );
		for(i in 1:length(p_cut))
		{
			pv2.order <- round ( dim(r.perm)[1] * p_cut[i]);
			if ( pv2.order>0)
			{
				pv2[i,1] <- p_cut[i];
				pv2[i,2] <- sort(r.perm[,1], decreasing=TRUE )[ pv2.order ];
			}
		}
	}
	
	return(list(pcut.05=pcut.05, pcut.01=pcut.01, pcut.tab=pv2));
}

fg.permutation<-function( fg.dat, n.perm, file.rdata.perm=NULL, options=list( cache.reset=T, n.loop=5, debug=F) )
{
	check_FG_DAT(fg.dat);
	check_n_perm(n.perm);
	if ( !missing(file.rdata.permp) )  check_rdata_perm(file.rdata.perm);

	if (missing(options)) 
		options <- get_default_options()
	else	
	{
		options0 <- get_default_options();
        	options0[names(options)] <- options;
        	options <- options0;
    	}

	cat(" PERM START n.perm=", n.perm, "\n");
	r.perm <- c();
	for(i in 1:n.perm)
	{
		p0 <- fg_permu_core( fg.dat, n.loop );
		r.perm <- rbind(r.perm, c( p0$maxLR2, p0$perm.sample ) ) ;

		cat(" PERM [", i, "] max.LR2=", p0$maxLR2, "\n");

		if(!is.null(file.rdata.perm)) 
			save( r.perm, file = file.rdata.perm );


	}

	p_cut <- fg_permu_cutoff( r.perm );

	if (!is.null(p_cut$pcut.05) )
		cat(" PERM cutoff(pv-0.05) = ", p_cut$pcut.05, "\n");
	if (!is.null(p_cut$pcut.01) )
		cat(" PERM cutoff(pv-0.01) = ", p_cut$pcut.01, "\n");
	
	fg.perm <- list(error=F, pcut.05=p_cut$pcut.05, pcut.01=p_cut$pcut.01, pcut.tab=p_cut$pcut.tab, max.lr2=r.perm[,1], full=r.perm)
	class(fg.perm)="FGWAS.PERM";
	
	return(fg.perm);
}

fg.perm.merge<-function( fg.perm, file.rdata.perm )
{
	check_FG_PERM(fg.perm);
	check_rdata_perm(file.rdata.perm);

	if( !file.exists(file.rdata.perm) )
	{
	 	warning("RData File is not existing\n")
	 	return(list(error=T, err.info=paste("RData File is not existing(",file.rdata.perm,")\n",sep="")))
	}
	
	r <- try( load(file.rdata.perm), FALSE);
	if (class(r)=="try-error")
		return(list(error=T, err.info=paste("RData File is not existing(",file.rdata.perm,")\n",sep="")));
		
	fg.perm$full <- rbind( fg.perm$full, r.perm ) ;

	p.cut <- fg_permu_cutoff( fg.perm$full );
	
	r.newperm <- list(error = F, 
			pcut.05 = p.cut$pcut.05, 
			pcut.01 = p.cut$pcut.01, 
			pcut.tab= p.cut$pcut.tab, 
			max.lr2 = fg.perm$full[,1], 
			full    = fg.perm$full );
	class( r.newperm )="FGWAS.PERM";

	return( r.newperm );
}