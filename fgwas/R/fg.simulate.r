fg.simulate<-function( file.prefix, curve, covariance, n.obs, n.snp, time.points, par0=NA, par1=NA, par2=NA, par.covar=NA, par.X=NA,
		phe.missing=0.03, snp.missing=0.03, sig.pos=NA, plink.format=FALSE )
{
	library(mvtnorm);

	check_n_obs( n.obs );
	check_n_snp( n.snp );
	check_no_curve( no.curve );
	check_no_covar( no.covar);
	check_time_points( time.points );

	if(!missing(par.qq)) check_par_qq( par.qq );
	if(!missing(par.Qq)) check_par_Qq( par.Qq );
	if(!missing(par.QQ)) check_par_QQ( par.QQ );
	if(!missing(par.covar)) check_par_covar( par.covar );
	if(!missing(par.X)) check_par_X( par.X );

	if ( !missing(file.prefix)) check_file_data( file.prefix );
	if ( !missing(sig.pos)) check_sig_pos( sig.pos );
	
	if (is.null(sig.pos))
	{
		sig.pos <- round( runif(1, n.snp*0.25, n.snp*0.75) );
		cat(" * A significant SNP is randomly specified to location(", sig.pos, ")\n" );
	}
	
	if(class(curve)=="fgwas.curve")
		fg_curve <- curve
	else
		fg_curve <- fg.getCovar( curve );

	if( missing(par0) && missing(par1) && missing(par2) )
	{
		par0 <- fg_curve@par_simu[1,]; 
		par1 <- fg_curve@par_simu[2,]; 
		par2 <- fg_curve@par_simu[3,]; 
	}
	
	if(class(covariance)=="fgwas.covar")
		fg_covar <- covariance
	else
		fg_covar <- fg.getCovar( covariance );
	if(missing(par.covar)) par.covar <- fg_covar@par_simu; 
	
	fg.dat <- proc_dat_simu( n.obs, n.snp, par.X, fg_curve, fg_covar, time.points, sig.pos, snp.missing, phe.missing );
	
	fg.dat$file.phe.out <- paste( file.prefix, ".phe.csv", sep="" );
	write.csv(data.frame(ID = fg.dat$ids, fg.dat$pheX, fg.dat$pheY), file = fg.dat$file.phe.out, quote=F, row.names=F );

	if ( !plink.format )
	{
		fg.dat$file.geno.dat <- paste(file.prefix, ".geno.dat", sep="");
		write.table(data.frame(fg.dat$snp.info, fg.dat$gen), file=fg.dat$file.geno.dat, quote=F, row.names=F );

		return(list(err=err,  
			file.simple.phe = fg.dat$file.phe.out,
    	   	file.simple.snp = fg.dat$file.snp.out));
	}
	else
	{
		r <- convert_simpe_to_plink( data.frame(fg.dat$snp.info, fg.dat$gen), file.snp.out );

		return(list(err=err,  
			file.simple.phe = fg.dat$file.phe.out,
			file.plink.bed = r$file.plink.bed,
   	    	file.plink.bim = r$file.plink.bim,
   	    	file.plink.fam = r$file.plink.fam));
	}
}

# str(data)
# $ file.geno.dat: chr ".poplarMF65-201204-geno-v2.dat"
# $ file.phe.csv : chr "phe-4.csv"
# $ n.obs  : int 66
# $ n.snp   : int 156362
# $ times   : int [1:24] 1 2 3 4 5 6 7 8 9 10 ...
# $ ids     : Factor w/ 66 levels "NO.45","NO.6586",..: 66 1 3 63 56 51 43 39 33 31 ...
# $ snp.info:'data.frame':       156362 obs. of  4 variables:
#  ..$ scaffoldId: int [1:156362] 1 1 1 1 1 1 1 1 1 1 ...
#  ..$ Loci      : int [1:156362] 72041 72124 72263 72359 74654 74813 75205 77129 77223 78064 ...
#  ..$ RefBase   : Factor w/ 4 levels "A","C","G","T": 2 1 1 3 1 3 3 1 4 1 ...
#  ..$ AltBase   : Factor w/ 4 levels "A","C","G","T": 4 2 3 1 2 1 4 4 3 4 ...
# $ gen     :'data.frame':       156362 obs. of  66 variables:
#  ..$ NO.69  : Factor w/ 16 levels "AA","AC","AG",..: 8 2 3 9 2 1 12 4 15 4 ...
#  ..$ NO.45  : Factor w/ 16 levels "AA","AC","AG",..: 6 2 1 9 2 9 12 4 15 4 ...
#  ..$ NO.6592: Factor w/ 25 levels "..",".A",".C",..: 13 13 9 7 7 17 19 10 24 7 ...
# $ phe     :'data.frame':       66 obs. of  24 variables:
#  ..$ X1 : num [1:66] 1 0.9 4.4 8.4 4.7 2.9 3.6 3 2.2 2.4 ...
#  ..$ X2 : num [1:66] 2.7 2.5 6.5 10.5 7.7 5.5 5.6 5 3.7 4.7 ...
#  ..$ X3 : num [1:66] 4.6 4.2 8.1 12 10.1 7.9 7.2 6.7 5.1 6.8 ...
#  ..$ X4 : num [1:66] 6.5 5.9 9.5 13.2 12.1 10 8.7 8.3 6.3 8.7 ...
#  ..$ X5 : num [1:66] 8.4 7.5 10.8 14.2 13.8 11.9 10 9.7 7.4 10.4 ...
#  ..$ X6 : num [1:66] 10.1 9 11.9 15.1 15.3 13.5 11.2 11 8.5 12 ...
#  ..$ X7 : num [1:66] 11.7 10.4 12.9 15.9 16.6 15 12.3 12.2 9.5 13.4 ...
#  ..$ X8 : num [1:66] 13.2 11.6 13.8 16.6 17.8 16.4 13.4 13.3 10.5 14.8 ...
#  ..$ X9 : num [1:66] 14.6 12.7 14.7 17.2 18.8 17.6 14.5 14.3 11.5 16 ...
# $ error   : logi FALSE


proc_dat_simu<-function( n.obs, n.snp, par.X, f.curve, f.covar, times, sig.idx, snp.missing, phe.missing )
{
	dat <- list(n.obs=n.obs, n.snp=n.snp,times=times);

	#generate SNPs
	dat$gen <- proc_simu_geno( n.obs, n.snp, sig.idx, prob.miss = snp.missing )
	dat$snp.info <- data.frame(chr=1, pos=1:n.snp, RefBase="A", AltBase="B");
	dat$ids <- paste("N", 1:n.obs, sep="");

	pheX <- matrix( rep(1, n.obs), ncol=1);
	if(length(par.X)>0)
	{
		for(i in 1:length(par.X) )
		{
			if(i==1)
				pheX <- round( runif(n.obs, 1, 2) )
			else			
				pheX <- cbind( pheX, runif(n.obs, -1, 1) ); 
		}
		colnames(pheX) <- paste("X", 1:length(par.X), sep="_");
		
	}
	else
		par.X<-0;

	pheY <- array( 0, dim=c(n.obs, length(times))); 
	
	colnames( pheY ) <- paste("Y", times, sep="_");
	rownames( pheY ) <- dat$ids;
	
	#generate traits
	sim.covar<-  f.covar$func( f.covar$simu_par, times )
	sim.mu   <-  f.curve$func( f.curve$simu_par$QQ, times );
	sim.mu   <-  rbind(sim.mu, f.curve$func( f.curve$simu_par$Qq, times ) );
	sim.mu   <-  rbind(sim.mu, f.curve$func( f.curve$simu_par$qq, times ) );
	
	if(is.null(f.covar$simu_par))
		sim.mu   <-  rbind(sim.mu, f.curve$func( f.curve$qq.par, times ) );

	d.gen <- dat$gen[sig.idx,];
	d.g <- array( 9, n.obs );
	if(1)
	{
		d.gen2 <- as.character(unlist(d.gen));
		snpB <- as.character(unlist(dat$snp.info[sig.idx,4]));
		snpA <- as.character(unlist(dat$snp.info[sig.idx,3]));
		
		QQ2<- paste(snpB, snpB, sep=""); 
		qq0<- paste(snpA, snpA, sep="");
		Qq1<- c(paste(snpA, snpB, sep=""), paste( snpB, snpA, sep="") ) ;
	
		d.g[which(d.gen2==QQ2)]<-2;
		d.g[which(d.gen2==qq0)]<-0;
		d.g[which(d.gen2==Qq1[1])]<-1;
		d.g[which(d.gen2==Qq1[2])]<-1;
	}

	for (i in 1:n.obs)
	{
		 if (d.g[i]==9) d.g[i] <- round(runif(1, 0, 2));
		 
		 y <- rmvnorm(1, sim.mu[ d.g[i] + 1, ], sim.covar );
		 pheY[i, ] <- y + sum( pheX[i,] * par.X );
	}
	
	dat$pheY    <- pheY;
	dat$pheX    <- pheX;
	dat$f.curve <- f.curve;
	dat$f.covar <- f.covar;
	dat$error   <- F;
	
	cat("Data simulation is done![Sig=", sig.idx, "]\n");
	
	return(dat);
}

#--------------------------------------------------------------
# private: fin.generate_bc_marker;
#
# genarate N Backcross Markers from marker disttance (cM): dist.
#
# input
#     dist : the vector for the marker distances
#   samp_N : the sample size
#--------------------------------------------------------------
fin.generate_bc_marker<-function( n.obs, dist )
{
	if (dist[1] != 0)
		cm=c(0, dist)/100
	else
		cm=dist/100;

	n  <- length(cm);
	rs <- 1/2*( exp(2*cm) - exp(-2*cm) ) / (exp(2*cm)+exp(-2*cm));
	mk <- array( 0, dim=c( n.obs, n ) );

	for (j in 1:n.obs)
		mk[j,1] <- ( runif(1)>0.5 );

	for (i in 2:n)
		for ( j in 1:n.obs )
		{
			if (mk[j,i-1]==1)
				mk[j,i] <- ( runif(1)>rs[i] )
			else
				mk[j,i] <- ( runif(1)<rs[i] );
    	}
    
	return(mk);
}

proc_simu_geno<-function( n.obs, n.snp, sig.idx, prob.miss=0.03 )
{
	dist <- cumsum( runif(n.snp, 0.05, 0.12 ) );
	snp1 <- t( fin.generate_bc_marker( n.obs,  dist) );
	snp2 <- t( fin.generate_bc_marker( n.obs,  dist) );
	for(i in 1:n.snp)
	{
		nmiss <- runif( 1, 0, n.obs*prob.miss );
		if (nmiss>=1)
			snp1[ i, sample(n.obs)[1:round(nmiss)] ] <- 9;
	}

	for(i in 1:n.snp)
	{
		nmiss <- runif( 1, 0, n.obs*prob.miss );
		if (nmiss>=1)
			snp2[ i, sample(n.obs)[1:round(nmiss)] ] <- 9;
	}

	cors <- c();
	snpx <- snp1+snp2;
	for(i in 1:n.snp)
		cors <- c(cors, cor(snpx[i,], snpx[sig.idx, ]) );
	
	cor.high <- which( abs(cors)>0.75 );
	if ( length( cor.high )>=2)
	{
		for( i in 1:length(cor.high) )
			if( cor.high[i]!=sig.idx)
			{
				snp1[ cor.high[i], ] <- snp1[ cor.high[i] -1, ];
				snp2[ cor.high[i], ] <- snp2[ cor.high[i] -1, ];
			}
	}

	snp1.s <- c(as.character(snp1))
	snp2.s <- c(as.character(snp2))

	if( length( which(snp1.s=="0") ) >0 ) snp1.s[ which(snp1.s=="0") ] <- "A";
	if( length( which(snp1.s=="1") ) >0 ) snp1.s[ which(snp1.s=="1") ] <- "B";
	if( length( which(snp1.s=="9") ) >0 ) snp1.s[ which(snp1.s=="9") ] <- ".";

	if( length( which(snp2.s=="0") ) >0 ) snp2.s[ which(snp2.s=="0") ] <- "A";
	if( length( which(snp2.s=="1") ) >0 ) snp2.s[ which(snp2.s=="1") ] <- "B";
	if( length( which(snp2.s=="9") ) >0 ) snp2.s[ which(snp2.s=="9") ] <- ".";
	
	snp.s <- paste(snp1.s, snp2.s, sep="");
	gen<-array(snp.s, dim=c(n.snp, n.obs))
	colnames(gen)<-paste("N",1:n.obs, sep="");
	rownames(gen)<-paste("S",1:n.snp, sep="");
	
	return(gen);
}
