gls.fgwas <- function( phe.mat, snp.mat, Y.prefix, Z.prefix, covar.names=NULL, op.cpu=0)
{
	if ( !require(nlme) )
		return(list( error=T, err.info="The package of nlme is not installed.") )

	sample.ids <- intersect( rownames(phe.mat), colnames(snp.mat)[-c(1:2)] );
	if(length(sample.ids)==0)
		return(list(error=T, err.info="No same population data in the phenotypic data and genotypic data."))

	phe.mat <- phe.mat[ sample.ids, , drop=F]
	snp.mat <- cbind( snp.mat[,c(1,2)], snp.mat[ , sample.ids, drop=F]);
	
	y.p <- grep( Y.prefix, colnames(phe.mat) ); 
	z.p <- grep( Z.prefix, colnames(phe.mat) );
	x.p <- c(); 
	if(!is.null(covar.names))
		x.p <- match( covar.names, colnames(phe.mat) );

	Y.sub <- substring( colnames(phe.mat)[y.p], nchar(Y.prefix)+1);
	Z.names <- paste( Z.prefix, Y.sub, sep="");
	z.p <- match( Z.names, colnames(phe.mat) );
	if( length( which(is.na(z.p)) )>0 )
		return(list(error=T, err.info="The Z columns are not matched with Y columns."))

	len.x <- length(x.p);
	len.y <- length(y.p);
	len.z <- length(z.p);
	if ( len.y != len.z)
		return(list(error=T, err.info="The Z columns are not matched with Y columns."))

	if(length(which(is.na(x.p)))>0)
		return(list(error=T, err.info="The covariate names are not matched with phenotypical data."))
	
	phe.mat <- cbind(phe.mat, ID = c(1:NROW(phe.mat)) );
	id.p <- NCOL(phe.mat);

	phe.gls.mat <- array(0,dim=c(0, 1 + len.x + 2 ))
	for(i in 1:len.y)
	{
		phe.tmp <- phe.mat[, c( ID=id.p, Y=y.p[i], Z=z.p[i], x.p), drop=F ];
		
		phe.tmp.col <- colnames(phe.tmp);
		phe.tmp.col[ c(1:3) ] <- c("ID", "Y", "Z");
		colnames(phe.tmp) <- phe.tmp.col;
		
		phe.tmp.missing <- which( is.na(phe.tmp$Z) | is.na(phe.tmp$Y) );
		if(length(phe.tmp.missing)>0)
			phe.tmp <- phe.tmp[-phe.tmp.missing, , drop=F ];
		phe.gls.mat <- rbind( phe.gls.mat, phe.tmp );
	}

	reg.str0 <- ifelse( len.x>0, paste("Y ~ ", paste(covar.names,collapse= "+")),  "Y ~ 1" );
	reg.str1 <- ifelse( len.x>0, paste("Y ~ ", paste(covar.names,collapse= "+"), "+ as.factor(SNP)") ,  "Y ~ 1 + as.factor(SNP)" );

	cat("* H0 =", as.character(reg.str0), "\n" );
	cat("* H1 =", as.character(reg.str1), "\n" );

	r0 <- try( gls( as.formula(reg.str0), phe.gls.mat, correlation = corAR1(form = ~ Z | ID ), method="ML" ) );	
	if(class(r0)=="try-error")
		return(list( error=T, err.info="Failed to call gls() method.") )
	
	cpu.fun<-function( sect )
	{
		library(nlme);
		
		if( (sect-1)*n.percpu+1 > NROW(snp.mat) )
			return(NULL);

		range.fr <- (sect-1)*n.percpu+1;
		range.to <- sect*n.percpu;
		range.to <- ifelse( range.to > NROW(snp.mat),  NROW(snp.mat), range.to );

		r.gls <- array( NA, dim = c( (range.to-range.fr+1), 5 ) );
		r.gls[,1]<-c(range.fr:range.to);
	
		reg.f0 <<- as.formula(reg.str0);
		reg.f1 <<- as.formula(reg.str1);

		r0 <- try( gls( reg.f0, phe.gls.mat, correlation = corAR1(form = ~ Z | ID ), method="ML" ) );	
		for(i in c(range.fr:range.to) )
		{
			snp <- unlist( snp.mat[i, c(3:NCOL(snp.mat)) ] ); 
			snp.gls <- snp[ phe.gls.mat$ID ];
			gls.mat <- cbind( phe.gls.mat, SNP=snp.gls );

			r1 <- try( gls( reg.f1, gls.mat, correlation = corAR1(form = ~ Z | ID ), method="ML" ) );    
			if(class(r1)=="try-error")
			{
				r.gls[i-range.to+1,] <- c( i, snp.mat[i,1], snp.mat[i,2], NA, NA );
				next;	
			}

			r <- try( anova( r0,r1 ) );                   
			if(class(r)=="try-error")
			{
				r.gls[i-range.to+1,] <- c( i, snp.mat[i,1], snp.mat[i,2], NA, NA );
				next;	
			}

			r.gls[(i-range.fr+1),] <- c( i, snp.mat[i,1], snp.mat[i,2], r[2,8], r[2,9]);
		}
		
		return(r.gls);
	}
	

	r.gls <- c();
	if( op.cpu>1 && require(snowfall) )
	{
		cat("Starting parallel computing, snowfall/snow......\n"); 
		sfInit(parallel = TRUE, cpus = op.cpu, type = "SOCK")

		n.percpu <- ceiling( NROW(snp.mat) / op.cpu );
		sfExport("n.percpu", "phe.gls.mat", "snp.mat", "covar.names", "reg.str0", "reg.str1" );
		
		gls.cluster <- sfClusterApplyLB( 1:op.cpu, cpu.fun);
		sfStop();

		cat("Stopping parallel computing......\n");
		for(i in 1:length(gls.cluster))
		{
			if (!is.null(gls.cluster[[i]]))
				r.gls <- rbind( r.gls, gls.cluster[[i]] );
		}
	}		
	else
	{
		cat("Starting the fGWAS estimate for each SNP......\n");
		n.percpu <- NROW(snp.mat);	
		r.gls <- cpu.fun(1);
	}

	colnames(r.gls) <- c("SNP.ID", "CHR", "POS", "L.Ratio", "pv");
	rownames(r.gls) <- rownames(snp.mat) 
	
	return(list(error=F, r=r.gls));
}

bls.fgwas <- function( phe.mat, snp.mat, Y.name, covar.names=NULL, op.cpu=0)
{
	if ( !require(nlme) )
		return(list( error=T, err.info="! The package of nlme is not installed.") )
	
	sample.ids <- intersect( rownames(phe.mat), colnames(snp.mat)[-c(1:2)] );
	if(length(sample.ids)==0)
		return(list(error=T, err.info="! No same population data in the phenotypic data and genotypic data."))

	phe.mat <- phe.mat[ sample.ids, , drop=F]
	phe.mat <- cbind(ID = c(1:NROW(phe.mat)), phe.mat);
	snp.mat <- cbind( snp.mat[,c(1,2)], snp.mat[ , sample.ids, drop=F]);
	
	len.x <- 0;
	x.p   <- c();
	if( !is.null( covar.names ) )
	{
		x.p <- match( covar.names, colnames(phe.mat) );
		len.x <- length(x.p);
	}
	
	phe.missing <- which( is.na(phe.mat[, dim(phe.mat)[2]]));
	if(length(phe.missing)>0)
		phe.mat <- phe.mat[-phe.missing, , drop=F ];

	reg.str0 <- ifelse( len.x>0, paste("Y ~ ", paste(covar.names,collapse= "+")),  "Y ~ 1" );
	reg.str1 <- ifelse( len.x>0, paste("Y ~ ", paste(covar.names,collapse= "+"), "+ as.factor(SNP)") ,  "Y ~ 1 + as.factor(SNP)" );

	cat("* H0 =", as.character(reg.str0), "\n" );
	cat("* H1 =", as.character(reg.str1), "\n" );
	
	r0 <- try( gls( as.formula(reg.str0), phe.mat, method="ML" ) );	
	if(class(r0)=="try-error")
		return(list( error=T, err.info="! Failed to call gls() method.") )
		
	cpu.fun<-function( sect )
	{
		library(nlme);
		
		if( (sect-1)*n.percpu+1 > NROW(snp.mat) )
			return(NULL);

		range.fr <- (sect-1)*n.percpu+1;
		range.to <- sect*n.percpu;
		range.to <- ifelse( range.to > NROW(snp.mat),  NROW(snp.mat), range.to );

		r.bls <- array( NA, dim = c( (range.to-range.fr+1), 5 ) );
		r.bls[,1]<-c(range.fr:range.to);
	
		reg.f0 <<- as.formula(reg.str0);
		reg.f1 <<- as.formula(reg.str1);
		
		r0 <- try( gls( reg.f0, phe.mat, method="ML" ) );	
		for(i in c(range.fr:range.to) )
		{
			snp <- unlist( snp.mat[i, c(3:NCOL(snp.mat)) ] ); 
			snp.bls <- snp[ phe.mat$ID ];
			bls.mat <- cbind( phe.mat, SNP=snp.bls );
			
			r1 <- try( gls( reg.f1, bls.mat, method="ML" ) );    
			if(class(r1)[1]=="try-error")
			{
				r.bls[i-range.to+1,] <- c( i, snp.mat[i,1], snp.mat[i,2], NA, NA );
				next;	
			}

			r <- try( anova( r0,r1 ) );                   
			if(class(r)[1]=="try-error")
			{
				r.bls[i-range.to+1,] <- c( i, snp.mat[i,1], snp.mat[i,2], NA, NA );
				next;	
			}

			r.bls[(i-range.fr+1),] <- c( i, snp.mat[i,1], snp.mat[i,2], r[2,8], r[2,9]);
		}
		
		return(r.bls);
	}
	
	r.bls <- c();
	if( op.cpu>1 && require(snowfall) )
	{
		cat("\n  Starting parallel computing, snowfall/snow......\n"); 
		sfInit(parallel = TRUE, cpus = op.cpu, type = "SOCK")

		n.percpu <- ceiling( NROW(snp.mat) / op.cpu );
		sfExport("n.percpu", "phe.mat", "snp.mat", "covar.names", "reg.str0", "reg.str1" );
		
		bls.cluster <- sfClusterApplyLB( 1:op.cpu, cpu.fun);
		sfStop();

		cat("  Stopping parallel computing......\n\n");
		for(i in 1:length(bls.cluster))
		{
			if (!is.null(bls.cluster[[i]]))
				r.bls <- rbind( r.bls, bls.cluster[[i]] );
		}
	}		
	else
	{
		cat("  Starting the fGWAS estimate for each SNP......\n");
		n.percpu <- NROW(snp.mat);	
		r.bls <- cpu.fun(1);
	}
	
	colnames(r.bls) <- c("SNP.ID", "CHR", "POS", "L.Ratio", "pv");
	rownames(r.bls) <- rownames(snp.mat) 
	
	return(list(error=F, r=r.bls));
	
}


get_sigsnp_nomulti_correction<-function( snp.mat, r.fgwas, fgwas.cutoff=0.05 )
{
	n.sample <- NCOL(snp.mat)-2;
	n.snp    <- NROW(snp.mat);
	if( n.snp <= n.sample)
		return(list(error=F, snp.mat=snp.mat));
		
	pv.sort  <- sort.int( r.fgwas[,5], decreasing=F, index.return=T);
	if ( pv.sort$x[n.sample] >= fgwas.cutoff)
	{
		sel.p <- pv.sort$ix[1:n.sample];
		sel.id <- r.fgwas[ sel.p, 1 ];
	}
	else
	{
		sel.p <- min( which(pv.sort$x > fgwas.cutoff) );
		sel.id <- r.fgwas[ pv.sort$ix[1:sel.p], 1 ];
	}
	
	snp.mat <- snp.mat[ sel.id, ,drop=F ]
	return(list(error=F, snp.mat=snp.mat));
}

fgwas_filter<-function( n.snp, n.ind, f_get_snpmat, snp.obj, phe.mat, Y.name, Z.name, covar.names, op.cpu=1, fgwas.cutoff=0.05, lasso_method="BLS" )
{
	cat( "SNP Filtering by fGWAS method......\n");

	cat("* SNP Count =", n.snp, "\n" );
	cat("* Sample Count =", n.ind, "\n" );
	cat("* p-value Threshold =", fgwas.cutoff, "\n" );
		
	snp.sect0 <- seq( 1, n.snp, 20000 );
	snp.sect1 <- c( snp.sect0[-1]-1, n.snp );
	
	r.fgwas <- c();
	snp.mat <- c();
	
	for(i in 1:length(snp.sect0))
	{
		snp.mat0 <- f_get_snpmat( snp.obj, snp.sect0[i]:snp.sect1[i] );

		cat("  Calculated SNP Range =", snp.sect0[i], snp.sect1[i], "\n" );

		r.fgwas0 <- list();
		if(lasso_method=="BLS")
			r.fgwas0 <- bls.fgwas( phe.mat, snp.mat0, Y.name, covar.names, op.cpu )
		else
			r.fgwas0 <- gls.fgwas( phe.mat, snp.mat0, Y.name, Z.name, covar.names, op.cpu )
		
		if( r.fgwas0$error )
			stop( r.fgwas0$err.info );
		
		r.fgwas <- rbind( r.fgwas, r.fgwas0$r);

		r.filter0 <- get_sigsnp_nomulti_correction( snp.mat0, r.fgwas0$r, fgwas.cutoff );
		if( r.filter0$error )
			stop( r.filter0$err.info );

		snp.mat <- rbind( snp.mat, r.filter0$snp.mat);

	}
	
	f.per <- round( NROW(snp.mat)/n.snp*100, digits=2);
	cat("*", NROW(snp.mat), "SNPs(", f.per ,"%) are left after fGWAS filtering.\n" );

	return(list(error=F, snp.mat=snp.mat, r.fgwas=r.fgwas) );
}


plink_fgwas_filter<-function( pd, Y.name, Z.name, covar.names, op.cpu=1, fgwas.cutoff=0.05, lasso_method="BLS")
{
	n.snp <- NCOL( pd$snp.mat$genotypes );
	n.ind <- NROW( pd$snp.mat$genotypes );
		
	get_sub_snpmat<- function(pd.obj, idx.snp)
	{
		snp.sub <- get_sub_snp( pd$snp.mat, idx.snp );
		# Append Chr and pos. information to SNP.MAT
		snp.mat <- cbind( snp.sub$info[,c(2,3)], snp.sub$snp )
		return(snp.mat);
	}
	
	r <- fgwas_filter( n.snp, n.ind, get_sub_snpmat, pd, pd$phe.mat, Y.name, Z.name, covar.names, op.cpu, fgwas.cutoff, lasso_method );
	return(r);
}

snpmat_fgwas_filter<-function( phe.mat, snp.mat, Y.name, Z.name, covar.names, op.cpu=1, fgwas.cutoff=0.05, lasso_method="BLS")
{
	n.snp <- NROW( snp.mat );
	n.ind <- NCOL( snp.mat ) -2 ;
		
	get_sub_snpmat<- function(snpmat.big, idx.snp)
	{
		return(snpmat.big[idx.snp,,drop=F]);
	}
	
	r <- fgwas_filter( n.snp, n.ind, get_sub_snpmat, snp.mat, phe.mat, Y.name, Z.name, covar.names, op.cpu, fgwas.cutoff, lasso_method );
	return(r);
}

