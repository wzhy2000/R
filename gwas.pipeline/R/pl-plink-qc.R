do.plink.qc1 <- function( plink.obj, options=NULL )
{
	subDir  <- "qc1"
	mainDir <- plink.obj$main.path;

	dir.create( file.path(mainDir, subDir), showWarnings = FALSE );

	options0 <- list( p.mind=0.05, 
			p.geno=0.05, 
			p.maf=0.01, 
			p.hwe=0.000001, 
			p.me=c(0.05, 0.05) );
	
	if(!is.null(options)) options0[names(options)] <- options;
	options <- options0;

	plink.path <- plink.obj$plink.path;
	plink.bfile.nobed <- plink.obj$genotype$plink.bfile.nobed;
	
	plink.out.bfile <- paste("qc1/", plink.bfile.nobed, "-qc1", sep=""); 
	plink.dat.fam   <- paste("qc1/", plink.bfile.nobed, "-qc1.fam", sep=""); 
	plink.out.log   <- paste("qc1/", plink.bfile.nobed, "-qc1-R.log", sep=""); 
	
	plink.out.gen   <- paste( plink.out.bfile, ".gen", sep="");
	plink.out.m10   <- paste( plink.out.bfile, ".m10", sep="");
	
	cat("[ PLINK QC1 ...]\n");
	cat("[ PLINK QC1 ...]\n", file=plink.out.log, append = FALSE);
	
	if( !file.exists( plink.dat.fam  )) 
	{
		t0 <- plink_command(plink.path, c("--noweb", 
										  "--bfile", plink.bfile.nobed, 
										  "--maf",   options$p.maf,
										  "--mind",  options$p.mind,
										  "--geno",  options$p.geno, 
										  "--hwe",   options$p.hwe,
										  "--me",    options$p.me[1], options$p.me[2], 
										  "--make-bed",
										  "--out", plink.out.bfile ) );
		cat(t0, file=plink.out.log, append = TRUE);
	}
	else
		cat(plink.dat.fam, "is found.\n");
		
	file.genome.gz <- paste(plink.out.gen, ".genome.gz", sep="");
	
	if( !file.exists(file.genome.gz  )) 
	{
		t1 <- plink_command(plink.path, c("--noweb", 
										  "--bfile",  plink.out.bfile, 
										  "--Z-genome", 
										  "--out", plink.out.gen) );
		cat(t1, file=plink.out.log, append = TRUE);
	}
	else
		cat(file.genome.gz, "is found.\n");
	
	
	file.out.mds <- paste( plink.out.m10, ".mds", sep="" )
	if( !file.exists( file.out.mds ) ) 
	{
		t2 <- plink_command(plink.path, c("--noweb", 
						"--bfile", plink.out.bfile, 
						"--read-genome",  file.genome.gz, 
						"--mds-plot", "10",
						"--out", plink.out.m10 ));
		cat(t2, file=plink.out.log, append = TRUE);	
	}
	else
		cat(file.out.mds, "is found.\n");
	
	t3 <- draw_pca_plot( file.out.mds, ".qc1", "QC1" );
	cat(t3, file=plink.out.log, append = TRUE);	
	
	plink.obj$qc1 <- list( 
				folder = "qc1",
				p.mind = options$p.mind,
				p.geno = options$p.geno,
				p.maf  = options$p.maf,
				p.hwe  = options$p.hwe,
				p.me   = options$p.me,
				file.out.mds   = file.out.mds,
				file.out.bfile = plink.out.bfile );
	
	return( plink.obj );
}

do.post.qc1 <- function( plink.obj, options=NULL  )
{
	cat("[ PLINK post QC1 ...]\n");

	plink.out.mds <- plink.obj$qc1$file.out.mds
	
	file.qc1.rem.id <- paste(plink.out.mds, ".rem", sep="");
	
	r <- remove_outlier_pca(  plink.out.mds, file.qc1.rem.id, expression( which( abs(tb$C3) > 0.05 ) ) );
	if(!is.null(r))
	{	
		plink.obj$post.qc1 <- list( file.qc1.rem.id = file.qc1.rem.id );
		return( plink.obj );
	}
	else
		return(NULL);
}

do.plink.qc2 <- function( plink.obj, options=NULL )
{
	subDir  <- "qc2/"
	mainDir <- plink.obj$main.path;

	dir.create( file.path(mainDir, subDir), showWarnings = FALSE );

	plink.path 		  <- plink.obj$plink.path;
	plink.bfile.nobed <- plink.obj$qc1$file.out.bfile;
	file.qc1.rem.id   <- plink.obj$post.qc1$file.qc1.rem.id;

	plink.out.bfile <- paste( subDir, plink.obj$genotype$plink.bfile.nobed, "-qc2", sep=""); 
	plink.out.log   <- paste( subDir, plink.obj$genotype$plink.bfile.nobed, "-qc2-R.log", sep=""); 

	plink.dat.fam    <- paste( plink.out.bfile, ".fam",   sep="");
	plink.out.mds    <- paste( plink.out.bfile, ".mds",   sep="");
	plink.out.mds10  <- paste( plink.out.bfile, ".m10", sep="");
	plink.out.log10  <- paste( plink.out.bfile, ".log10", sep="");
	plink.out.line   <- paste( plink.out.bfile, ".line",  sep="");

	cat("[ PLINK QC2 ...]\n");
	cat("[ PLINK QC2 ...]\n", file=plink.out.log, append = FALSE);

	if( !file.exists( plink.dat.fam )) 
	{
		cat("Create", plink.dat.fam, "...\n");
		t1 <- plink_extract_snp( plink.path, plink.out.bfile, plink.bfile.nobed, file.rem.id = file.qc1.rem.id );
		cat( t1, file = plink.out.log, append = TRUE);
	}
	else
		cat(plink.dat.fam, "is found.\n");
		
	file.genome.gz <- paste( plink.out.mds, "genome.gz", sep="." )	
	if( !file.exists( file.genome.gz ) )
	{
		cat("Create",  file.genome.gz , "...\n");
		t2 <- plink_command(plink.path, c("--noweb", 
					"--bfile", plink.out.bfile, 
					"--Z-genome", 
					"--out", plink.out.mds ));
		cat(t2, file = plink.out.log, append = TRUE);	
	}
	else
		cat(file.genome.gz, "is found.\n");
	
	plink.out.mds <- paste( plink.out.mds10, "mds", sep="." );
	if( !file.exists( plink.out.mds ) )
	{
		cat("Create", plink.out.mds, "...\n");
		t3 <- plink_command(plink.path, c("--noweb", 
					"--bfile", plink.out.bfile, 
						"--read-genome",  file.genome.gz, 
						"--mds-plot", "10",
						"--out", plink.out.mds10 ));
		cat(t3, file = plink.out.log, append = TRUE);	
	}
	else
		cat(plink.out.mds, "is found.\n");
	
	t4 <- draw_pca_plot( plink.out.mds, ".qc2", "QC2" );
	cat(t4, file = plink.out.log, append = TRUE);	

	plink.obj$qc2 <- list(
		folder = "qc2",
		plink.out.mds   = plink.out.mds, 
		plink.out.bfile = plink.out.bfile );
		
	return( plink.obj );
}

do.post.qc2 <- function( plink.obj, options=list(mds.cov=5) )
{
	cat("[ PLINK post QC2 ...]\n");

	subDir  <- "qc2/"
	mainDir <- plink.obj$main.path;

	dir.create( file.path(mainDir, subDir), showWarnings = FALSE );

	plink.dat.fam   <- paste( plink.obj$qc2$plink.out.bfile, "fam", sep=".");
	plink.dat.mds   <- plink.obj$qc2$plink.out.mds;
	file.phe.long   <- plink.obj$phenotype$file.phe.long
    file.pca.cov    <- paste(plink.obj$qc2$plink.out.bfile, ".mds.cov.csv", sep="");
	plink.out.log   <- paste( subDir, plink.obj$genotype$plink.bfile.nobed, "-post-qc2-R.log", sep=""); 
 
	t1 <- replenish_family_phe( plink.obj$phenotype$file.phe.long,  plink.dat.fam,  plink.obj$phenotype$cols.pheno );
	cat(t1, file = plink.out.log, append = FALSE);	

	if( !file.exists( file.pca.cov ) )
	{
		cat("Create", file.pca.cov, "...\n");
		t2 <- make_ibs_cov( file.pca.cov, 
					      file.phe.long, plink.dat.mds,  
						  plink.obj$phenotype$cols.cov,
						  plink.obj$phenotype$cols.time, 
						  3+c(1:options$mds.cov) );
		cat(t2, file = plink.out.log, append = TRUE);	
	}

	plink.obj$qc2$file.pca.cov <- file.pca.cov;
	return(plink.obj);
}

draw_pca_plot<-function( file.mds, qc.id, title )
{
	tb<-read.table( file.mds, header=T)
	
	file.c1c2 <- paste(file.mds, qc.id,".C1C2.pdf", sep="");
	pdf( file.c1c2, width=4);
	plot(tb$C1, tb$C2, type="p", xlab="PCA1", ylab="PCA2", main=title);
	dev.off();
	
	file.c1c3 <- paste(file.mds, qc.id, ".C1C3.pdf", sep="");
	pdf( file.c1c3, width=4);
	plot(tb$C1, tb$C3, type="p", xlab="PCA1", ylab="PCA3", main=title);
	dev.off();

	file.c2c3 <- paste(file.mds, qc.id, ".C2C3.pdf", sep="");
	pdf( file.c2c3, width=4);
	plot(tb$C2, tb$C3, type="p", xlab="PCA2", ylab="PCA3", main=title);
	dev.off();
	
	t <- paste("PCA plot: ", file.c1c2, file.c1c3, file.c2c3, sep=""); 
	return(t);
}

plink_extract_snp<-function(plink.path, plink.out.bfile, plink.bfile.nobed, file.keep.id=NULL,  file.rem.id=NULL)
{
	plink.cmd <- c(plink.path, 
		"--noweb",
		"--bfile",  plink.bfile.nobed, 
		"--make-bed", 
		"--out",   plink.out.bfile )
		
	if (!is.null(file.keep.id))
		plink.cmd <-c( plink.cmd, "--keep",  file.keep.id ); 
	if (!is.null(file.rem.id))
		plink.cmd <-c( plink.cmd, "--remove",  file.rem.id ); 

	t1 <- try(system(paste(plink.cmd, collapse=" "), intern = TRUE))
	return(t1);
}

plink_command<-function(plink.path, plink.parms)
{
	t1 <- try(system(paste(c(plink.path, plink.parms), collapse=" "), intern = TRUE))
	show(t1);
	
	return(t1)
}

remove_outlier_pca<-function( file.out.mds, file.mds.rem, exp.filter )
{
	tb<-read.table( file.out.mds, header=T)
	
	tb.rem <- tb[ eval( exp.filter ), c(1:3), drop=F];

	cat("The following will be removed from PLINK data set.\n");
	show(tb.rem);
	
	tb.rem <- tb.rem[,c(1,2),drop=F];
	if (NROW(tb.rem)>0)
	{
		write.table(tb.rem, file=file.mds.rem, quote=F, row.names=F, col.names=F);
		return( file.mds.rem );
	}
	else
		return( NULL );
}

replenish_family_phe<-function( file.phe.long,  plink.dat.fam, phe.long.cols )
{
	tb.phe      <- read.csv(file.phe.long);
	tb.phe.long <- tb.phe[, phe.long.cols];
	tb.mean     <- rowMeans( tb.phe.long, na.rm=T);
	tb_phe      <- data.frame( shareid=tb.phe[,1], mean=tb.mean );
	
	fam     <- read.table( plink.dat.fam );
	fam_new <- sqldf("select fam.V1, fam.V2, fam.V3, fam.V4, fam.V5, tb_phe.mean from fam left join tb_phe on fam.V2==tb_phe.shareid");

	mean.miss <- which(is.na(fam_new[,6]));
	cat("New family file is generated with ",length(mean.miss), " missing data\n");
	
	write.table(fam_new, file=plink.dat.fam, quote=F, sep=" ", row.names=F, col.names=F);
	
	return(plink.dat.fam);
}


make_ibs_cov<-function( file.out.cov, file.phe.long, plink.file.mds, col.phe.cov, col.phe.age, col.mds)
{
	tb.mds <- read.table( plink.file.mds, header=T)
	tb.phe <- read.csv(file.phe.long)  

	# dont use tb.mds[,2] %in% tb.phe[,1] 
	idx <- match( as.character(tb.mds[,2]), as.character(tb.phe[,1]) );

	if ( length( which(is.na(idx)) )>0)
	{
		tb.mds <- tb.mds[ - which(is.na(idx)),];
		cat("! ", length( which(is.na(idx)) ), "Individuals( from MDS file) can not be found in the phenotype file\n");
	}
	
	t1 <- paste( plink.file.mds, "NROW=", NROW(tb.mds), "NCOL=", NCOL(tb.mds), sep="  " );	
	t2 <- paste( file.phe.long, "NROW=", NROW(tb.phe), "NCOL=", NCOL(tb.phe), sep="  " );
	t3 <- paste( "MATCHED=", length(idx), "MISSING=", NROW(tb.phe)-length(idx), sep="  " );
	                                   
	sex <- tb.phe[idx, col.phe.cov, drop=F];
	age <- rowMeans( tb.phe[idx, col.phe.age, drop=F], na.rm=T);
	
	df <- data.frame( shareid=tb.mds[,c(1,2)], sex=sex, age=age, tb.mds[ , col.mds, drop=F]);
	write.csv( df, file=file.out.cov, row.names=F, quote=F);

	t4 <- paste( "File=", file.out.cov, "NROW=", NROW(df), "NCOL=", NCOL(df), sep="  " );

	t <- c( t1, t2, t3, t4 );

	return(t);
}