do.blasso<-function( plink.obj, options=list() )
{
	library(gwas.lasso)
	
	subDir  <- "blasso"
	mainDir <- plink.obj$main.path;
	dir.create( file.path(mainDir, subDir), showWarnings = FALSE );
	
	y <- c();
	if ( file.exists(plink.obj$phenotype$file.phe.mean ) )
	{
		file.phe.mean  <- plink.obj$phenotype$file.phe.mean
		y <- read.csv(file.phe.mean)
	}
	else
	{
		tb <- read.csv( plink.obj$phenotype$file.phe.long, header = T );
		bmi.mean <- rowMeans( tb[,plink.obj$phenotype$cols.pheno], na.rm=T );
		age.mean <- rowMeans( tb[,plink.obj$phenotype$cols.time], na.rm=T );
	
		y <- cbind(tb[,1] , bmi.mean);
	}


	if ( !is.null(plink.obj$qc2$file.pca.cov ) )
		file.pca.cov  <- plink.obj$qc2$file.pca.cov;

	n.cov   <- options$covariate.count;
	FAM_Cov <- read.csv( file.pca.cov );
	FAM_Cov <- FAM_Cov[, c(2:(n.cov+2)), drop=F];
	colnames(FAM_Cov) <- c("shareid", paste("COV", 1:n.cov, sep="") );

	y0 <- y [ match( FAM_Cov$shareid, y[,1]), ];
	colnames(y0) <- c("shareid", "Y")
	
	newphe <- cbind(y0, FAM_Cov[,-1, drop=F])
	newphe.csv <- "blasso/temp-phe-cov.csv";
	write.csv( newphe, file=newphe.csv , row.names=F, quote=F);

	show(head(newphe));

	file.plink.bed <- paste( plink.obj$qc2$plink.out.bfile, "bed", sep="." );
	file.plink.bim <- paste( plink.obj$qc2$plink.out.bfile, "bim", sep="." );
	file.plink.fam <- paste( plink.obj$qc2$plink.out.bfile, "fam", sep="." );
	file.ret.rdata <- "blasso/blasso-ret.rdata";

	ret <- bls.plink( newphe.csv, 
			file.plink.bed, 
			file.plink.bim, 
			file.plink.fam, 
			"Y", 
			paste("COV", 1:n.cov, sep=""), 
			refit = TRUE, 
			add.used = T, 
			dom.used = T, 
			fgwas.filter = T,       
            options=list(nParallel.cpu = 7) )

	save(ret, file=file.ret.rdata);
	summary(ret);
	plot(ret);
	
	plink.obj$blasso <- list(rdata=file.ret.rdata);
	return( plink.obj );
}
