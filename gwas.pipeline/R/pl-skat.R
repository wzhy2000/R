do.SKAT<-function( plink.obj, options=list(qc.method="qc2") )
{
	cat("[ SKAT ...]\n");

	library(SKAT)
	
	subDir  <- "skat"
	mainDir <- plink.obj$main.path;

	dir.create( file.path(mainDir, subDir), showWarnings = FALSE );

	plink.bfile <- plink.obj$genotype$plink.bfile.nobed;
	if(options$qc.method=="qc2") plink.obj$genotype$qc2
	if(options$qc.method=="impute") plink.obj$genotype$impute
		
	file.plink.bed <- paste( plink.bfile, "bed", sep="." );
	file.plink.bim <- paste( plink.bfile, "bim", sep="." );
	file.plink.fam <- paste( plink.bfile, "fam", sep="." );

	if ( !is.null( plink.obj$gene$file.gene.hg19 ) )
		file.gene.set <- plink.obj$gene$file.gene.hg19
	else
	{
		cat("Gene set file is not set.");
		return(NULL);
	}	
		
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

	n.cov  <- options$covariate.count;
	FAM_Cov <- read.csv( file.pca.cov );
	FAM_Cov <- FAM_Cov[, c(2:(n.cov+2)), drop=F];
	colnames(FAM_Cov) <- c("shareid", paste("COV", 1:n.cov, sep="") );

	y0 <- y [ match( FAM_Cov$shareid, y[,1]), 2 ];

	# show first 6 rows
	show(head( FAM_Cov ));

	obj <- SKAT_Null_Model(y0 ~., data=cbind(y0=y0, FAM_Cov[,-1, drop=F]), out_type="C");

	SKAT.file.SSD       <- "skat/skat-gen-hg19.SSD"
	SKAT.file.Info      <- "skat/skat-gen-hg19.SSD.info"
	file.ret.rdata      <- "skat/skat-ret.rdata"

	Generate_SSD_SetID( file.plink.bed, file.plink.bim, file.plink.fam, file.gene.set, SKAT.file.SSD, SKAT.file.Info)

	FAM <- Read_Plink_FAM( file.plink.fam, Is.binary=FALSE)

	SSD.INFO<-Open_SSD( SKAT.file.SSD, SKAT.file.Info)
	
	cat("sample=", SSD.INFO$nSample, "nSets=", SSD.INFO$nSets);

	summary(obj);

	out2 <- SKAT_CommonRare.SSD.All(SSD.INFO, obj);

	summary(out2);

	save( out2, file=file.ret.rdata );
	
	Close_SSD();
	
	plink.obj$skat <- list(rdata=file.ret.rdata);
	return( plink.obj );
}

draw_skat<-function( plink.obj )
{
	library(LSKAT);

	rdata.lskat <- "lskat-bmi-c1c2.rdata"
	rdata.skat  <- "skat-bmi-c1c2.rdata";
	pdf.lskat   <- "lskat-bmi-c1c2.pdf";
	pdf.skat    <- "skat-bmi-c1c2.pdf";

	load(rdata.lskat);
	summary(ret);
	plot(ret, pdf.file=pdf.lskat);

	library(sqldf);

	rlskat <- ret$gene;
	colnames(rlskat)<-c("id", "name", "chr", "min_pos", "snp", "rare", "Q", "pv");

	load(rdata.skat);
	rskat <- out2$results[,c(1,2)];
	colnames(rskat)<-c("name", "pv");

	rlskat <- sqldf("select rlskat.id, rlskat.name, rlskat.chr, rlskat.min_pos, rlskat.snp, rlskat.rare, rlskat.Q, rskat.pv from rlskat left join rskat on rskat.name=rlskat.name");    

	colnames(rlskat)<-c("id", "name", "chr", "min.pos", "snp", "rare", "Q", "pv");
	ret$gene <- rlskat;

	plot(ret, pdf.file=pdf.skat);
	summary(ret);
}