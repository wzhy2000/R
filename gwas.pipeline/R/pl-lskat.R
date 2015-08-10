do.LSKAT <- function(plink.obj, options=NULL)
{
	library(LSKAT);
	cat("[ LSKAT ...]\n");
	
	subDir  <- "lskat"
	mainDir <- plink.obj$main.path;

	dir.create( file.path(mainDir, subDir), showWarnings = FALSE );

	file.long.csv  <- plink.obj$phenotype$file.phe.long;
	file.phe.long  <- "lskat/temp-phenos-long.csv"
	file.phe.time  <- "lskat/temp-phenos-time.csv"

	tb <- read.csv( file.long.csv, header=T );
	write.csv(tb[,c(1,plink.obj$phenotype$cols.pheno)], file=file.phe.long, quote=F, row.names=F);
	write.csv(tb[,c(1,plink.obj$phenotype$cols.time)],  file=file.phe.time, quote=F, row.names=F);

	file.phe.cov   <- plink.obj$qc2$file.pca.cov;
	file.gene.set  <- plink.obj$gene$file.gene.hg19;
	
	file.plink.bed <- paste( plink.obj$qc2$plink.out.bfile, "bed", sep=".");
	file.plink.bim <- paste( plink.obj$qc2$plink.out.bfile, "bim", sep=".");
	file.plink.fam <- paste( plink.obj$qc2$plink.out.bfile, "fam", sep=".");

	file.ret.rdata <- "lskat/lskat-ret.rdata"

	options0=list(y.cov.count= NA, y.cov.time= 0, g.maxiter = 10, weights.common=c(0.5,0.5), weights.rare=c(1,25), run.cpp=F, debug=T, n.cpu=7);
	if(!is.null( options ) )
		options0[names(options) ] <- options

	ret<-longskat_gene_plink( file.plink.bed, file.plink.bim, file.plink.fam, file.gene.set, file.phe.long, file.phe.cov, file.phe.time, gene.range = NULL,  options=options0);
	
	save( ret, file=file.ret.rdata );
	
	plink.obj$lskat <- list(rdata=file.ret.rdata);

	return(plink.obj);	
}


