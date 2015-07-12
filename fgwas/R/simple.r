fg.simple.load<-function( file.simple.geno, file.pheY.csv, file.pheX.csv=NULL, file.pheT.csv=NULL, no.curve=NULL, no.covar=NULL, time.points = NULL, snp.range = NULL )
{
	check_geno_dat( file.simple.geno );
	check_phe_csv( file.pheY.csv);

	if( !missing(file.pheX.csv)) check_phe_csv( file.pheX.csv);
	if( !missing(file.pheT.csv)) check_phe_csv( file.pheT.csv);
	if( !missing(time.points)) check_time_points( time.points );
	if( !missing(snp.range)) check_snp_range( snp.range );
	if( !missing(no.curve)) check_no_curve(no.curve );
	if( !missing(no.covar)) check_no_covar( no.covar );

	pheT.same <- ( is.null(file.pheT.csv) || file.pheT.csv=="" );

	phe.mat <- fg_read_phe( file.pheY.csv, file.pheX.csv, file.pheT.csv, time.points );    
	if ( is.null(phe.mat) )  
	     stop( paste("Can't find phenotype file[", file.pheY.csv, file.pheX.csv, file.pheT.csv, "]") );

	r.dat <- fg_read_simple( file.simple.geno, phe.mat, snp.range );
	
	r.phe <- fg_split_phemat(r.dat$phe.mat, pheT.same, time.points );

	r.est <- fg_dat_est( r.phe$pheY, r.phe$pheX, r.phe$pheT, no.curve, no.covar )
	if ( r.est$error ) 
	     stop( "Can not estimate the parameter of mean vector or covariance." );

	r.dat$phe.mat         <- r.dat$phe.mat;
	r.dat$snp.mat         <- r.dat$snp.mat;
	r.dat$snp.range       <- r.dat$snp.range;
	
	r.dat$phe.est         <- r.est;
	r.dat$times           <- time.points;
	r.dat$file.geno.dat   <- file.simple.geno;
	r.dat$file.pheY.csv   <- file.pheY.csv;
	r.dat$file.pheT.csv   <- file.pheT.csv;
	r.dat$file.pheX.csv   <- file.pheX.csv;
	
	r.dat$pheY	      <- r.phe$pheY;
	r.dat$pheX	      <- r.phe$pheX;
	r.dat$pheT	      <- r.phe$pheT;
	r.dat$pheT.same	      <- pheT.same;
	
	class(r.dat)          <- "FGWAS.SIMPLE.DAT"

	return( r.dat );
}

fg_read_phe <- function(file.pheY.csv, file.pheX.csv, file.pheT.csv, time.points )
{
	phe.mat <- try( read.csv( file.pheY.csv, header=T) );
	if(class(phe.mat)=="try-error")
		return(NULL);
	
	colnames(phe.mat) <- c("ID", paste("Y", 1:(NCOL(phe.mat)-1), sep="_"));
	rownames(phe.mat) <- phe.mat[,1];
	
	if(is.null(time.points)) time.points <- c(1:(NCOL(phe.mat)-1));
	
	if(!is.null(file.pheX.csv) )
	{
		pheX <- try( read.csv( file.pheX.csv, header=T) );
		if(class(pheX)=="try-error")
			return(NULL);
		
		colnames(pheX) <- c("IDX", paste("X", 1:(NCOL(pheX)-1), sep="_"));
		ordY <- match( phe.mat[,1], pheX[,1] );
		phe.mat <- cbind( phe.mat, pheX[ordY, , drop=F] );
		if ( length( which(is.na(phe.mat[,c("IDX")])))>0 )
			phe.mat <- phe.mat[-which(is.na(phe.mat[,c("IDX")])),,drop=F ];
			
		phe.mat$IDX <- NULL;	
	}
	if(!is.null(file.pheT.csv) )
	{
		pheT <- try( read.csv( file.pheT.csv, header=T) );
		if(class(pheX)=="try-error")
			return(NULL);

		colnames(pheT) <- c("IDT", paste("T", 1:(NCOL(pheT)-1), sep="_"));
		ordT <- match( phe.mat[,1], pheT[,1] );
		phe.mat <- cbind( phe.mat, pheT[ordT, , drop=F] );
		if ( length( which(is.na(phe.mat[,c("IDT")])))>0 )
			phe.mat <- phe.mat[-which(is.na(phe.mat[,c("IDT")])),,drop=F ];
		phe.mat$IDT <- NULL;	
	}
	else
	{
		pheT <- t(replicate( NROW(phe.mat), time.points));
		colnames(pheT) <- c( paste("T", 1:(NCOL(pheT) ), sep="_") );
		phe.mat <- cbind( phe.mat, pheT );
	}
	
	return(phe.mat);
}

#scaffoldId, Loci, RefBase, AltBase, P1.A, P1.B, P2.A, P2.B....
fg_read_simple <- function( file.geno.dat, phe.mat, snp.range )
{
	tb.gen <- read.table( file.geno.dat, header = T);

	cat("gen.table[", dim(tb.gen), "]\n");
	
	gen.ids <- colnames(tb.gen)[-c(1:4)];
	phe.ids <- phe.mat[,1];
	temp.idx <- match( gen.ids, phe.ids );
	if(length(which(is.na(temp.idx)))>0)
	{
		cat("* No phenotype data can be found, remove genotype data,", gen.ids[which(is.na(temp.idx))], ".\n");
		tb.gen <- tb.gen[, -which(is.na(temp.idx))+4, drop=F ];
	}
	
	temp.idx <- match( phe.ids, gen.ids );
	if(length(which(is.na(temp.idx)))>0) 
	{
		cat("* No genotype data can be found, remove phenotype data,", phe.ids[which(is.na(temp.idx))], ".\n");
		phe.mat<- phe.mat[-which(is.na(temp.idx)),, drop=F ];
	}

	gen.ids <- colnames(tb.gen)[-c(1:4)];
	phe.ids <- phe.mat[,1];
	temp.idx <- match( gen.ids, phe.ids );
	phe.mat <- phe.mat[ temp.idx, ,drop=F]
	
	if (!is.null(snp.range))
	{
		if ( is.na(snp.range[1]) || snp.range[1] <1 ) snp.range[1] <- 1;
		if ( is.na(snp.range[2]) || snp.range[1] > NROW(tb.gen) ) snp.range[2] <- NROW(tb.gen);
		snp.range <- sort(snp.range)

		cat("* SNP range,", snp.range, ".\n");
		tb.gen <- tb.gen[snprange[1]:snprange[2]]
	}	
	else
		snp.range<-c(1, NROW(tb.gen) );
		
	n.snp    <- NROW(tb.gen)
	n.obs    <- NROW(phe.mat);
	snp.info <- tb.gen[, 1:4];
	
	get_geno_code<-function( d.snp )
	{
		snpB <- as.character(unlist(d.snp[4]));
		snpA <- as.character(unlist(d.snp[3]));
		d.gen2 <- d.snp[-c(1:4)];
		
		QQ2<- paste(snpB, snpB, sep=""); 
		qq0<- paste(snpA, snpA, sep="");
		Qq1<- c(paste(snpA, snpB, sep=""), paste( snpB, snpA, sep="") ) ;
		
		d.g <- rep( NA, length(d.gen2) );
		d.g[which(d.gen2==QQ2)]<-2;
		d.g[which(d.gen2==qq0)]<-0;
		d.g[which(d.gen2==Qq1[1])]<-1;
		d.g[which(d.gen2==Qq1[2])]<-1;
		
		if( mean(d.g, na.rm=T)/2 > 0.5) d.g <- 2 - d.g;
		
		return(d.g);
	}
	
	snp.mat <- apply( tb.gen, 1, get_geno_code );
	# HERE snp.mat is inversed
	snp.mat <- t(snp.mat);
	colnames(snp.mat) <- colnames( tb.gen )[-c(1:4)];
	snp.mat <- cbind( snp.info[, c(1,2), drop=F ], snp.mat ); 	
	
	return( list( n.obs=n.obs, n.snp=n.snp, snp.info=snp.info, snp.mat=snp.mat, phe.mat=phe.mat, snp.range = snp.range) )
}

fg_split_phemat<-function(phe.mat, pheT.same, time.points)
{
	y.sub <- grep( "Y_", colnames(phe.mat) ) ;
	pheY  <- as.matrix( phe.mat[, y.sub, drop = F] );
	rownames(pheY) <- rownames(phe.mat);
	
	pheX  <- NULL;
	if( length( grep( "X_", colnames(phe.mat) ) )>0)
	{
		x.sub <- grep( "X_", colnames(phe.mat) ) ;
		pheX <- as.matrix( phe.mat[, x.sub, drop = F] );
		rownames(pheX) <- rownames(phe.mat);
	}
	
	pheT  <- NULL;
	if(length( grep( "T_", colnames(phe.mat) ) )>0)
	{
		t.sub <- grep( "T_", colnames(phe.mat) ) ;
		pheT <- as.matrix( phe.mat[, t.sub, drop = F] );
		rownames(pheT) <- rownames(phe.mat);
	}
	
	if( pheT.same )
	{
		if(!is.null( time.points ))
			pheT <- time.points
		else
			pheT <- c( 1:NCOL(pheY) );
	}
	
	return(list(pheY=pheY, pheX=pheX,pheT=pheT));
}

summary.FGWAS.SIMPLE.DAT <- function(fg.dat)
{

#todo

}
