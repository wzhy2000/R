load_plink_binary<-function(file.plink.bed,  file.plink.bim, file.plink.fam, file.phe.long )
{
	plink <- try( read.plink( file.plink.bed,  file.plink.bim, file.plink.fam) );
	if(class(plink)=="try-error")
		return(NULL);
	
	phe.log <- try( read.csv(file.phe.long, header=T) );
	if(class(phe.log )=="try-error")
		return(NULL);
	
	gen.ids <- rownames(plink$fam) ;

	# !!! First Column must be individual ID.
	phe.ids <- phe.log[,1];
	phe.log <- phe.log[, -1, drop=F];
	rownames( phe.log ) <- phe.ids;
	
	phe.idx <- c();
	gen.rem.idx <- c();
	
	for(i in 1:length(gen.ids))
	{
		m.idx <- which(gen.ids[i] == phe.ids );
		if (length(m.idx)==1)
			phe.idx <- c(phe.idx, m.idx[1])
		else
		{
			gen.rem.idx <- c(gen.rem.idx, i);
			cat("*ID:", gen.ids[i], " dont have longitudinal phenos.\n" );
		}
	}

	phe.log <- phe.log[ phe.idx, ,drop=F];
	if(length(gen.rem.idx)>0)
	{
		cat("!", length(gen.rem.idx), "individuals are removed from genotype data due to no phenotypic data.\n")
		plink$genotypes<- plink$genotypes[-gen.rem.idx,]
		plink$fam      <- plink$fam[-gen.rem.idx ,]
	}
	
	phe.rem.idx <- c();
	for(i in 1:NROW(phe.log))
	{
		if (all(is.na(phe.log[m.idx[1],])))
			phe.rem.idx <-c(phe.rem.idx, i);
	}				

	if(length(phe.rem.idx)>0)
	{
		cat("!", length(phe.rem.idx), "individuals are removed from pheotype and genotype data due to invalid measurement.\n")
		phe.log <- phe.log[ -phe.rem.idx, ,drop=F]
		plink$genotypes<- plink$genotypes[-phe.rem.idx,]
		plink$fam      <- plink$fam[-phe.rem.idx ,]
	}
	
	return(list(snp.mat=plink, phe.mat=phe.log));
}

get_sub_snp<-function(snp.mat, snp.set.idx)
{
	s.mat <- as( snp.mat$genotypes[, snp.set.idx, drop=F ], "numeric");
	snp.imp <-c();
	snp.maf <- c();
	snp.names <- c();

	f.impute<-function(s.mat.i )
	{
		s.miss <- which( is.na(s.mat.i) );

		if (length(s.miss)>0)
		{
			n.AA <- length( which( s.mat.i == 0 ) );
			n.Aa <- length( which( s.mat.i == 1 ) );
			n.aa <- length( which( s.mat.i == 2 ) );
			n.s  <- n.AA + n.Aa + n.aa;
			
			r.miss <- runif( length(s.miss) );
			r.snp  <- rep(2, length(s.miss));
			r.snp [ r.miss <= n.AA/n.s ]<-0;
			r.snp [ r.miss <= (n.AA + n.Aa)/n.s ]<-1;
			s.mat.i[s.miss] <- r.snp;
		}
		
		if (mean(s.mat.i)/2>0.5) s.mat.i <- 2 - s.mat.i;

		return( s.mat.i );
	}
	

	snp.imp <- t( apply(s.mat, 2, f.impute) );
	snp.maf <- rowMeans(snp.imp) /2;

	map <- snp.mat$map[snp.set.idx, ,drop=F];
	rownames(snp.imp) <-rownames( map );

	return(list(maf=snp.maf, snp=snp.imp, info=map[,c(2,1,4)]) );
}

impute_simple_snp <- function( snpmat )
{
	f_imputed <- function( snp )
	{
		s.miss <- which( is.na( snp[3:length(snp)] ) );
		if ( length(s.miss)>0 )
		{
			n.AA <- length(which( snp[3:length(snp)] == 0 ) );
			n.Aa <- length(which( snp[3:length(snp)] == 1 ) );
			n.aa <- length(which( snp[3:length(snp)] == 2 ) );
			n.s  <- n.AA + n.Aa + n.aa;
				
			r.miss <- runif( length(s.miss) );
			r.snp  <- rep(2, length( s.miss ));
			r.snp [ r.miss <= n.AA/n.s ]<-0;
			r.snp [ r.miss <= (n.AA + n.Aa)/n.s ]<-1;
			snp [ s.miss + 2 ] <- r.snp;
		}
		
		if (mean(  snp[3:length(snp)] )/2>0.5)  
			snp[ 3:length(snp) ] <- 2 -  snp[ 3:length(snp) ];

		return(snp)
	}

	total_miss <- length(which(is.na(snpmat)));
	
	# !!! Applying matrix[m,n] will be [n.m]!!!
	snpmat <- t ( apply(snpmat, 1, f_imputed) ) ;
	
	cat("* Missing SNPs are imputed(", total_miss, "SNPs).\n");
		
	return(snpmat); 
}


convert_simpe_to_plink <- function( snp.mat, snp.file.base )
{	
	chromosome <- snp.mat[,1];
	position <- snp.mat[,2];

	# PLINK raw data: 1/2/3==> AA,AB,BB, 0==>NA
	snp.mat <- snp.mat[,-c(1,2),drop=F] + 1;

	sub.name <- colnames(snp.mat);
	snp.name <- rownames(snp.mat);
	
	###snps
	dim.snps <- dim(snp.mat);

	snps <- as.raw( as.matrix(snp.mat ) );
	snps <- array(snps, dim=dim.snps);
	colnames(snps) <- sub.name;
	rownames(snps) <- snp.name;
	class(snps) <- "SnpMatrix";
	
	r <- write.plink( file.base=snp.file.base, snp.major = F, snps=t(snps), id=sub.name, chromosome=chromosome, position= position, na.code=0);

	cat("Genotype files have been converted into PLINK binary format(bed/bim/fam)\n");

	return;
}
