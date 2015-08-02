read_gen_dataset<-function( file.set, file.bim )
{
	tb.bim <- read.table(file.bim);
	snp <- tb.bim$V2;
	
	tb.gen <- read.table(file.set, sep=" ", header=F);
	idx <- match(tb.bim$V2, tb.gen$V2)
	
	genes <- unique(tb.gen[idx,1]);
	
	return(list(len=length(genes), ids=genes, snps=tb.gen[idx,]));
}

get_gen_group<-function(gen.list, idx)
{
	gen.name <- gen.list$ids[idx];
	snp.idx <- which(gen.list$snps[,1]==gen.name);
	return(list(name=gen.name, snps=gen.list$snps[snp.idx,2]))
}

get_gen_family<-function(gen.lib, gen.name)
{
	snp.idx <- which(gen.lib$snps[,1]==gen.name);
	if (length(snp.idx)==0)
		return(NULL)
	else
		return(list(name=gen.name, snps=gen.lib$snps[snp.idx,2]));
}

get_snp_plink_info<-function(idx, snp.mat, gen.tb=NA)
{
	s <- snp.mat$genotypes[, idx, drop=F ]

	s.mat <- as(s, "numeric");
	s.mat.i <- s.mat[,1];
	s.miss <- which( is.na(s.mat.i) );
	s0 <- which( s.mat.i == 0 );
	s1 <- which( s.mat.i == 1 );
	s2 <- which( s.mat.i == 2 );

	if (length(s.miss)>0)
		s.mat.i <- s.mat.i[-s.miss];
		
	if ( mean(s.mat.i) > 1 ) s.mat.i <- 2 -s.mat.i;
	snp.imp <- as.matrix( s.mat.i, dim=c(length(s.mat.i),1) ); 
	snp.maf <- sum(snp.imp)/(length(s.mat.i)*2);

	snp.map <- snp.mat$map[idx,]
	gene.name <- "";
	if (!all(is.na(gen.tb)))
	{
		gen.idx <- which( gen.tb[,2]==snp.map$snp.name)
		if (length(gen.idx)>0)
			gene.name <-  gen.tb[gen.idx[1],1];
	}
	
	return(list(snp=snp.imp, maf=snp.maf, name=snp.map$snp.name, chr=snp.map$chromosome, loc=snp.map$position, gene=gene.name, nmiss=length(s.miss), miss=s.miss ) );
}

get_snp_mat<-function(snp.mat, gen.info, snp.impute="mean")
{
	snps <- match(as.character(gen.info$snps), as.character(snp.mat$map[,2]));
	if (length(which(is.na(snps)))>0)
		snps <- snps[-which(is.na(snps))];

	if(length(snps)==0) return(NULL);

	s.mat <- as( snp.mat$genotypes[, snps, drop=F ], "numeric");
	snp.imp <-c();
	snp.maf <- c();
	snp.names <- c();
	
	for(i in 1:dim(s.mat)[2])
	{
		s.mat.i <- s.mat[,i] ;
		s.miss <- which( is.na(s.mat.i) );

		if (length(s.miss)>0)
		{
			if(snp.impute=="mean")
			{
				s.mat.i[s.miss] <- mean(s.mat.i, na.rm=T);
			}
			else
			{
				n.s0 <- length( which( s.mat.i == 0 ) );
				n.s1 <- length( which( s.mat.i == 1 ) );
				n.s2 <- length( which( s.mat.i == 2 ) );
				n.s  <- length(s.mat.i)

				r.miss<- runif( length(s.miss) );
				r.snp <- rep(2, length(s.miss));
				r.snp[r.miss <= n.s0/n.s ]<-0;
				r.snp[r.miss <= (n.s0 + n.s1)/n.s ]<-1;
				s.mat.i[s.miss] <- r.snp;
			}
		}
		
		if (mean(s.mat.i)/2>0.5) s.mat.i <- 2 - s.mat.i;
	   
		snp.imp <- rbind( snp.imp, s.mat.i );
		snp.maf <- c(snp.maf, mean(s.mat.i)/2);
		snp.names <- c(snp.names, gen.info$snps[i]);
	}
	
	rownames(snp.imp) <- snp.names;
	
	map <- snp.mat$map[snps, ,drop=F];

	return(list(maf=snp.maf, snp=snp.imp, info=map[,c(2,1,4)]) );
}

colSds<-function(mat, na.rm=T)
{
	r<-c();
	for(i in 1:dim(mat)[2])
		r <- c(r, sd(mat[,i], na.rm=na.rm));
	return(r);
}

read_gen_phe_cov<-function(file.plink.bed, file.plink.bim, file.plink.fam, file.phe.long, file.phe.time, file.phe.cov)
{
	library(snpStats);

	snp.mat <- read.plink( file.plink.bed,  file.plink.bim, file.plink.fam);

	phe.long <- read.csv(file.phe.long, header=T);
	cat("  PHE LONG =", file.phe.long, "\n");
	cat("* Individuals =", NROW(phe.long), "\n");
	cat("* Times =", NCOL(phe.long)-1, "\n");
	cat("* Mean =",  colMeans(phe.long[,-1], na.rm=T), "\n");
	cat("* SD =",    colSds(phe.long[,-1], na.rm=T),"\n");
	idx.na <- which( rowSums(is.na(phe.long[,-1]))==NCOL(phe.long)-1);
	if( length(idx.na)>0) phe.long <- phe.long[ -idx.na, ];

	phe.time <- NULL;
	if (!is.null(file.phe.time))
	{
		phe.time <- read.csv(file.phe.time, header=T);
		cat("  PHE TIME =", file.phe.time, "\n");
		cat("* Individuals =", NROW(phe.time), "\n");
		cat("* Times =", NCOL(phe.time)-1, "\n");
		cat("* Mean =",  colMeans(phe.time[,-1], na.rm=T), "\n");
		cat("* SD =",    colSds(phe.time[,-1], na.rm=T),"\n");
		idx.na <- which( rowSums(is.na(phe.time[,-1]))==NCOL(phe.time)-1);
		if( length(idx.na)>0) phe.time <- phe.time[ -idx.na, ];
	}
	
	phe.cov <- read.csv(file.phe.cov, header=T);
	cat("  PHE COV =", file.phe.cov, "\n");
	cat("* Individuals =", NROW(phe.cov), "\n");
	cat("* Covariate =", NCOL(phe.cov)-1, "\n");
	cat("* Mean =",  colMeans(phe.cov[,-c(1)], na.rm=T), "\n");
	cat("* SD =",    colSds(phe.cov[,-c(1)], na.rm=T), "\n");
	
	ids.fam<-as.character(snp.mat$fam$member);
	cat("  GENO FAM =", file.plink.fam, "\n");
	cat("* Individuals =", NROW(phe.cov), "\n");

	ids.phe <- intersect(as.character(phe.long[,1]), as.character(phe.cov[,1]) );
	if(!is.null(phe.time))
		ids.phe <- intersect(ids.phe, as.character(phe.time[,1]) );
	
	ids.set <- intersect(ids.phe, ids.fam);
	cat("  COMMON Individuals=", length(ids.set), "\n");
	
	# eg. c(10:1)[match(c(4, 6,8,2,3), c(10:1))]
	
	idx.long <- match( ids.set, as.character(phe.long[,1]) );
	phe.long <- phe.long[idx.long, ];
	
	idx.cov <- match( ids.set, as.character(phe.cov[,1]) );
	phe.cov <- phe.cov[idx.cov, ];

	if(!is.null(phe.time))
	{
		idx.time <- match( ids.set, as.character(phe.time[,1]) );
		phe.time <- phe.time[idx.time, ];
	}

	if (!all(ids.set==ids.fam) )
	{
		idx.fam <- match( ids.set, ids.fam );
		snp.mat$genotypes<- snp.mat$genotypes[idx.fam,]
		snp.mat$fam      <- snp.mat$fam[idx.fam,]
		ids.fam <- as.character(snp.mat$fam$member);
	}

	if( !is.null(phe.time) && !all( phe.long[,1] == phe.time[,1]) )
		stop("! ID MATCH ERROR between PHE.LONG and PHE.TIME. \n");

	if (!( all(phe.long[,1]==phe.cov[,1]) && all(phe.long[,1]==ids.fam) ) )
		stop("! ID MATCH ERROR among 3 files( PHE.LONG, PHE.COV, PLINK.FAM). \n");

	return(list(snp.mat=snp.mat, phe.long=phe.long, phe.time=phe.time, phe.cov = phe.cov));
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

	r <- write.plink( file.base=snp.file.base, snp.major = F, snps=t(snps), 
	    	id=sub.name, 
	    	father=rep(0,dim.snps[2]), 
	    	mother=rep(0,dim.snps[2]), 
	    	sex=rep(0,dim.snps[2]), 
	    	phenotype=rep(-9,dim.snps[2]), 
			chromosome=chromosome, 
			genetic.distance=position, 
			position= position, 
			allele.1 = rep("A",dim.snps[1]), 
			allele.2 = rep("B",dim.snps[1]), 
			na.code=0);

	cat("Genotype files have been converted into PLINK binary format(bed/bim/fam)\n");

	return(list(file.plink.bed = paste(snp.file.base, ".bed", sep=""),
   	    	file.plink.bim = paste(snp.file.base, ".bim", sep=""),
   	    	file.plink.fam = paste(snp.file.base, ".fam", sep="")));
}

shrink_snpmat<-function(snp.mat, gen.list, gene.range )
{
	snp.mat0 <- snp.mat; 
	
	#snp.idx  <- unlist( lapply( gene.range, function(i) { which( gen.list$snps[,1] == gen.list$ids[i] ) } ) );
	snp.idx <- which(!is.na(match(gen.list$snp[,1], gen.list$id[gene.range])))
	snp.name <- unique( gen.list$snps[snp.idx,2] );
	
	snp.idx0 <- match( as.character(snp.name), as.character(snp.mat$map[,2]));
	if (length(which(is.na(snp.idx0)))>0)
		snp.idx0 <- snp.idx0[-which(is.na(snp.idx0))];

	if(length(snp.idx0)==0) return(NULL);

	snp.mat0$genotypes <- snp.mat$genotypes[, snp.idx0, drop=F ];
	snp.mat0$map <- snp.mat$map[snp.idx0,];
	
	return( snp.mat0 );
}
