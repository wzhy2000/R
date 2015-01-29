read_gen_dataset<-function( file.set )
{
	tb <- read.table(file.set, sep=" ", header=F);
	genes <- unique(tb[,1]);
	return(list(len=length(genes), ids=genes, snps=tb));
}

get_gen_group<-function(gen.lib, idx)
{
	gen.name <- gen.lib$ids[idx];
	snp.idx <- which(gen.lib$snps[,1]==gen.name);
	return(list(name=gen.name, snps=gen.lib$snps[snp.idx,2]))
}

get_gen_family<-function(gen.lib, gen.name)
{
	snp.idx <- which(gen.lib$snps[,1]==gen.name);
	if (length(snp.idx)==0)
		return(NULL)
	else
		return(list(name=gen.name, snps=gen.lib$snps[snp.idx,2]));
}

get_snp_info<-function(idx, snp.mat, gen.tb=NA)
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

get_snp_mat<-function(snp.mat, gen.info)
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

	phe.long <- read.csv(file.phe.long);
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
		phe.time <- read.csv(file.phe.time);
		cat("  PHE TIME =", file.phe.time, "\n");
		cat("* Individuals =", NROW(phe.time), "\n");
		cat("* Times =", NCOL(phe.time)-1, "\n");
		cat("* Mean =",  colMeans(phe.time[,-1], na.rm=T), "\n");
		cat("* SD =",    colSds(phe.time[,-1], na.rm=T),"\n");
		idx.na <- which( rowSums(is.na(phe.time[,-1]))==NCOL(phe.time)-1);
		if( length(idx.na)>0) phe.time <- phe.time[ -idx.na, ];
	}
	
	phe.cov <- read.table(file.phe.cov, sep=" ", header=F);
	cat("  PHE COV =", file.phe.cov, "\n");
	cat("* Individuals =", NROW(phe.cov), "\n");
	cat("* Covariate =", NCOL(phe.cov)-2, "\n");
	cat("* Mean =",  colMeans(phe.cov[,-c(1,2)], na.rm=T), "\n");
	cat("* SD =",    colSds(phe.cov[,-c(1,2)], na.rm=T), "\n");
	
	ids.fam<-as.integer(rownames(snp.mat$fam));
	cat("  GENO FAM =", file.plink.fam, "\n");
	cat("* Individuals =", NROW(phe.cov), "\n");

	ids.phe <- intersect(phe.long[,1], phe.cov[,2]);
	if(!is.null(phe.time))
		ids.phe <- intersect(ids.phe, phe.time[,1]);
	
	ids.set <- intersect(ids.phe, ids.fam);
	cat("  COMMON Individuals=", length(ids.set), "\n");
	
	# eg. c(10:1)[match(c(4, 6,8,2,3), c(10:1))]
	
	idx.long <- match( ids.set, phe.long[,1] );
	phe.long <- phe.long[idx.long, ];
	
	idx.cov <- match( ids.set, phe.cov[,2] );
	phe.cov <- phe.cov[idx.cov, ];

	if(!is.null(phe.time))
	{
		idx.time <- match( ids.set, phe.time[,1] );
		phe.time <- phe.time[idx.time, ];
	}

	if (!all(ids.set==ids.fam) )
	{
		idx.fam <- match( ids.set, ids.fam );
		snp.mat$genotypes<- snp.mat$genotypes[idx.fam,]
		snp.mat$fam      <- snp.mat$fam[idx.fam,]
	}
	
	if( !is.null(phe.time) && !all(phe.long[,1]==phe.time[,1]) )
		stop("! ID MATCH ERROR between PHE.LONG and PHE.TIME. \n");
	
	if (!( all(phe.long[,1]==phe.cov[,2]) && all(phe.long[,1]==rownames(snp.mat$fam)) ) )
		stop("! ID MATCH ERROR among 3 files( PHE.LONG, PHE.COV, PLINK.FAM). \n");

	return(list(snp.mat=snp.mat, phe.long=phe.long, phe.time=phe.time, phe.cov = phe.cov));
}

check_plink_file<-function( file.plink.bed, file.plink.bim, file.plink.fam )
{
	cat("Checking PLINK file......\n");
	cat("* BED file =", file.plink.bed, "\n");
	cat("* BIM file =", file.plink.bim, "\n");
	cat("* FAM file =", file.plink.fam, "\n");
	
	library(snpStats);

	snp.mat <- try( read.plink( file.plink.bed,  file.plink.bim, file.plink.fam) );
	if(class(snp.mat)=="try-error")
	{
		return(list(bSuccess=F));
	}
	
	n.idv <- NROW(snp.mat$fam)
	n.snp <- NCOL(snp.mat$genotypes)
	cat("* Individuals =", n.idv, "SNPs=", n.snp, "\n")
	cat("* PLINK loading successfully.\n")
	
	return(list(bSuccess=T, family=snp.mat$fam));
}
	
check_pheno_file<-function( file.phe.long, file.phe.time, chk.family ) 
{
	cat("Checking phenotype file......\n");
	cat("* PHE.LONG =", file.phe.long , "\n");
	cat("* PHE.TIME =", file.phe.time , "\n");

	phe.long <- try( read.csv(file.phe.long) );
	if (class(phe.long)=="try-error")
	{
		cat("! Can not open file(", file.phe.long, ")\n");
		return(list(bSuccess=F));
	}
	
	phe.time <- NULL;
	if(!is.null(file.phe.time))
	{
		phe.time <- try( read.csv(file.phe.time) );
		if (class(phe.time)=="try-error")
		{
			cat("! Can not open file(", file.phe.time, ")\n");
			return(list(bSuccess=F));
		}
	}
	
	if(!is.null(phe.time))
	{
		idx.inter <- intersect( phe.long[,1], phe.time[,1] );
		if( !( length(idx.inter)==length(phe.long[,1]) && length(idx.inter)==length(phe.time[,1] ) ) )
		{
			cat("! PHE.LONG don't have consistent IDs with PHE.TIME.\n" );
			return(list(bSuccess=F));
		}		
	}

	ids <- as.integer(rownames(chk.family));
	phe.idx <- c();
	for(i in 1:length(ids))
	{
		m.idx <- which(phe.long[,1]==ids[i]);
		if (length(m.idx)==1)
			phe.idx<-c(phe.idx, m.idx[1])
		else
		{
			cat("! ID:", ids[i], "can not be found in the phenotype file.\n" );
			next;
		}
		
		if (all(is.na(phe.long[ m.idx[1],-1] )))
			cat("! ID:", ids[i], "all time points are missed.\n" );
	}

	return(list(bSuccess=T))
}

check_covariate_file<-function( file.phe.cov, chk.family, n.cov ) 
{
	cat("Checking covariate file......\n");
	cat("* COV.FILE =", file.phe.cov , "\n");
	
	phe.cov <- try( read.table(file.phe.cov, sep=" ", header=F) );
	if (class(phe.cov)=="try-error")
		return(list(bSuccess=F));

	if( NCOL(phe.cov) < 2 + n.cov)
	{
		cat("! Insufficient covariates in the covariate file, ", NCOL(phe.cov), "<", 2 + n.cov, ".\n" );
		return(list(bSuccess=F));
	}	
	
	ids <- as.integer(rownames(chk.family));
	phe.idx<-c();
	for(i in 1:length(ids))
	{
		m.idx <- which( phe.cov[,2] == ids[i] );
		if (length(m.idx)>=1)
			phe.idx<-c(phe.idx, m.idx[1])
		else
		{
			m.idx <- which( phe.cov[,1] == ids[i] );
			if (length(m.idx)>=1)
				phe.idx<-c(phe.idx, m.idx[1])
			else
			{
				cat("! ID:", ids[i], "can not be found in the covariate file.\n" );
				return(list(bSuccess=F));
			}
		}
	}

	cat("* Individuals =", NROW(phe.cov), "\n");
	cat("* COV.LEN =", NCOL(phe.cov)-2, "\n");
	cat("* COV.MEAN =", colMeans(phe.cov[,-c(1:2), drop=F], na.rm=T), "\n");
	cat("* COV.SD =", colSds(phe.cov[,-c(1:2), drop=F], na.rm=T), "\n");	

	return(list(bSuccess=T))
}

check_geneset_file<-function( file.gene.set ) 
{
	cat("Checking gene definition file......\n");
	cat("* GEN.SET.FILE =", file.gene.set , "\n");

	tb <- try( read.table(file.gene.set, sep=" ", header=F) );
	if(class(tb)=="try-error")
		return(list(bSuccess=F));

	genes <- unique(tb[,1]);
	
	cat("* GENEs =", length(genes), "\n");
	cat("* SNPs =", NROW(tb), "\n");
	
	return(list(bSuccess=T))
}