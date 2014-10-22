read_gen_dataset<-function( file.set )
{
	tb <- read.table(file.set, sep=" ", header=F);
	genes <- unique(tb[,1]);
	return(list(len=length(genes), ids=genes, snps=tb));
}

get_gen_group<-function(gen, idx)
{
	gen.names <- gen$ids[idx];
	snp.idx <- which(gen$snps[,1]==gen.names);
	return(list(name=gen.names, snps=gen$snps[snp.idx,2]))
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
	if(length(snps)==0) return(NA);

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

read_gen_phe_cov<-function(file.plink.bed, file.plink.bim, file.plink.fam, file.phe.long, file.phe.cov)
{
	library(snpStats);

	snp.mat <- read.plink( file.plink.bed,  file.plink.bim, file.plink.fam);

	phe.log <- read.csv(file.phe.long);
	ids<-as.integer(rownames(snp.mat$fam));
	phe.idx <- c();
	rem.idx <- c();
	for(i in 1:length(ids))
	{
		m.idx <- which(phe.log[,1]==ids[i]);
		if (length(m.idx)==1)
			phe.idx<-c(phe.idx, m.idx[1])
		else
			stop("*ID:", ids[i], "dont have longitudinal phenos.\n" );
				
		if (all(is.na(phe.log[m.idx[1],-1])))
			rem.idx <-c(rem.idx, i);
	}
	phe.log <- phe.log[phe.idx, ];

	phe.cov <- read.table(file.phe.cov, sep=" ", header=F);
	phe.idx<-c();
	for(i in 1:length(ids))
	{
		m.idx <- which(phe.cov[,2]==ids[i]);
		if (length(m.idx)>=1)
			phe.idx<-c(phe.idx, m.idx[1])
		else
		{
			m.idx <- which( phe.cov[,1] == ids[i] );
			if (length(m.idx)>=1)
				phe.idx<-c(phe.idx, m.idx[1])
			else
				cat("! ID:", ids[i], "can not be found in the covariate file.\n" );
		}
	}
	phe.cov <- phe.cov[phe.idx, ]
	
	if(length(rem.idx)>0)
	{
		cat("!", length(rem.idx), "individuals are removed from pheotype and genotype data.\n")
		phe.log <- phe.log[-rem.idx, ]
		phe.cov <- phe.cov[-rem.idx, ]
		snp.mat$genotypes<- snp.mat$genotypes[-rem.idx,]
		snp.mat$fam      <- snp.mat$fam[-rem.idx,]
	}

	return(list(snp.mat=snp.mat, phe.log=phe.log, phe.cov = phe.cov));
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
	
check_pheno_file<-function( file.phe.long, chk.family ) 
{
	cat("Checking phenotype file......\n");
	cat("* PHE.FILE =", file.phe.long , "\n");

	phe.log <- try( read.csv(file.phe.long) );
	if (class(phe.log)=="try-error")
		return(list(bSuccess=F));
	
	ids <- as.integer(rownames(chk.family));
	phe.idx <- c();
	rem.idx <- c();
	for(i in 1:length(ids))
	{
		m.idx <- which(phe.log[,1]==ids[i]);
		if (length(m.idx)==1)
			phe.idx<-c(phe.idx, m.idx[1])
		else
		{
			cat("! ID:", ids[i], "can not be found in the phenotype file.\n" );
			return(list(bSuccess=F));
		}
		
		if (all(is.na(phe.log[m.idx[1],-1])))
			cat("! ID:", ids[i], "all time points are missed.\n" );
	}

	phe.log <- phe.log[phe.idx, ];

	cat("* Individuals =", NROW(phe.log), "\n");
	cat("* PHE.TIMES =", NCOL(phe.log)-1, "\n");
	cat("* PHE.MEAN =", colMeans(phe.log[,-1], na.rm=T), "\n");
	cat("* PHE.SD =", colSds(phe.log[,-1], na.rm=T),"\n");

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