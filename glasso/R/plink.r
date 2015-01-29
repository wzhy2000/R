load_plink_binary<-function(file.plink.bed,  file.plink.bim, file.plink.fam, file.phe.long )
{
	snp.mat <- try( read.plink( file.plink.bed,  file.plink.bim, file.plink.fam) );
	if(class(snp.mat)=="try-error")
		return(NULL);
	
	phe.log <- try( read.csv(file.phe.long, header=T) );
	if(class(phe.log )=="try-error")
		return(NULL);
	
	gen.ids <- rownames(snp.mat$fam) ;

	# !!! First Column must be individual ID.
	phe.ids <- c(phe.log[,1]);
	phe.log <- phe.log[, -1, drop=F];
	rownames( phe.log ) <- phe.ids
	
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
		snp.mat$genotypes<- snp.mat$genotypes[-gen.rem.idx,]
		snp.mat$fam      <- snp.mat$fam[-gen.rem.idx ,]
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
		snp.mat$genotypes<- snp.mat$genotypes[-phe.rem.idx,]
		snp.mat$fam      <- snp.mat$fam[-phe.rem.idx ,]
	}
	
	return(list(snp.mat=snp.mat, phe.mat=phe.log));
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

		return( s.mat.i );
	}
	

	snp.imp <- t( apply(s.mat, 2, f.impute) );
	snp.maf <- rowMeans(snp.imp) /2;

	map <- snp.mat$map[snp.set.idx, ,drop=F];
	rownames(snp.imp) <-rownames( map );

	return(list(maf=snp.maf, snp=snp.imp, info=map[,c(2,1,4)]) );
}