
library(LSKAT);
library(SKAT);
library(mvtnorm);

g.snp.hap1  <- "skat-test-1.hap";
g.snp.pos2  <- "skat-test-1.pos"

get_con_param<-function(parm.id)
{
        for (e in commandArgs())
        {
                ta = strsplit(e,"=", fixed=TRUE);
                if(! is.na( ta[[1]][2]))
                {
                        temp = ta[[1]][2];
                        if( ta[[1]][1] == parm.id) {
                                temp = as.character(temp);
                                return (temp);
                        }
                }
        }

        return(NA);
}

check.type1err<-function( n.loop, n.rep, n.sample, snprange, w.rare, w.common, par, f.simu)
{
	n.minsect<- snprange[1];
	n.maxsect<- snprange[2];
	snp.mat  <- get_snp_mat(g.snp.hap1, g.snp.pos2, n.minsect, n.maxsect, n.sample, n.rep)

	rs <- c();
	for(i in 1:n.loop)
	{
		phe <- simu.long.phe.none( n.sample, par, f.simu);

		h0 <- longskat_est_model( phe$y, phe$cov, bTime=0, g.maxiter=10 );
		
		par_null <- c(h0$sig_a, h0$sig_b, h0$sig_e, h0$rho)

		for(j in 1:n.rep)
		{

			ls <- longskat_gene_run( h0$y.delt, phe$cov, t(snp.mat[[j]]), par_null, weights.rare=w.rare, weights.common=w.common);

cat("TYPE1.i=", i, j, ls$pv, n.sample, g.phe.len, ls$snp.total, ls$snp.rare,  h0$sig_a, h0$sig_b, h0$sig_e, h0$rho,  h0$cov, h0$val, "\n\n");

			rs <- rbind(rs, c(i, j, ls$pv, n.sample, g.phe.len, ls$snp.total, ls$snp.rare, h0$sig_a, h0$sig_b, h0$sig_e, h0$rho,  h0$cov, h0$val, ls$qv ) );
		}

	}
	
	colnames(rs)<-c("i","j", "pv","sample","phe","snp","rare","sig_a", "sig_b","sig_e", "rho", "a", "b", "LR", "Q")
	
	save(rs, file=paste("simu-type1-", "0.rdata", sep="") );
	
	return(rs);
}

check.power<-function( n.loop, n.sample, snprange, w.rare, w.common, par, f.simu)
{
	n.minsect<- snprange[1];
	n.maxsect<- snprange[2];
	snp.mat <- get_snp_mat(g.snp.hap1, g.snp.pos2, n.minsect, n.maxsect, n.sample, n.loop)

	rs <- c();
	for(i in 1:n.loop)
	{
		phe <- simu.long.phe.power( n.sample, snp.mat[[i]], par, f.simu);

		# LONG-SKAT( Gene ) 
		h0 <- longskat_est_model( phe$y, phe$cov, bTime=0, g.maxiter=2 );
		par_null <- c(h0$sig_a, h0$sig_b, h0$sig_e, h0$rho);
		
		ls <- longskat_gene_run( h0$y.delt, phe$cov, t(snp.mat[[i]]), par_null, weights.rare=w.rare, weights.common=w.common );

		# LONG-SKAT( SNP ) 
		#ls0 <- check.power.lskat.snp( phe, snp.mat[[i]], h0, par_null, w.rare, w.common );
		ls0 <- ls;

		# SKAT( Baseline ) 
		phe.bl <- data.frame(phe$y[,1], phe$cov) ;
		colnames(phe.bl)<-c("Y","X1","X2");
		obj.bl <- SKAT_Null_Model(Y ~ 1 + X1 + X2, data=phe.bl, out_type="C")
		r.bl   <- SKAT( as.matrix(snp.mat[[i]]), obj.bl );

		# SKAT( MEAN ) 
		mu <- apply( phe$y, 1, function(x){mean(x)} );
		phe.mu <- data.frame(mu, phe$cov) ;
		colnames(phe.mu)<-c("Y","X1","X2");
		obj.mu<-SKAT_Null_Model(Y ~ 1 + X1 + X2, data=phe.mu, out_type="C")
		r.mu<-SKAT(as.matrix(snp.mat[[i]]), obj.mu );

cat("POWER.i=", i,  ls$pv, ls0$pv, r.bl$p.value, r.mu$p.value, "par=", h0$sig_a, h0$sig_b, h0$sig_e, h0$rho,  h0$u, h0$cov, "\n\n");

		rs <- rbind(rs, c(i, ls$pv, ls0$pv, r.bl$p.value, r.mu$p.value, n.sample, g.phe.len, ls$snp.total, ls$snp.rare,  
			h0$sig_a, h0$sig_b, h0$sig_e, h0$rho,  h0$u, h0$cov, h0$val, ls$qv ) );
	}
	
	colnames(rs)<-c("i", "lskat-gene", "lskat-snp", "skat-bl", "skat-mu", "sample","phe","snp","rare","sig_a", "sig_b","sig_e", "rho", "u", "a", "b", "LR", "Q")
	
	return(rs);
}

check.power.lskat.snp<-function( phe, snp.mat, h0, par_null, w.rare, w.common )
{
	snp.res <- c();
	for(j in 1:dim(snp.mat)[2])
	{
		snp.vec<- snp.mat[,j,drop=F];
		s.miss <- which(snp.vec!=0 & snp.vec!=1 & snp.vec!=2)
		snp1   <- list( snp=snp.vec, maf=mean(snp.vec[,1])/2, name=j, chr=1, loc=j, gene=1, nmiss=length(s.miss), miss=s.miss );

		ls0 <- longskat_snp_run( h0$y.delt, phe$y, phe$cov, snp1, par_null, 2, weights.rare=w.rare, weights.common=w.common )
		snp.res <- rbind( snp.res, c(j, snp1$maf, snp1$nmiss, ls0$snp.rare, ls0$qv, ls0$pv));
	}

	sig <- which.min(snp.res[,6]);		
	
	return(list(pv=snp.res[sig,6], qv=snp.res[sig,5]));
}

get_snp_mat<-function(snp.file1, snp.file2, min.sect, max.sect, n.sample, mat.count)
{
	snp.hap <- read.table(snp.file1, header=F);
	snp.pos <- read.table(snp.file2, header=T);
	snp.maxpos <- max(snp.pos$CHROM_POS);
	snp.minpos <- min(snp.pos$CHROM_POS);

	rare.cutoff <- 1/sqrt(2*n.sample);
	
	mat.ret <- list();	
	n.mat <- 1;
	while(n.mat<=mat.count)
	{
		snp.start <- as.integer(runif(1, min=snp.minpos, max=snp.maxpos));
		snp.sect  <- as.integer(runif(1, min=min.sect,   max=max.sect));
		
		p.sel<-which( snp.pos$CHROM_POS>snp.start & snp.pos$CHROM_POS<(snp.start+snp.sect) );

		snp.pair <- sample(1:(dim(snp.hap)[2]-2));
		snp.mat1 <- snp.hap[ snp.pair[1:n.sample], p.sel + 2 ]; 
		snp.mat2 <- snp.hap[ snp.pair[(n.sample+1):(2*n.sample)], p.sel + 2 ]; 
		snp.mat <- snp.mat1 + snp.mat2 -2;
	
		maf <- colMeans(snp.mat)/2;
		m.same <- which( maf==1 | maf==0 );
		if (length(m.same)>0)
			snp.mat <- snp.mat[, -m.same, drop=F ];
		
		#make sure of 2: Minor Allels, 0:Major Allels
		maf <- colMeans(snp.mat)/2;
		m.minor <- which( maf>0.5 );
		if (length(m.minor)>0)
			for(i in 1:length(m.minor))
				snp.mat[,m.minor[i]] <- 2 - snp.mat[,m.minor[i]];
			
		maf <- colMeans(snp.mat)/2;
		m.smallmaf <- which( maf < 5/n.sample);
		if (length(m.smallmaf)>0)
			snp.mat <- snp.mat[, -m.smallmaf,drop=F ];

		if (dim(snp.mat)[2]==1)
			next;

		maf <- colMeans(snp.mat)/2;
		n.rare <- length(which(maf<rare.cutoff));
		if (n.rare<1)
			next;
		
		if ( n.rare>=1 && n.rare>length(maf)*1/3 )
		{
			n.rare0 <- as.integer(length(maf)*1/3)+1;
			if (n.rare0>1)
			{
				rare.rm <- which(maf<rare.cutoff)[c(n.rare0:n.rare)];
				if (length(rare.rm)>0) snp.mat <- snp.mat[, -rare.rm, drop=F ];
			}
		}
		
		mat.ret[[n.mat]] <- snp.mat;
		n.mat <- n.mat + 1;
	}
	
	return(mat.ret);
}

simu.long.phe.none<-function( n.sample, par, f.simu )
{
	par_cov <- cbind(1, rnorm(n.sample, 0, 1 ), ifelse(runif(n.sample)>0.5, 0, 1) );

	y <-  f.simu(n.sample, par) + par_cov%*%c( 1, par$a, par$b )%*%t(rep(1, par$times));

	return(list(y=y, cov=par_cov[,-1]))
}

simu.long.phe.power<-function( n.sample, snp.mat, par, f.simu )
{
	par_cov <- cbind(1, rnorm(n.sample, 0, 1 ), ifelse(runif(n.sample)>0.5, 0, 1) );

	y <-  f.simu(n.sample, par) + par_cov%*%c(1, par$a, par$b )%*%t(rep(1, par$times));

	rare.cutoff <- 1/sqrt(2*n.sample);
	maf <- colMeans(snp.mat)/2;
	n.snp <- length(maf);
	snp.rare.total <- length(which(maf<rare.cutoff))
	
	beta.ratio <- cumsum( par$beta.effect );
	snp.comm <- 0;
	snp.rare <- 0;
	for( i in 1:n.snp)
	{
		if (maf[i]<rare.cutoff)
		{
			beta <- par$c1[ which.max(snp.rare/snp.rare.total <= beta.ratio) ]
			if (is.na(beta)) beta<-runif(1, 0.15, 0.2461);
			y <- y + snp.mat[,i] * beta;
			snp.rare <- snp.rare + 1;
		}
		else
		{
			if (snp.comm < 2 )
			{
				y <- y + snp.mat[,i] * 0.27;
				snp.comm <- snp.comm+1;
			}
		}
	}
	return(list(y=y, cov=par_cov[,-1]))
}

#par1 : rho
#par2 : sigma
f.mn.ar1<-function( sample, par)
{
	ncol <- par$times;
	AR.1 <- array(0,dim=c(ncol,ncol));
	for(i in 1:ncol)
	for(j in 1:ncol)
		AR.1[i,j] <- par$par1^abs(i-j);
		
	sigma.mn1 <- par$sig_a^2 + par$par2^2*AR.1;
	sigma.mn2 <- diag(par$sig_e^2, ncol) 
	r <- rmvnorm( sample,  rep(0, ncol), sigma.mn1 + sigma.mn2 )  ;
	
	sigma <- sigma.mn1 + sigma.mn2
cat("MN.AR1---", min(r), max(r), mean(r), min(sigma), max(sigma),mean(sigma), "\n");	
	
	return(r);
}

#par1 : phi
#par2 : v
f.mn.sad<-function( sample, par )
{
	phi <- par$par1;
	v   <- par$par2;
	
	ncol <- par$times;
	sad.1 <- array(0,dim=c(ncol,ncol));
	for(i in 1:ncol)
	for(j in i:ncol)
	{
		sad.1[i,j] <- (1 - phi^(2*(j-i)))/(1-phi^2)*phi^abs(j-i);
		sad.1[j,i] <- sad.1[i,j];
	}
	
	sigma.mn1 <- par$sig_a^2 + sad.1*v^2;
	sigma.mn2 <- diag(par$sig_e^2, ncol); 
	r <- rmvnorm( sample,  rep(0, ncol), sigma.mn1 + sigma.mn2 );

	sigma <- sigma.mn1 + sigma.mn2
cat("MN.SAD---",min(r), max(r), mean(r), min(sigma), max(sigma),mean(sigma), "\n");	

	return(r);
}

#par1 : rho
#par2 : sigma
f.mn.cm<-function( sample, par )
{
	ncol <- par$times;
	AR.1 <- array(par$par1,dim=c(ncol,ncol));
	for(i in 1:ncol)
		AR.1[i,i] <- 1;
		
	sigma.mn1 <- par$sig_a^2 + par$par2^2*AR.1;
	sigma.mn2 <- diag(par$sig_e^2, ncol);
	r <- rmvnorm( sample,  rep(0, ncol), sigma.mn1 + sigma.mn2 );

	sigma <- sigma.mn1 + sigma.mn2
cat("MN.CM---", min(r), max(r), mean(r), min(sigma), max(sigma),mean(sigma), "\n");	

	return(r);
}

f.mt.ar1<-function( sample, par )
{
	ncol <- par$times;
	AR.1 <- array(0,dim=c(ncol,ncol));
	for(i in 1:ncol)
	for(j in 1:ncol)
		AR.1[i,j] <- par$par1^abs(i-j);
		
	sigma.mn <- par$sig_a^2 + par$par2^2*AR.1;
	sigma.t  <- diag(par$sig_e^2, ncol); 
	r <- rmvnorm( sample,  rep(0, ncol), sigma.mn ) + rmvt( sample,  delt=rep(0, ncol), sigma=sigma.t, df=10 );
	
	sigma <- sigma.mn  + sigma.t
cat("MT.AR1---", min(r), max(r), mean(r), min(sigma), max(sigma),mean(sigma), "\n");	

	return(r);

}

f.sn.ar1<-function( sample, par )
{
	library(sn);
	
	ncol <- par$times;
	AR.1 <- array(0,dim=c(ncol,ncol));
	for(i in 1:ncol)
	for(j in 1:ncol)
		AR.1[i,j] <- par$par1^abs(i-j);
		
	sigma.mn <- par$sig_a^2 + par$par2^2*AR.1;
	sigma.sn <- diag(par$sig_e^2, ncol);
	r <- rmvnorm( sample,  rep(0, ncol), sigma.mn ) + rmsn( sample,  xi=rep(0, ncol), sigma.sn, alpha = runif(ncol)*3 );

	sigma <- sigma.mn  + sigma.sn
cat("SN.AR1---", min(r), max(r), mean(r), min(sigma), max(sigma),mean(sigma), "\n");	
	return(r);
}

#par1 : rho1
#par2 : rho2
#par3 : sigma1
#par4 : sigma2
#par5 : ratio
f.mmn.ar1<-function( sample, par )
{
	ncol <- par$times;

	AR1.1 <- array(0,dim=c(ncol,ncol));
	for(i in 1:ncol)
	for(j in 1:ncol)
		AR1.1[i,j] <- par$par1^abs(i-j);
		
	AR1.2 <- array(0,dim=c(ncol,ncol));
	for(i in 1:ncol)
	for(j in 1:ncol)
		AR1.2[i,j] <- par$par2^abs(i-j);

	sigma.mn  <- par$sig_a^2 + AR1.1*par$par3^2*par$par5 + AR1.2*par$par4^2*(1-par$par5);
	sigma.mmn <- diag(par$sig_e^2, ncol); 
	r <-  rmvnorm( sample,  rep(0, ncol), sigma.mn ) + rmvnorm( sample,  rep(0, ncol), sigma.mmn );
	
	sigma <- sigma.mmn  + sigma.mn
cat("MMN.AR1---", min(r), max(r), min(sigma), max(sigma),"\n");		
	
	return(r);
}



