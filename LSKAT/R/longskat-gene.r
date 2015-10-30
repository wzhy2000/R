library(mvtnorm);

est_gen_Q.R<-function(Y.delt, Y.time, maf, Z, X, par_null, y.cov.time )
{
	sig_a2  <- par_null[1]^2;
	sig_b2  <- par_null[2]^2;
	sig_e2  <- par_null[3]^2;
	par_rho <- par_null[4];

	n <- dim(Y.delt)[1];
	m <- dim(Y.delt)[2];
	k <- NCOL(Z);
	
	# SKAT   qv,pv
	Q.v<-0;
	# burden qu,pu
	Q.u<-0;
	Q.w<-c();

	AR.1 <- matrix(rep(1:m,m),ncol=m, byrow=T)
	AR.1 <- par_rho^abs(t(AR.1)-AR.1)

	V.j <- array(1, dim=c(m,m)) * sig_a2 +  AR.1 * sig_b2 + diag(1, m) * sig_e2
	V.j_1 <- solve(V.j);

	if(0)
	{
		V_1 <- kronecker( diag(1,n), V.j_1 );

		K1 <-  kronecker( Z, array(1, dim=c(m,1)) ); 
		KY <-  array( c(t(Y.delt)), dim=c(m*n,1));
		Q.v<- t( KY ) %*% V_1 %*% K1 %*% t(K1) %*%V_1 %*% KY /2;

		kr.Z<- kronecker( Z, array(1, dim=c(m,1)) );
		kr.X<- kronecker( X, array(1, dim=c(m,1)) );
		Q.w <- t(kr.Z)%*% V_1 %*% kr.Z - (t(kr.Z)%*% V_1 %*%kr.X) %*% solve(t(kr.X)%*% V_1 %*% kr.X) %*% (t(kr.X)%*% V_1 %*%kr.Z);
	}
	else
	{
		V.j_x <- list();
		Y.j_x <- list();
		M.j_x <- c();

		for(i in 1:n){
			y.i <- Y.delt[i,];
			y.j_i <- V.j[!is.na(y.i), !is.na(y.i)];
			V.j_x[[i]] <- solve(y.j_i);
			Y.j_x[[i]] <- y.i[!is.na(y.i)];
			M.j_x[i] <- length(Y.j_x[[i]]);
		}

		for(i in 1:k){
			Q.i <- 0;
			for(j in 1:n)
				Q.i <- Q.i + Z[j,i]*t(rep(1, M.j_x[j] )) %*% V.j_x[[j]] %*% (Y.j_x[[j]])
			Q.v <- Q.v + (Q.i)^2;
			Q.u <- Q.u + Q.i;
		}
		
		Q.v <- Q.v/2;
		Q.u <- Q.u^2/2;

		n.x <- NCOL(X) + y.cov.time ;
		W0 <- array(0, dim=c(k,k));
		W1 <- array(0, dim=c(k,n.x));
		W2 <- array(0, dim=c(n.x,n.x));
		W3 <- array(0, dim=c(n.x,k));
		
		y.time.max <- max( Y.time, na.rm=T);	
		y.time.min <- min( Y.time, na.rm=T);	
		Y.time <- (Y.time - y.time.min)/(y.time.max - y.time.min)*2 - 1;

		for (i in 1:n){
			m.i <- M.j_x[i] ;

			kr.X<- kronecker( X[i,,drop=F], array(1, dim=c(m.i,1)) )
			
			if( y.cov.time>0 )	
			{
				y.t.i <- c( Y.time[i, ] );
				y.t.i <- y.t.i[ !is.na(y.t.i) ];
				
				#if ( length( which(is.na(y.t.i)) )> 0) 
				#	y.t.i <- y.t.i[ -which(is.na(y.t.i)) ];

				if(y.cov.time >= 1)
					kr.X<- cbind( kr.X, y.t.i);
				if(y.cov.time >= 2)
					kr.X<- cbind( kr.X, 0.5 * (3*y.t.i^2 - 1) );
				if(y.cov.time >= 3)
					kr.X<- cbind( kr.X, 0.5 * (5*y.t.i^3 - 3*y.t.i) );
				if(y.cov.time >= 4)
					kr.X<- cbind( kr.X, 1/8 * (35*y.t.i^4 - 30*y.t.i^2 + 3 ) );

				#for(t in 1:y.cov.time)
				#	kr.X<- cbind(kr.X, y.t.i^t);
			}
			
			kr.Z<- kronecker( Z[i,,drop=F], array(1, dim=c( m.i, 1 )) )
			W0 <- W0 + t(kr.Z) %*% V.j_x[[i]] %*% kr.Z;
			W1 <- W1 + t(kr.Z) %*% V.j_x[[i]] %*% kr.X;
			W2 <- W2 + t(kr.X) %*% V.j_x[[i]] %*% kr.X;
			W3 <- W3 + t(kr.X) %*% V.j_x[[i]] %*% kr.Z;
		}
		
		Q.w <- W0 - W1 %*% solve(W2) %*% W3;
	}

#cat("v=", Q.v, "u=", Q.u, "\n");
#show(Q.w/2 );

	return(list(v=Q.v, w=Q.w/2, u=Q.u));
}

est_gen_Q<-function(Y.delt, Y.time, Z.scale.skat, Z.scale.burden, X, par_null, y.cov.time, run.cpp=T)
{
	t0 <- proc.time()

	r0 <- NA;
	if(!run.cpp)
	{
		r.skat<- est_gen_Q.R( Y.delt, 
						Y.time, 
						Z.scale.skat$maf, 
						Z.scale.skat$new, 
						X, 
						par_null, 
						y.cov.time);

		r.burden<- est_gen_Q.R( 
						Y.delt, 
						Y.time, 
						Z.scale.burden$maf, 
						Z.scale.burden$new, 
						X, 
						par_null, 
						y.cov.time)

		r <- list( v.skat = r.skat$v, 
						w.skat   = r.skat$w, 
						v.burden = r.burden$u, 
						w.burden = r.burden$w);
		
	}	
	else	
		r0<- .Call( "est_gen_Q_C", 
						as.matrix( Y.delt*1.0 ), 
						as.matrix( Z*1.0 ),  
						as.matrix( X*1.0 ), 
						as.vector(maf*1.0), 
						as.vector( par_null*1.0) );

	t1 <- proc.time();

	# r1<- est_gen_Q.R(Y.delt, maf, Z, X, par_null)
	# if( any(c(r1$v) != c(r0$v) ) ) cat("V1<>V0");
	# if( any(c(r1$w) != c(r0$w) ) ) cat("W1<>W0");
	# cat("EST_GEN_Q Runtime:", (t1-t0)[3], "\n");
	
	return(r);		
}

#public:
longskat_gene_run<-function( r.model, snp.mat, weights.common=c(0.5,0.5), weights.rare=c(1,25), run.cpp=T, rare.cutoff=NULL, r.corr=2, debug=F ) 
{
	Y.delt   <- r.model$y.delt; 
	Y.time   <- r.model$y.time; 
	par_null <- c( r.model$par$sig_a,  r.model$par$sig_b, r.model$par$sig_e, r.model$par$rho );
	X <- cbind(1, r.model$y.cov);

	## SKAT Mode
 	Z.scale.skat <- SKAT_Scale_Genotypes( X, 
 					snp.mat, 
 					weights.common = weights.common, 
 					weights.rare   = weights.rare, 
 					rare.cutoff    = rare.cutoff, 
 					r.corr.common  = 0, 
 					r.corr.rare    = 0 )

	## Burden Mode
 	Z.scale.burden <- SKAT_Scale_Genotypes( X, 
 					snp.mat, 
 					weights.common = weights.common, 
 					weights.rare   = weights.rare, 
 					rare.cutoff    = rare.cutoff, 
 					r.corr.common  = 1, 
 					r.corr.rare    = 1 )

	Q <- try( est_gen_Q( Y.delt, 	
					Y.time, 
					Z.scale.skat, 
					Z.scale.burden, 
					X, 
					par_null, 
					r.model$y.cov.time, 
					run.cpp = run.cpp ) );

	if (class(Q)=="try-error")
		return(NULL);
		
	p.lskat <- get_Qv_pvalue(Q$v.skat, Q$w.skat);
	p.burden <- get_Qu_pvalue(Q$v.burden, Q$w.burden);

	if(debug) cat( " MAF=", (Z.scale.skat$maf), "\n");
	
	r.lskat <- list( mle      = r.model, 
					snp.total = length(Z.scale.skat$maf), 
					snp.rare  = Z.scale.skat$rare, 
					q.lskat    = Q$v.skat, 
					p.lskat    = p.lskat$p.value, 
					q.burden  = Q$v.burden, 
					p.burden  = p.burden$p.value);
	class(r.lskat) <- "LSKAT.gen.ret";
	
	return(r.lskat)
}

#public:
print.LSKAT.gen.ret<-function(r.lskat, useS4=FALSE)
{
	cat("  LSKAT/Gene Summary\n");	
	cat("[1] MLE Results:\n");	
	cat("* SIGMA_A =", 		   r.lskat$mle$par$sig_a, "\n");
	cat("* SIGMA_B =",         r.lskat$mle$par$sig_b, "\n");
	cat("* SIGMA_E =",         r.lskat$mle$par$sig_e, "\n");
	cat("* RHO =",             r.lskat$mle$par$rho, "\n");
	cat("* MU =",              r.lskat$mle$par$mu, "\n");
	cat("* Beta(Cov)=",        r.lskat$mle$par$par_cov, "\n");
	cat("* Beta(Time)=",       r.lskat$mle$par$par_t, "\n");
	cat("* L(min) =",          r.lskat$mle$likelihood, "\n");

	cat("[2] LSKAT:\n");	
	cat("* SNP total = ",      r.lskat$snp.total, "\n");	
	cat("* Rare SNPs = ",      r.lskat$snp.rare, "\n");	
	cat("* Q @ SKAT  = ",      r.lskat$q.lskat, "\n");	
	cat("* p-value   = ",      r.lskat$p.lskat, "\n");	
	cat("* Q @ burden= ",      r.lskat$q.burden, "\n");	
	cat("* p-value   = ",      r.lskat$p.burden, "\n");	
}

#private
longskat_gene_task<-function(r.model, gene.range, PF, gen.list, weights.common, weights.rare, run.cpp, snp.impute, rare.cutoff, debug) 
{
	rs.name <-c();
	rs.lst <- vector("list", length(which(!is.na(gene.range))) );
	rs.lst.idx <- 0;

	if(debug) cat("Shrinking genotype data...\n");

	snpmat <- shrink_snpmat(PF$snp.mat, gen.list, gene.range );
	if(is.null(snpmat))
		stop("No SNP selected from PLINK files.");

	if(debug) cat("Shrinked genotype data!\n");

	for(i in gene.range )
	{
		if (is.na(i)) next;
		
		gen <- get_gen_group( gen.list, i );
		if (debug) cat("  Finding", length(gen$snps), "SNPs...\n");

		gen.mat <- try( get_snp_mat( snpmat, gen, snp.impute ) );
		if( is.null(gen.mat) || class(gen.mat)=="try-error" || length(gen.mat$maf)==0 ) 
		{
			if (debug) cat("! No SNPS for Gene[", i, "]=", as.character(gen$name), "\n");
			next;
		}
	
		#gen.mat$snp:[P,N]=>[N,P]
		ls <- try( longskat_gene_run( r.model, 
						t(gen.mat$snp), 
						weights.common = weights.common, 
						weights.rare   = weights.rare, 
						run.cpp        = run.cpp, 
						rare.cutoff    = rare.cutoff, 
						debug          = debug) );

		rs.lst.idx <- rs.lst.idx + 1;
		rs.name <- c(rs.name, as.character(gen$name));
		
		if(is.null(ls) || class(ls) == "try-error" )
		{
			rs.lst[[rs.lst.idx]] <- c( i, 
						min(gen.mat$info[,2]), 
						min(gen.mat$info[,3]), 
						NROW(gen.mat$mat), 
						NA, NA, NA, NA, NA );
			cat("! Failed to calculate Gene[", i, "]=", as.character(gen$name), "\n");
		}	
		else
		{
			rs.lst[[rs.lst.idx]] <- c(i, 
						min(gen.mat$info[,2]), 
						min(gen.mat$info[,3]), 
						ls$snp.total, 
						ls$snp.rare, 
						ls$q.lskat, 
						ls$p.lskat, 
						ls$q.burden, 
						ls$p.burden );
						
			if(debug) cat(" [", i, "]", 
						as.character(gen$name), 
						ls$snp.total,
						ls$snp.rare, 
						ls$q.lskat, 
						ls$p.lskat, 
						ls$q.burden, 
						ls$p.burden, "\n");
		}
	}	

	rs <- do.call( "rbind", rs.lst );
	ret <- data.frame( id    = rs[,1], 
					name     = rs.name, 
					chr      = rs[,2], 
					min.pos  = rs[,3], 
					snp      = rs[,4], 
					rare     = rs[,5], 
					Q.lskat  = rs[,6], 
					p.lskat  = rs[,7], 
					Q.burden = rs[,8], 
					p.burden = rs[,9] );

	return(ret);
}

#public:
longskat_gene_plink<-function( file.plink.bed, file.plink.bim, file.plink.fam, 
				file.gene.set, 
				file.phe.long, 
				file.phe.cov, 
				file.phe.time=NULL, 
				gene.range = NULL,  
				options=list() )
{
	cat( "[ LONGSKAT_GENE_PLINK ] Procedure.\n");
	cat( "Checking the optional items......\n");

	if (missing(options)) 
		options <- get_default_options()
	else	
	{
		options0 <- get_default_options();
		options0[names(options)] <- options;
		options <- options0;
	}
	
	show_options( options );

	chk.genset <- check_geneset_file( file.gene.set )
	if ( !chk.genset$bSuccess )
		stop("Gene defintion file  is failed to load.")

	chk.phe <- check_pheno_file( file.phe.long, file.phe.time, file.plink.fam ) 
	if ( !chk.phe$bSuccess )
		stop("Phenotypic data file is failed to load.")

	chk.cov <- check_covariate_file( file.phe.cov, file.plink.fam, options$y.cov.count )
	if ( !chk.cov$bSuccess )
		stop("Covariate data file is failed to load.")

	chk.plink <- check_plink_file( file.plink.bed, file.plink.bim, file.plink.fam )
	if ( !chk.plink$bSuccess )
		stop("PLINK file can not be loaded by the snpStats package.")

	cat( "Starting to load all data files......\n");
	
	PF <- read_gen_phe_cov ( file.plink.bed, file.plink.bim, file.plink.fam, 
					file.phe.long, 
					file.phe.time, 
					file.phe.cov ); 
	
	PF.par <- list(file.plink.bed = file.plink.bed, 
					file.plink.bim = file.plink.bim, 
					file.plink.fam = file.plink.fam, 
					file.phe.long  = file.phe.long, 
					file.phe.time  = file.phe.time, 
					file.phe.cov   = file.phe.cov,
					file.gene.set  = file.gene.set,
					y.cov.count    = options$y.cov.count, 
					y.cov.time     = options$y.cov.time,
					weights.common = options$weights.common, 
					weights.rare   = options$weights.rare,
					rare.cutoff    = options$rare.cutoff);

	if( is.na(options$y.cov.count) )
		options$y.cov.count <- NCOL(PF$phe.cov)-1;

	cat( "Starting to estimate the SIGMA_A, SIGMA_B, SIGMA_E and other parameters......\n");

	r.model <- longskat_est_model(PF$phe.long, 
					PF$phe.cov, 
					PF$phe.time, 
					y.cov.time = options$y.cov.time, 
					g.maxiter  = options$g.maxiter, 
					debug      = options$debug);
	
	if( class(r.model) != "LSKAT.null.model" ) 
		stop("! Failed to estimate the parameters of Covariance Compoment.");

	cat("* SIGMA_A =",   r.model$par$sig_a, "\n");
	cat("* SIGMA_B =",   r.model$par$sig_b, "\n");
	cat("* SIGMA_E =",   r.model$par$sig_e, "\n");
	cat("* RHO =",       r.model$par$rho, "\n");
	cat("* MU =",        r.model$par$mu, "\n");
	cat("* Beta(Cov) =", r.model$par$par_cov, "\n");
	cat("* Beta(Time) =",r.model$par$par_t, "\n");
	cat("* L(min) =",    r.model$likelihood, "\n");

	gen.list <- read_gen_dataset( file.gene.set, file.plink.bim );
	
	if (is.null(gene.range)) 
		gene.range<- c(1:gen.list$len)
		
	if( length(which( gene.range > gen.list$len))>0 )
	{
		warning("The gene range is out of data set.");
		gene.range <- gene.range[- which( gene.range > gen.list$len ) ];
	}

	cat("* GENE.RANGE =", min(gene.range),"-", max(gene.range), "[", length(gene.range), "]\n");

	tm.start <- proc.time();
	cpu.fun<-function( sect )
	{
		g.range0 <- gene.range[((sect-1)*n.percpu+1):(sect*n.percpu)];
		
		PF <- read_gen_phe_cov ( PF.par$file.plink.bed, 
						PF.par$file.plink.bim, 
						PF.par$file.plink.fam, 
						PF.par$file.phe.long, 
						PF.par$file.phe.time, 
						PF.par$file.phe.cov );
		
		res.cluster <- longskat_gene_task( r.model, 
						g.range0, 
						PF, 
						gen.list, 
						options$weights.common, 
						options$weights.rare, 
						options$run.cpp, 
						options$snp.impute, 
						options$rare.cutoff, 
						options$debug );
		
		return(res.cluster);
	}

	lskat.ret<-c();
	if( options$n.cpu>1 && require(snowfall) )
	{
		cat("Starting parallel computing, snowfall/snow......\n"); 
		sfInit(parallel=TRUE, cpus=options$n.cpu, type="SOCK")

		n.percpu <- ceiling( length(gene.range)/options$n.cpu );
		sfExport("n.percpu", "gene.range", "r.model", "gen.list", "options", "PF.par" );
		
		res.cluster <- sfClusterApplyLB( 1:options$n.cpu, cpu.fun);
		sfStop();

		cat("Stopping parallel computing......\n"); 

		lskat.ret <- do.call("rbind", res.cluster);
	}
	else
	{
		cat("Starting the LSKAT estimate for each gene......\n");
		lskat.ret <- longskat_gene_task( r.model, 
						gene.range, 
						PF, 
						gen.list, 
						options$weights.common, 
						options$weights.rare, 
						options$run.cpp, 
						options$snp.impute, 
						options$rare.cutoff, 
						options$debug );
	}
	
	tm <- proc.time() - tm.start;
	cat( "* RUNTIME =", tm[3], "seconds \n");

	ret = list( gene=lskat.ret, par=PF.par, mle=r.model);
	class(ret) <- "LSKAT.gen.plink";

	return(ret);
}

#public:
print.LSKAT.gen.plink<-function(r.lskat, useS4=FALSE)
{
	summary.LSKAT.gen.plink(r.lskat);
}

#public:
summary.LSKAT.gen.plink<-function(r.lskat)
{
	cat("  LSKAT/Gene Summary\n");	
	cat("[1] Paramaters:\n");	
	cat("* PHE.LOG.FILE = ",   r.lskat$par$file.phe.long, "\n");	
	cat("* PHE.LOG.TIME = ",   r.lskat$par$file.phe.time, "\n");	
	cat("* PHE.COV.FILE = ",   r.lskat$par$file.phe.cov, "\n");	
	cat("* Covariate Count = ",r.lskat$par$y.cov.count, "\n");	
	cat("* Time Effect = ",    r.lskat$par$y.cov.time, "\n");	
	cat("* Weight Rare = ",    r.lskat$par$weights.rare, "\n");	
	cat("* Weight Common = ",  r.lskat$par$weights.common, "\n");	

	cat("[2] MLE Results:\n");	
	cat("* SIGMA_A =", 		   r.lskat$mle$par$sig_a, "\n");
	cat("* SIGMA_B =",         r.lskat$mle$par$sig_b, "\n");
	cat("* SIGMA_E =",         r.lskat$mle$par$sig_e, "\n");
	cat("* RHO =",             r.lskat$mle$par$rho, "\n");
	cat("* MU =",              r.lskat$mle$par$mu, "\n");
	cat("* Beta(Cov)=",        r.lskat$mle$par$par_cov, "\n");
	cat("* Beta(Time)=",       r.lskat$mle$par$par_t, "\n");
	cat("* L(min) =",          r.lskat$mle$likelihood, "\n");

	#r.lskat$gene
	# [1] id 
	# [2] name
	# [3] chr
	# [4] min.pos
	# [5] SNP count
	# [6] Rare Counts
	# [7] Q
	# [8] p-value
	df.gene <- cbind( r.lskat$gene[, c(2,3,4,5,6,8),drop=F], bonferroni=p.adjust( r.lskat$gene[, 8], method="bonferroni"));
	df.gene <- df.gene[ order(df.gene[,6]), ];
	colnames(df.gene) <- c("Gene", "Chr.", "Pos.", "Total", "Rare", "p-value", "Bonferroni");
	
	n.sig <- c();
	n.sig.i <- c();
	for( i in 20:4 )
	{
		n <- length( which( df.gene$bonferroni <10^(-i)) );
		if(n>0)
		{
			n.sig   <- c(n.sig, n);
			n.sig.i <- c(n.sig.i, i);
		}	
	}

	cat("[3] Longitudinal SKAT Results:\n");
	cat("* Gene Count = ", NROW(r.lskat$gene), "\n");	
	cat("* Total SNP = ", sum(r.lskat$gene[,5]), "\n");	
	cat("* Rare SNP = ",  sum(r.lskat$gene[,6]), "\n");	
	
	if (length(n.sig.i)>0)
	{
		cat("* Significant Genes (Bonferroni Correct):\n")	
		for(i in 1:length(n.sig.i))
			cat("* <= 10^(-", n.sig.i[i], ") Levels:", n.sig[i], "\n")	
	}
	
	cat("* Top 20 Genes (Bonferroni Correct):\n");
	show(df.gene[c(1:20),]);
}

#public:
plot.LSKAT.gen.plink<-function( r.lskat, pdf.file=NA, title="", y.max=NA, bonferroni=F )
{
	if(is.na(pdf.file)) 
		pdf.file<-paste(r.lskat$par$file.phe.long, ".pdf", sep="");

	#CHR, POS, NAME, PV
	df.gene <- r.lskat$gene[,c(3,4,2,8)];
	df.gene <- df.gene[ order(df.gene[,1], df.gene[,2]), ]; 

	if(bonferroni) df.gene[,4] <- df.gene[,4] * NROW(df.gene);

	pdf( pdf.file, width=6, height=3.5)

	par(mar=c(3.2, 3.5,  2.5, 1))
	par(mgp=c(1.4, 0.3, 0), tck=-0.03)
	
	draw_manhattan( df.gene[,c(1,3,4)], 
					map.title = title, 
					ifelse(bonferroni, 0.0001, 0.05)/NROW(df.gene), 
					0.7, 
					y.max = y.max );

	dev.off();

	cat("* Manhattan plot is save into ", pdf.file, "\n");	
}


get_default_options<-function()
{
 	options <- list( rare.cutoff   = NULL, 
					y.cov.count   = NA, 
					y.cov.time    = 0, 
					g.maxiter     = 20, 
					weights.common= c(0.5,0.5), 
					weights.rare  = c(1,25), 
					run.cpp       = T, 
					debug         = F, 
					n.cpu         = 1, 
					snp.impute    = "mean" );

	return(options);	
}

show_options<-function(options)
{
	cat( "* Covariate Count: ",  options$y.cov.count,  "\n");
	cat( "* Covariate Time Effect: ",  options$y.cov.time, "\n");
	cat( "* Parallel Computing: ", ifelse( options$n.cpu>1, "Yes,", "No,"), options$n.cpu,  "CPU(s)\n");
	cat( "* Debug Output: ", ifelse( options$debug, "Yes", "No"),"\n");
	cat( "* SNP Impute: ",  options$snp.impute, "\n");
	cat( "* C/C++ Module Used Output: ", ifelse( options$run.cpp, "Yes", "No"), "\n");
	cat( "* Beta Weights for Common SNPs: ",  options$weights.common[1], options$weights.common[2], "\n");
	cat( "* Beta Weights for Rare SNPs: ",  options$weights.rare[1], options$weights.rare[2], "\n");
	cat( "* Common-Rare cutoff:", options$rare.cutoff, "\n");
}


  