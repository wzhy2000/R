library(mvtnorm);

est_gen_Q.R<-function(Y.delt, Y.t, maf, Z, X, par_null, time.effect)
{
	sig_a2  <- par_null[1]^2
	sig_b2  <- par_null[2]^2
	sig_e2  <- par_null[3]^2
	par_rho <- par_null[4]

	n <- dim(Y.delt)[1];
	m <- dim(Y.delt)[2];
	k <- length(maf);
	
	Q.v<-0;
	Q.w<-c();

	AR.1 <- matrix(rep(1:m,m),ncol=m, byrow=T)
	AR.1 <- par_rho^abs(t(AR.1)-AR.1)

	V.j <- array(1, dim=c(m,m)) * sig_a2 +  AR.1 * sig_b2 + diag(1, m) * sig_e2
	V.j_1 <- solve(V.j);

	if(0)
	#if (dim(Z)[2]<13)
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

		for(i in 1:n)
		{
			y.i <- Y.delt[i,];
			y.j_i <- V.j[!is.na(y.i), !is.na(y.i)];
			V.j_x[[i]] <- solve(y.j_i);
			Y.j_x[[i]] <- y.i[!is.na(y.i)];
			M.j_x[i] <- length(Y.j_x[[i]]);
		}	

		for(i in 1:k)
		{
			Q.i <- 0;
			for(j in 1:n)
				Q.i <- Q.i + Z[j,i]*t(rep(1, M.j_x[j] )) %*% V.j_x[[j]] %*% (Y.j_x[[j]])
			Q.v <- Q.v + (Q.i)^2;
		}

		Q.v <- Q.v/2;

		n.x <- NCOL(X) + time.effect ;
		W0 <- array(0, dim=c(k,k));
		W1 <- array(0, dim=c(k,n.x));
		W2 <- array(0, dim=c(n.x,n.x));
		W3 <- array(0, dim=c(n.x,k));
		for(i in 1:n)
		{
			m.i <- M.j_x[i] ;

			kr.X<- kronecker( X[i,,drop=F], array(1, dim=c(m.i,1)) )
			
			if( time.effect>0 )	
			{
				t.i <- c(Y.t[i, ]);
				if ( length( which(is.na(t.i)) )> 0) 
					t.i <- t.i[ -which(is.na(t.i)) ];

				for(t in 1:time.effect)
					kr.X<- cbind(kr.X, t.i^t);
			}
			
			kr.Z<- kronecker( Z[i,,drop=F], array(1, dim=c( m.i, 1 )) )
			W0 <- W0 + t(kr.Z) %*% V.j_x[[i]] %*% kr.Z;
			W1 <- W1 + t(kr.Z) %*% V.j_x[[i]] %*% kr.X;
			W2 <- W2 + t(kr.X) %*% V.j_x[[i]] %*% kr.X;
			W3 <- W3 + t(kr.X) %*% V.j_x[[i]] %*% kr.Z;
		}
		
		Q.w <- W0 - W1 %*% solve(W2) %*% W3;
	}
	
	return(list(v=Q.v, w=Q.w/2));
}

est_gen_Q<-function(Y.delt, Y.t, maf, Z, X, par_null, time.effect, run.cpp=T)
{
	r0 <- NA;

	t0 <- proc.time()

	if(!run.cpp)
		r0<- est_gen_Q.R( Y.delt, Y.t, maf, Z, X, par_null, time.effect)
	else	
		r0<- .Call( "est_gen_Q_C", as.matrix( Y.delt*1.0 ) , as.matrix( Z*1.0 ),  as.matrix( X*1.0 ), as.vector(maf*1.0), as.vector( par_null*1.0) );

	t1 <- proc.time();

	# r1<- est_gen_Q.R(Y.delt, maf, Z, X, par_null)
	# if( any(c(r1$v) != c(r0$v) ) ) cat("V1<>V0");
	# if( any(c(r1$w) != c(r0$w) ) ) cat("W1<>W0");
	# cat("EST_GEN_Q Runtime:", (t1-t0)[3], "\n");
	
	return(r0);		
}

longskat_gene_run<-function(Y.delt, Y.t, X, snp.mat, par_null, time.effect, weights.rare, weights.common, run.cpp=T, debug=F ) 
{
	Z.scale <- SKAT_Scale_Genotypes(X, t(snp.mat), weights.rare=weights.rare, weights.common=weights.common )

	Q <- try( est_gen_Q( Y.delt, Y.t, Z.scale$maf, Z.scale$new, X, par_null, time.effect, run.cpp=run.cpp ) );

	if (class(Q)=="try-error")
		return(NULL);
		
	P <- get_Q_pvale(Q$v, Q$w);

	if(debug) cat( " MAF=", (Z.scale$maf), "\n");
	
	return( list(snp.total=length(Z.scale$maf), snp.rare=Z.scale$rare, qv=Q$v, pv=P$p.value ) )
}

longskat_gene_task<-function(gene.range, PF, Y.delt, Y.t, X, par_null, gen.lst, time.effect, weights.common, weights.rare, run.cpp, debug) 
{
	rs.name <-c();
	rs.lst <- vector("list", length(which(!is.na(gene.range))) );
	rs.lst.idx <- 0;

	for(i in gene.range )
	{
		if (is.na(i)) next;
		
		gen <- get_gen_group( gen.lst, i );
		if (debug) cat(" Finding", length(gen$snps), "SNPs...\n");

		gen.mat <- get_snp_mat( PF$snp.mat, gen );
		if(is.null(gen.mat)) next;

		if (length(which(gen.mat$maf>0.5))>0) browser();

		#snp.filter <- which( gen.mat$maf > 0.025 | gen.mat$maf < 0.01 );
		#if(length(snp.filter)>0)
		#{
		#	gen.mat$maf <- gen.mat$maf[ -snp.filter ];
		#	gen.mat$snp <- gen.mat$snp[ -snp.filter, ,drop=F];
		#	cat("Remove", length(snp.filter), "SNP due to 0.025-0.01\n");
		#}	

		if (length(gen.mat$maf)==0) next;

		ls <- longskat_gene_run(Y.delt, Y.t, X, gen.mat$snp, par_null, 
				    time.effect = time.effect,	
				    weights.rare=weights.rare, 
				    weights.common=weights.common, 
				    run.cpp=run.cpp,
				    debug = debug);

		
		if(is.null(ls))
		{
			cat("! Failed to calculate Gene[", i, "]=", as.character(gen$name), "\n");
			next;
		}	
		
		rs.lst.idx <- rs.lst.idx + 1;
		rs.lst[[rs.lst.idx]] <- c(i, min(gen.mat$info[,2]), min(gen.mat$info[,3]), ls$snp.total, ls$snp.rare, ls$qv, ls$pv );
		rs.name <- c(rs.name, as.character(gen$name));

		if(debug) cat(" [", i, "]", as.character(gen$name), ls$snp.total,ls$snp.rare, ls$qv, ls$pv, "\n");
	}	

	rs <- do.call( "rbind", rs.lst );
	ret <- data.frame(id=rs[,1], name=rs.name, chr=rs[,2], min.pos=rs[,3], snp=rs[,4], rare=rs[,5], Q=rs[,6], pv = rs[,7] );

	return(ret);
}

longskat_gene_plink<-function( file.plink.bed, file.plink.bim, file.plink.fam, file.phe.long, file.phe.time=NULL, file.phe.cov, file.gene.set, gene.range = NULL,  
	options=list(g.cov= 2, g.ntime= 0, g.maxiter = 20, w.common=c(0.5,0.5), w.rare=c(1,25), run.cpp=T, debug=F, n.cpu=1) )
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
    	
	cat( "* Covariate Count: ",  options$g.cov,  "\n");
	cat( "* Covariate Time Effect: ",  options$g.ntime, "\n");
	cat( "* Parallel Computing: ", ifelse( options$n.cpu>1, "Yes,", "No,"), options$n.cpu,  "CPU(s)\n");
	cat( "* Debug Output: ", ifelse( options$debug, "Yes", "No"),"\n");
	cat( "* C/C++ Module Used Output: ", ifelse( options$run.cpp, "Yes", "No"), "\n");
	cat( "* Beta Weights for Common SNPs: ",  options$w.common[1], options$w.common[2], "\n");
	cat( "* Beta Weights for Rare SNPs: ",  options$w.rare[1], options$w.rare[2], "\n");

	chk.plink <- check_plink_file( file.plink.bed, file.plink.bim, file.plink.fam )
	if ( !chk.plink$bSuccess )
		stop("PLINK file can not be loaded by the snpStats package.")

	chk.phe <- check_pheno_file( file.phe.long, file.phe.time, chk.plink$family ) 
	if ( !chk.phe$bSuccess )
		stop("Phenotypic data file is failed to load.")

	chk.cov <- check_covariate_file( file.phe.cov, chk.plink$family, options$g.cov )
	if ( !chk.cov$bSuccess )
		stop("Covariate data file is failed to load.")

	chk.genset <- check_geneset_file( file.gene.set )
	if ( !chk.genset$bSuccess )
		stop("Gene defintion file  is failed to load.")

	cat( "Starting to load all data files......\n");
	
	PF <- read_gen_phe_cov ( file.plink.bed, file.plink.bim, file.plink.fam, file.phe.long, file.phe.time, file.phe.cov ); 
	PF.pars <- list(file.plink.bed = file.plink.bed, 
			file.plink.bim = file.plink.bim, 
			file.plink.fam = file.plink.fam, 
			file.phe.long  = file.phe.long, 
			file.phe.time  = file.phe.time, 
			file.phe.cov   = file.phe.cov,
			file.gene.set  = file.gene.set,
			g.cov          = options$g.cov, 
			g.ntime        = options$g.ntime,
			w.common       = options$w.common, 
			w.rare         = options$w.rare );

	Y.ncol <- NCOL(PF$phe.long);

	Y <- matrix( as.integer( as.matrix( PF$phe.long[,-1] )), ncol=Y.ncol-1)
	Y.cov <- PF$phe.cov[, c( 3:(2 + options$g.cov))];
	Y.t <- NULL;
	if(!is.null(PF$phe.time))
		Y.t <- matrix( as.integer( as.matrix( PF$phe.time[,-1] )), ncol = Y.ncol-1)


	cat( "Starting to estimate the SIGMA_A, SIGMA_B, SIGMA_E and other parameters......\n");

	h0 <- longskat_est_model(Y, Y.t, Y.cov, nTime=options$g.ntime, g.maxiter=options$g.maxiter, debug=options$debug);
	if( !h0$bSuccess ) 
		stop("! Failed to estimate the parameters of Covariance Compoment.");

	cat("* SIGMA_A =", h0$sig_a, "\n");
	cat("* SIGMA_B =", h0$sig_b, "\n");
	cat("* SIGMA_E =", h0$sig_e, "\n");
	cat("* RHO =", h0$rho, "\n");
	cat("* MU =", h0$u, "\n");
	cat("* Beta(Cov) =", h0$par_cov, "\n");
	cat("* Beta(Time) =", h0$par_t, "\n");
	cat("* L(min) =", h0$val, "\n");

	Y.t    <- h0$y.t;
	Y.delt <- h0$y.delt;
	X <- as.matrix(cbind(1, PF$phe.cov[,c( 3:( 2+options$g.cov ))] ));
	par_null <- c( h0$sig_a,  h0$sig_b, h0$sig_e, h0$rho );
	
	gen.lst <- read_gen_dataset( file.gene.set );
	if (is.null(gene.range)) gene.range<- c(1:gen.lst$len)
	if( length(which( gene.range > gen.lst$len))>0 )
	{
		warning("The gene range is out of data set.");
		gene.range <- gene.range[- which( gene.range > gen.lst$len ) ];
	}

	cat("* GENE.RANGE =", min(gene.range),"-", max(gene.range), "[", length(gene.range), "]\n");

	tm.start <- proc.time();
	cpu.fun<-function( sect )
	{
		g.range0 <- gene.range[((sect-1)*n.percpu+1):(sect*n.percpu)];
		
		PF <- read_gen_phe_cov ( PF.pars$file.plink.bed, PF.pars$file.plink.bim, PF.pars$file.plink.fam, PF.pars$file.phe.long, PF.pars$file.phe.time, PF.pars$file.phe.cov );
		
		res.cluster <- longskat_gene_task( g.range0, PF, Y.delt, Y.t, X, par_null, gen.lst, options$g.ntime, options$w.common, options$w.rare, options$run.cpp, options$debug );
		
		return(res.cluster);
	}

	gen.ret<-c();
	if( options$n.cpu>1 && require(snowfall) )
	{
		cat("Starting parallel computing, snowfall/snow......\n"); 
		sfInit(parallel=TRUE, cpus=options$n.cpu, type="SOCK")

		n.percpu <- ceiling( length(gene.range)/options$n.cpu );
		sfExport("n.percpu", "gene.range", "Y.delt", "Y.t", "X", "par_null", "gen.lst", "options", "PF.pars" );
		
		res.cluster <- sfClusterApplyLB( 1:options$n.cpu, cpu.fun);
		sfStop();

		cat("Stopping parallel computing......\n"); 
		for(i in 1:length(res.cluster))
			gen.ret <- rbind(gen.ret, res.cluster[[i]]);
	}
	else
	{
		cat("Starting the LSKAT estimate for each gene......\n");
		gen.ret <- longskat_gene_task(gene.range, PF, Y.delt, Y.t, X, par_null, gen.lst, options$g.ntime, options$w.common, options$w.rare, options$run.cpp, options$debug );
	}
	
	tm <- proc.time() - tm.start;
	cat( "* RUNTIME =", tm[3], "seconds \n");

	
	r.mle <- list(sig_a   = h0$sig_a,
		   sig_b   = h0$sig_b,
		   sig_b   = h0$sig_e,
		   rho     = h0$rho,
		   mu      = h0$u,
		   beta.cov= h0$par_cov, 
		   beta.time= h0$par_t,
		   val     = h0$val );

	ret = list( gene=gen.ret, pars=PF.pars, mle=r.mle);
	class(ret) <- "LSKAT.gen.ret";

	return(ret);
}

summary.LSKAT.gen.ret<-function(rlskat)
{
	cat("  LSKAT/Gene Summary\n");	
	cat("[1] Paramaters:\n");	
	cat("* PHE.LOG.FILE = ",   rlskat$pars$file.phe.long, "\n");	
	cat("* PHE.LOG.TIME = ",   rlskat$pars$file.phe.time, "\n");	
	cat("* PHE.COV.FILE = ",   rlskat$pars$file.phe.cov, "\n");	
	cat("* Covariate Count = ",rlskat$pars$g.cov, "\n");	
	cat("* Time Effect = ",    rlskat$pars$g.ntime, "\n");	
	cat("* Weight Rare = ",    rlskat$pars$w.rare, "\n");	
	cat("* Weight Common = ",  rlskat$pars$w.common, "\n");	

	cat("[2] MLE Results:\n");	
	cat("* SIGMA_A =", rlskat$mle$sig_a, "\n");
	cat("* SIGMA_B =", rlskat$mle$sig_b, "\n");
	cat("* SIGMA_E =", rlskat$mle$sig_e, "\n");
	cat("* RHO =",     rlskat$mle$rho, "\n");
	cat("* MU =",      rlskat$mle$mu, "\n");
	cat("* Beta(Cov)=", rlskat$mle$par_cov, "\n");
	cat("* Beta(Time)=", rlskat$mle$par_t, "\n");
	cat("* L(min) =",  rlskat$mle$val, "\n");

	df.gene <- rlskat$gene[, c(2,3,4,5,6,8),drop=F];
	df.gene <- df.gene[ order(df.gene[,6]),,drop=F ]; 
	df.gene.bonf <- df.gene[,6] * NROW(df.gene);
	
	n.sig <- c();
	n.sig.i <- c();
	for( i in 20:4 )
	{
		n <- length( which( df.gene.bonf <10^(-i)) );
		if(n>0)
		{
			n.sig   <- c(n.sig, n);
			n.sig.i <- c(n.sig.i, i);
		}	
	}

	#rlskat$gene
	# [1] id 
	# [2] name
	# [3] chr
	# [4] min.pos
	# [5] SNP count
	# [6] Rare Counts
	# [7] Q
	# [8] p-value

	cat("[3] Longitudinal SKAT Results:\n");
	cat("* Gene Count = ", NROW(rlskat$gene), "\n");	
	cat("* Total SNP = ", sum(rlskat$gene[,5]), "\n");	
	cat("* Rare SNP = ",  sum(rlskat$gene[,6]), "\n");	
	
	if (length(n.sig.i)>0)
	{
		cat("* Significant Genes (Bonferroni Correct):\n")	
		for(i in 1:length(n.sig.i))
			cat("* <= 10^(-", n.sig.i[i], ") Levels:", n.sig[i], "\n")	
	}
	
	cat("* Top 20 Genes (Bonferroni Correct):\n");
	df.show <- df.gene[c(1:20), c(1,2,3,4,6)];
	df.show <- cbind(df.show, df.show[,5]*NROW(df.gene) )
	colnames(df.show) <- c("Gene", "Chr.", "Pos.", "SNPs", "p-value", "Bonferroni");

	show(df.show);
}

plot.LSKAT.gen.ret<-function( rlskat, pdf.file=NA, title="", y.max=NA, bonferroni=F )
{
	if(is.na(pdf.file)) pdf.file<-paste(rlskat$pars$file.phe.long, ".pdf", sep="");

	#CHR, POS, NAME, PV
	df.gene <- rlskat$gene[,c(3,4,2,8)];
	df.gene <- df.gene[ order(df.gene[,1], df.gene[,2]), ]; 

	if(bonferroni) df.gene[,4] <- df.gene[,4] * NROW(df.gene);

	pdf( pdf.file, width=6, height=3.5)

	par(mar=c(3.2, 3.5,  2.5, 1))
	par(mgp=c(1.4, 0.3, 0), tck=-0.03)
	draw_manhattan( df.gene[,c(1,3,4)], map.title=title, ifelse(bonferroni, 0.0001, 0.05)/NROW(df.gene), 0.7, y.max= y.max );

	dev.off();

	cat("* Manhattan plot is save into ", pdf.file, "\n");	
}


get_default_options<-function()
{
	options <- list(g.cov= 2, g.ntime= 0, g.maxiter = 20, w.common=c(0.5,0.5), w.rare=c(1,25), run.cpp=T, debug=F, n.cpu=1);
	return(options);	
}

longskat_gene_plink_profile<-function( file.plink.bed, file.plink.bim, file.plink.fam, file.phe.long, file.phe.time=NULL, file.phe.cov, file.gene.set, gene.names, options=list()) 
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
	
	cat( "* Covariate Count: ",  options$g.cov,  "\n");
	cat( "* Covariate Time Effect: ",  options$g.ntime, "\n");
	cat( "* Parallel Computing: ", ifelse( options$n.cpu>1, "Yes,", "No,"), options$n.cpu,  "CPU(s)\n");
	cat( "* Debug Output: ", ifelse( options$debug, "Yes", "No"),"\n");
	cat( "* C/C++ Module Used Output: ", ifelse( options$run.cpp, "Yes", "No"), "\n");
	cat( "* Beta Weights for Common SNPs: ",  options$w.common[1], options$w.common[2], "\n");
	cat( "* Beta Weights for Rare SNPs: ",  options$w.rare[1], options$w.rare[2], "\n");

	chk.plink <- check_plink_file( file.plink.bed, file.plink.bim, file.plink.fam )
	if ( !chk.plink$bSuccess )
		stop("PLINK file can not be loaded by the snpStats package.")

	chk.phe <- check_pheno_file( file.phe.long, file.phe.time, chk.plink$family ) 
	if ( !chk.phe$bSuccess )
		stop("Phenotypic data file is failed to load.")

	chk.cov <- check_covariate_file( file.phe.cov, chk.plink$family, options$g.cov )
	if ( !chk.cov$bSuccess )
		stop("Covariate data file is failed to load.")

	chk.genset <- check_geneset_file( file.gene.set )
	if ( !chk.genset$bSuccess )
		stop("Gene defintion file  is failed to load.")

	cat( "Starting to load all data files......\n");
	
	PF <- read_gen_phe_cov ( file.plink.bed, file.plink.bim, file.plink.fam, file.phe.long, file.phe.time, file.phe.cov ); 
	PF.pars <- list(file.plink.bed = file.plink.bed, 
			file.plink.bim = file.plink.bim, 
			file.plink.fam = file.plink.fam, 
			file.phe.long  = file.phe.long, 
			file.phe.time  = file.phe.time, 
			file.phe.cov   = file.phe.cov,
			file.gene.set  = file.gene.set,
			g.cov          = options$g.cov, 
			g.ntime        = options$g.ntime,
			w.common       = options$w.common, 
			w.rare         = options$w.rare );

	Y.ncol <- dim(PF$phe.long)[2];

	Y <- matrix( as.integer( as.matrix( PF$phe.long[,-1] )), ncol=Y.ncol-1)
	Y.cov <- PF$phe.cov[,c(3:(2+options$g.cov))];
	Y.t <- NULL;
	if(!is.null(PF$phe.time))
		Y.t <- matrix( as.integer( as.matrix( PF$phe.time[,-1] )), ncol=Y.ncol-1)
	
	cat( "Starting to estimate the SIGMA_A, SIGMA_B, SIGMA_E and other parameters......\n");

	h0 <- longskat_est_model(Y, Y.t, Y.cov, nTime=options$g.ntime, g.maxiter=options$g.maxiter, debug=options$debug);
	if( !h0$bSuccess ) 
		stop("! Failed to estimate the parameters of Covariance Compoment.");

	cat("* SIGMA_A =", h0$sig_a, "\n");
	cat("* SIGMA_B =", h0$sig_b, "\n");
	cat("* SIGMA_E =", h0$sig_e, "\n");
	cat("* RHO =", h0$rho, "\n");
	cat("* MU =", h0$u, "\n");
	cat("* Beta(Cov) =", h0$par_cov, "\n");
	cat("* Beta(Time) =", h0$par_t, "\n");
	cat("* L(min) =", h0$val, "\n");

	Y.delt <- h0$y.delt;
	Y.t    <- h0$y.time;
	X <- as.matrix(cbind(1, PF$phe.cov[,c(3:(2+options$g.cov))] ));
	par_null <- c(h0$sig_a, h0$sig_b, h0$sig_e, h0$rho);
	
	cat("* GENE.NAMES =", length(gene.names), "-", gene.names, "\n");
	cat("Starting the LSKAT estimate for each gene......\n");

	tm.start <- proc.time();
	
	gen.lst <- read_gen_dataset(file.gene.set);
	
	rs.lst <- list();
	rs.lst.idx <- 0;

	for(k in 1:length(gene.names) )
	{
		if (is.na(k)) next;
		
		gen <- get_gen_family( gen.lst, gene.names[k] );
		if (options$debug) cat(" Finding", length(gen$snps), "SNPs...\n");
		gen.mat <- get_snp_mat( PF$snp.mat, gen );

		if(is.null(gen.mat)) next;
		if (length(which(gen.mat$maf>0.5))>0) browser();

		if (length(gen.mat$maf)==0) next;
		
		ls <- longskat_gene_run(Y.delt, Y.t, X, gen.mat$snp, par_null, 
				    time.effect = options$g.ntime,
				    weights.rare  = options$w.rare, 
				    weights.common= options$w.common, 
				    run.cpp=options$run.cpp,
				    debug = options$debug);

		if(is.null(ls))
		{
			cat("! Failed to calculate Gene[", k, "]=", as.character(gen$name), "\n");
			next;
		}	

		rs.lst.idx <- rs.lst.idx + 1;

		rs.name  <- c(as.character(gen$name));
		rs.lskat <- c(0, min(gen.mat$info[,2]), min(gen.mat$info[,3]), ls$snp.total, ls$snp.rare, ls$qv, ls$pv );

		if(options$debug) cat(" [", k, "]", as.character(gen$name), ls$snp.total, ls$snp.rare, ls$qv, ls$pv, "\n");

		if(NROW(gen.mat$info)>1)
		{
			for(i in 1:NROW(gen.mat$info) )
			{
				snp.mat <- gen.mat$snp[-i,,drop=F]
				ls <- longskat_gene_run(Y.delt, Y.t, X, snp.mat, par_null, 
						    time.effect = options$g.ntime,
						    weights.rare  = options$w.rare, 
						    weights.common= options$w.common, 
						    run.cpp=options$run.cpp,
						    debug = options$debug);

				if(is.null(ls))
				{
					cat("! Failed to calculate SNP[", i, "]=", gen.mat$info[i,1], "\n");
					next;
				}	
	
				rs.name  <- c(rs.name, paste("-", gen.mat$info[i,1], sep=""));
				rs.lskat <- rbind(rs.lskat, c(i, min(gen.mat$info[,2]), min(gen.mat$info[,3]), ls$snp.total, ls$snp.rare, ls$qv, ls$pv ) );
			}
			
			if(NROW(gen.mat$info)>2)
			{
				pv <- c();
				for(i in 1:NROW(gen.mat$info) )
				{
					snp.mat <- gen.mat$snp[i,,drop=F]
					ls <- longskat_gene_run(Y.delt, X, snp.mat, par_null, 
							    time.effect = options$g.ntime,
							    weights.rare  = options$w.rare, 
							    weights.common= options$w.common, 
							    run.cpp=options$run.cpp,
							    debug = options$debug);

					rs.name  <- c(rs.name, paste("*", gen.mat$info[i,1], sep=""));
					rs.lskat <- rbind(rs.lskat, c(i, min(gen.mat$info[,2]), min(gen.mat$info[,3]), ls$snp.total, ls$snp.rare, ls$qv, ls$pv ) );
					pv <- c(pv, ls$pv);
				}
				
				r.pv <- sort.int(pv, decreasing=T, index.return=T);
				r.sel <- c();
				for(i in r.pv$ix )
				{
					r.sel <- c(r.sel, i);
					snp.mat <- gen.mat$snp[r.sel,,drop=F]
					ls <- longskat_gene_run(Y.delt, X, snp.mat, par_null, 
							    time.effect = options$g.ntime,
							    weights.rare  = options$w.rare, 
							    weights.common= options$w.common, 
							    run.cpp=options$run.cpp,
							    debug = options$debug);

					rs.name  <- c(rs.name, paste("*", gen.mat$info[r.sel,1], sep="", collapse="&"));
					rs.lskat <- rbind(rs.lskat, c(i, min(gen.mat$info[,2]), min(gen.mat$info[,3]), ls$snp.total, ls$snp.rare, ls$qv, ls$pv ) );
					pv <- c(pv, ls$pv);
				}

			}
			
			if(0)
			{
				if(NROW(gen.mat$info)>3)
				{
					p <- combn(1:NROW(gen.mat$info), 2);
					for(i in 1:NCOL(p) )
					{
						snp.mat <- gen.mat$snp[ p[,i],, drop=F ]
						ls <- longskat_gene_run(Y.delt, X, snp.mat, par_null, 
								    time.effect = options$g.ntime,
								    weights.rare  = options$w.rare, 
								    weights.common= options$w.common, 
								    run.cpp=options$run.cpp,
								    debug = options$debug);

						rs.name  <- c(rs.name, paste("*", gen.mat$info[p[1,i],1], "&", gen.mat$info[p[2,i],1], sep=""));
						rs.lskat <- rbind(rs.lskat, c(i, min(gen.mat$info[,2]), min(gen.mat$info[,3]), ls$snp.total, ls$snp.rare, ls$qv, ls$pv ) );
					}

					for(i in 1:NCOL(p) )
					{
						snp.mat <- gen.mat$snp[ -p[,i],, drop=F ]
						ls <- longskat_gene_run(Y.delt, X, snp.mat, par_null, 
								    time.effect = options$g.ntime,
								    weights.rare  = options$w.rare, 
								    weights.common= options$w.common, 
								    run.cpp=options$run.cpp,
								    debug = options$debug);

						rs.name  <- c(rs.name, paste("-", gen.mat$info[p[1,i],1], "&", gen.mat$info[p[2,i],1], sep=""));
						rs.lskat <- rbind(rs.lskat, c(i, min(gen.mat$info[,2]), min(gen.mat$info[,3]), ls$snp.total, ls$snp.rare, ls$qv, ls$pv ) );
					}
				}
			}
			
			rownames(rs.lskat)<-rs.name;
			colnames(rs.lskat)<-c("idx", "chr", "min.pos", "snp", "rare", "Q", "pv");
		}			
		
		rs.lst[[rs.lst.idx]] <- rs.lskat;
	}	

	tm <- proc.time() - tm.start;
	cat( "* RUNTIME =", tm[3], "seconds \n");

	return(rs.lst);
}
