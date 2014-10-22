est_snp_Q.R<-function(Y.delt, maf, Z, X, par_null )
{
	sig_a2  <- par_null[1]^2
	sig_b2  <- par_null[2]^2
	sig_e2  <- par_null[3]^2
	par_rho <- par_null[4]

	x.col <- dim(X)[2];
	n <- dim(Y.delt)[1];
	m <- dim(Y.delt)[2];
	k <- 1;

	AR.1 <- array(0,dim=c(m,m));
	for(i in 1:m)
	for(j in 1:m)
		AR.1[i,j] <- par_rho^abs(i-j);

	V.j <- array(1, dim=c(m,m)) * sig_a2 +  AR.1 * sig_b2 + diag(1, m) * sig_e2
	V.j_1 <- solve(V.j);

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

	Q.i <- 0;
	for(i in 1:n)
	#	Q.i <- Q.i + (Z[i,1]*t(rep(1,m))%*%V.j_x[[i]]%*%t(y.delt[i,,drop=F]))
		Q.i <- Q.i + Z[i,1]*t(rep(1, M.j_x[i] ))%*%V.j_x[[i]]%*% ( Y.j_x[[i]] )
		
		
	Q.v1 <- (Q.i)^2/2
	Q.v <- Q.v1[1,1]

	W0 <- array(0, dim=c(k, k));
	W1 <- array(0, dim=c(k, x.col));
	W2 <- array(0, dim=c(x.col, x.col));
	W3 <- array(0, dim=c(x.col, k));
	for(i in 1:n)
	{
		kr <- array(1, dim=c(1, M.j_x[i])) %*% V.j_x[[i]] %*% array(1, dim=c(M.j_x[i],1))
			
		W0 <- W0 + Z[i,1]^2 * kr;
		W1 <- W1 + Z[i,1] * kr %*% t( X[i,] );
		W2 <- W2 + kr[1,1] *  X[i,] %*%t ( X[i,] );	
		W3 <- W3 + Z[i,1] * X[i,] %*% kr;
	}
		
	Q.w1 <- (W0 - W1 %*% solve(W2) %*% W3)/2;
	Q.w <- Q.w1[1,1];
		
	#cat("Q.v==Q.v1", Q.v==Q.v1, "Q.w==Q.w1", Q.w==Q.w1, Q.v, Q.v1, Q.w, Q.w1, "\n");	
	
	r   <- sqrt(Q.v)/Q.w;
	sig_r2 <- 0;
	for(i in 1:n)
	{
		KY <- t(Y.j_x[[i]] - Z[i,1]*r)
		sig_ri <- KY%*%V.j_x[[i]]%*%t(KY);
		sig_r2 <- sig_r2 + sig_ri[1,1];
	}
	
	sig_r2 <- sig_r2/( n - 1- (x.col-1) -1 )
	p.v <-  pchisq(Q.v/(Q.w*sig_r2), df=1, lower.tail=F);
	
	return(list(v=Q.v, w=Q.w, r=r, sig_r2=sig_r2, chi2=Q.v/(Q.w*sig_r2), p.v=p.v ));
}

est_snp_Q<-function(Y.delt, maf, Z, X, par_null, run.cpp=T)
{
	r0 <- NA;
	t0<-Sys.time()

	#if(!run.cpp)
		r0<- est_snp_Q.R(Y.delt, maf, Z, X, par_null)
	#else	
	#	r0<-.Call("est_snp_Q_C", Y.delt, Z,  X, maf, par_null );

	t1<-Sys.time()
	
#if( r1$v!=r0$v ) cat("V1<>V0")
#if( all(r1$w!=r0$w) ) cat("W1<>W0")
#cat("R Time:", t1-t0, "C Time:", t2-t1 , "\n");

	return(r0);		
}

get_weights<-function(maf, n, beta.common=c(0.5, 0.5), beta.rare=c(1, 25) )
{
	rare.cutoff <- 1/sqrt(2*n);
	wj <- c();	
	for( i in 1:length(maf))
	{
		if ( maf[i] > rare.cutoff )
			wj<- c( wj, dbeta( maf[i], beta.common[1], beta.common[2] ) )
		else
			wj<- c( wj, dbeta( maf[i], beta.rare[1], beta.rare[2] ) )
	}

	return(wj);
}

SKAT_Scale_Genotypes_snp= function( Z, weights.rare=c(1,25), weights.common=c(0.5,0.5) )
{
	n<-dim(Z)[1];
	rare.cutoff <- 1/sqrt(2*n);
	Z.maf <- colMeans(Z)/2;
	
	wr <- get_weights(Z.maf, n, beta.common=weights.common, beta.rare=weights.rare);
	Z <- Z * wr;
	
	return( list( new=Z, maf=Z.maf, rare = ifelse(Z.maf<rare.cutoff,1,0) ) )
}

longskat_snp_run<-function( Y.delt, X, snp1, par_null,  weights.rare, weights.common, run.cpp=T, debug=debug)
{
	if (length(snp1$miss)>0)
	{
		if(debug) cat("! Missing SNP:", length(snp1$miss), "\n");
		Y.delt <- Y.delt[ -snp1$miss, ,drop=F];
		X <- X[-snp1$miss, ,drop=F];
	}

	Z.scale <- SKAT_Scale_Genotypes_snp(snp1$snp, weights.rare=weights.rare, weights.common=weights.common )

	Q <- est_snp_Q( Y.delt, Z.scale$maf, Z.scale$new, X, par_null, run.cpp=run.cpp );

	P <- get_Q_pvale(Q$v, Q$w);

	return( list(snp.total=length(Z.scale$maf), snp.rare=Z.scale$rare, qv=Q$v, pv=P$p.value) )
}

longskat_snp_task<-function(snp.range, file.gene.set, PF, Y.delt, X, par_null, weights.rare, weights.common, run.cpp=T, debug=debug)
{
	gen.tb <- read.table(file.gene.set, sep=" ", header=F);

	rs.name <-c();
	gene.name <-c();
	rs.lst <- vector("list", length(which(!is.na(snp.range))) );
	rs.lst.idx <- 0;
	
	for(i in snp.range )
	{
		if (is.na(i)) next;
		
		snp1 <- get_snp_info(i, PF$snp.mat, gen.tb );

		if (any(is.na(snp1)) ) next;
		if (length(which(snp1$maf>0.5))>0) browser();
		if (length(snp1$maf)==0) next;

		ls <- longskat_snp_run( Y.delt, X, snp1, par_null, weights.rare, weights.common, run.cpp, debug = debug);
		
		rs.lst.idx <- rs.lst.idx + 1;
		rs.lst[[rs.lst.idx]] <- c(i, snp1$chr, snp1$loc, snp1$maf, length(snp1$miss), ls$snp.total, ls$snp.rare, ls$qv, ls$pv )
		rs.name     <- c(rs.name, as.character(snp1$name) );
		gene.name   <- c(gene.name, as.character(snp1$gene) );

		if (debug) cat(" [", i, "]", snp1$chr, snp1$loc, as.character(snp1$name), as.character(snp1$gene), snp1$maf, length(snp1$miss), ls$snp.total, ls$snp.rare, ls$qv, ls$pv, "\n");
	}	

	rs <- do.call( "rbind", rs.lst );
	ret <- data.frame(id=rs[,1], chr=rs[,2], loc=rs[,3], name=rs.name, gene=gene.name, maf=rs[,4], miss=rs[,5], total=rs[,6], rare=rs[,7], Q=rs[,8], pv = rs[,9]);

	return(ret);
}

longskat_snp_plink<-function( file.plink.bed, file.plink.bim, file.plink.fam, file.phe.long, file.phe.cov, file.gene.set, snp.range=NA,  
	options=list( g.cov = 2, g.btime= 0, g.maxiter = 20, w.common=c(0.5,0.5), w.rare=c(1,25), run.cpp=T, debug=F, n.cpu=1 ) )
{	
	cat( "[ LONGSKAT_SNP_PLINK ] Procedure\n");
	cat( "Checking the optional items......\n");

	cat( "* Covariate Count: ",  options$g.cov,  "\n");
	cat( "* Covariate Time Effect: ",  options$g.btime, "\n");

	if (is.null(options$n.cpu)   || is.na(options$n.cpu)) options$n.cpu<-1;
	cat( "* Parallel Computing: ", ifelse( options$n.cpu>1, "Yes,", "No,"), options$n.cpu,  "CPU(s)\n");

	if (is.null(options$debug)   || is.na(options$debug)) options$debug<-FALSE;
	cat( "* Debug Output: ", ifelse( options$debug, "Yes", "No"),"\n");

	if (is.null(options$run.cpp) || is.na(options$run.cpp)) options$run.cpp<-TRUE;
	cat( "* C/C++ Module Used Output: ", ifelse( options$run.cpp, "Yes", "No"), "\n");

	if (is.null(options$w.common)|| is.na(options$w.common)) options$w.common<-c(0.5,0.5);
	cat( "* Beta Weights for Common SNPs: ",  options$w.common[1], options$w.common[2], "\n");
	
	if (is.null(options$w.rare)  || is.na(options$w.rare)) options$w.rare<-c(1,25);
	cat( "* Beta Weights for Rare SNPs: ",  options$w.rare[1], options$w.rare[2], "\n");

	chk.plink <- check_plink_file( file.plink.bed, file.plink.bim, file.plink.fam )
	if ( !chk.plink$bSuccess )
		stop("PLINK file can not be loaded by the snpStats package.")

	chk.phe <- check_pheno_file( file.phe.long, chk.plink$family ) 
	if ( !chk.phe$bSuccess )
		stop("Phenotypic data file is failed to load.")

	chk.cov <- check_covariate_file( file.phe.cov, chk.plink$family, options$g.cov )
	if ( !chk.cov$bSuccess )
		stop("Covariate data file is failed to load.")

	chk.genset <- check_geneset_file( file.gene.set )
	if ( !chk.genset$bSuccess )
		stop("Gene defintion file  is failed to load.")

	cat( "Starting to load all data files......\n");

	PF.pars <- list(file.plink.bed = file.plink.bed, 
			file.plink.bim = file.plink.bim, 
			file.plink.fam = file.plink.fam, 
			file.phe.long  = file.phe.long, 
			file.phe.cov   = file.phe.cov,
			file.gene.set  = file.gene.set);

	PF <- read_gen_phe_cov( file.plink.bed, file.plink.bim, file.plink.fam, file.phe.long, file.phe.cov ); 

	y.ncol<-dim(PF$phe.log)[2];
	Y <- matrix( as.integer( as.matrix( PF$phe.log[,c(2:y.ncol)] )), ncol = y.ncol - 1)
	C <- PF$phe.cov[, 3:(2+options$g.cov), drop=F];
	
	cat( "Starting to estimate the SIGMA_A, SIGMA_B, SIGMA_E and other parameters......\n");

	h0 <- longskat_est_model( Y, C, bTime=options$g.btime, g.maxiter=options$g.maxiter, debug=options$debug );
	if( !h0$bSuccess ) 
		stop("! Fialed to estimate the parameters of Covariance Compoment.");

	cat("* SIGMA_A =", h0$sig_a, "\n");
	cat("* SIGMA_B =", h0$sig_b, "\n");
	cat("* SIGMA_E =", h0$sig_e, "\n");
	cat("* RHO =", h0$rho, "\n");
	cat("* MU =", h0$u, "\n");
	cat("* COV =", h0$cov, "\n");
	cat("* L(min) =", h0$val, "\n");
	
	Y.delt <- h0$y.delt;
	X <- as.matrix( cbind(1, C) );
	par_null <- c(h0$sig_a, h0$sig_b, h0$sig_e, h0$rho)

	gen.len <- dim(PF$snp.mat$genotypes)[2];
	if (all(is.na(snp.range)))
		snp.range<- c(1:gen.len)

	cat("* SNP.RANGE =", min(snp.range),"-", max(snp.range), "[", length(snp.range), "]\n");
	
	cpu.fun<-function( sect )
	{
		s.range0 <- snp.range[((sect-1)*n.percpu+1):(sect*n.percpu)];
		
		PF <- read_gen_phe_cov ( PF.pars$file.plink.bed, PF.pars$file.plink.bim, PF.pars$file.plink.fam, PF.pars$file.phe.long, PF.pars$file.phe.cov );
		
		ret.cluster <- longskat_snp_task( s.range0, PF.pars$file.gene.set, PF, Y.delt, X, par_null, weights.rare=options$w.rare, weights.common=options$w.common, run.cpp=options$run.cpp, debug=options$debug );

		return(ret.cluster);
	}

	snp.ret<-c();
	tm.start <- proc.time();
	if( options$n.cpu>1 && require(snowfall) )
	{
		cat("Starting parallel computing, snowfall/snow......\n"); 
		sfInit(parallel=TRUE, cpus=options$n.cpu, type="SOCK")

		n.percpu <- ceiling( length(snp.range)/options$n.cpu );
		sfExport("n.percpu", "snp.range", "Y.delt", "X", "par_null", "options", "PF.pars" );
		
		ret.cluster <- sfClusterApplyLB( 1:options$n.cpu, cpu.fun);
		sfStop();

		cat("Stopping parallel computing......\n"); 
		for(i in 1:length(ret.cluster))
			snp.ret <- rbind( snp.ret, ret.cluster[[i]]);
	}
	else
	{
		cat("Starting the LSKAT estimate for each gene......\n");
		snp.ret <- longskat_snp_task( snp.range, file.gene.set, PF, Y.delt, X, par_null, weights.rare=options$w.rare, weights.common=options$w.common, run.cpp=options$run.cpp, debug=options$debug );
	}
	
	tm <- proc.time() - tm.start;
	cat( "* RUNTIME =", tm[3], "seconds \n");
	
	ret<-list(snp.ret=snp.ret, pars=PF.pars);
	class(ret) <- "LSKAT.snp.ret";

	return(ret);
}

summary.LSKAT.snp.ret<-function(r.snp)
{
	# 1st = Chr, 2nd= Loc, 3rd=SNP.name, 4th=Gene.name, 5th=MAF, 
	# 6th = NMISS, 7th=Total, 8th=Rare, 9th=Q, 10th= PV,
	df.snp <- r.snp$snp.ret[, c(1,2,3,4,5,10),drop=F];
	df.snp <- df.snp[ order(df.snp[,6]),,drop=F ]; 
	df.snp[,6] <- df.snp[,6]*NROW(df.snp);
	
	n.sig <- c();
	n.sig.i <- c();
	for( i in 20:4 )
	{
		n <- length( which(df.snp[,6]<10^(-i)) );
		if(n>0)
		{
			n.sig   <- c(n.sig, n);
			n.sig.i <- c(n.sig.i, i);
		}	
	}

	cat("  Summary of LSKAT result:\n");	
	cat("* PHE.LOG.FILE = ", r.snp$pars$file.phe.long, "\n")	
	cat("* PHE.COV.FILE = ", r.snp$pars$file.phe.cov, "\n")	
	if (length(n.sig.i)>0)
	{
		cat("* Significant Genes (Bonferroni Correct):\n")	
		for(i in 1:length(n.sig.i))
			cat("* <= 10^(-", n.sig.i[i], ") Levels:", n.sig[i], "\n")	
	}
	
	cat("* Top 20 SNPs (Bonferroni Correct):\n")	
	cat("* [", "X", "]", "SNP Name", "\t", "Gene Name", "\t", "Chromosome", "\t", "Position", "\t", "MAF", "\t", "p-Value","\n");	
	for(i in 1:20)
	{
		if( i > NROW(df.snp) ) next;
		cat("* [", i, "]", as.character(df.snp[i,3]), "\t", df.snp[i,4], "\t", df.snp[i,1], "\t", df.snp[i,2], "\t", df.snp[i,5],"\t", df.snp[i,6], "\n");	
	
	}
}

plot.LSKAT.snp.ret<-function( r.snp, pdf.file=NA, title="",  y.max=NA )
{
	if(is.na(pdf.file)) pdf.file<-paste(r.snp$pars$file.phe.long, ".pdf", sep="");

	#CHR, POS, NAME, PV
	df.snp <- r.snp$snp.ret[,c(1,2,3,10)];
	df.snp <- df.snp[ order(df.snp[,1], df.snp[,2]), ]; 
	df.snp[,4] <- df.snp[,4]*NROW(df.snp);

	pdf( pdf.file, width=6, height=3.5)

	par(mar=c(3.2, 3.5,  2.5, 1))
	par(mgp=c(1.4, 0.3, 0), tck=-0.03)
	draw_manhattan( df.snp[,c(1,3,4)], map.title=title, 0.0001/NROW(df.snp), 0.7, y.max= y.max );

	dev.off();

	cat("* Manhattan plot is save into ", pdf.file, "\n");	
}
