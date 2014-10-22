library(mvtnorm);

est_gen_Q.R<-function(y.delt, maf, Z, X, par_null)
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

	AR.1 <- array(0,dim=c(m,m));
	for(i in 1:m)
	for(j in 1:m)
		AR.1[i,j] <- par_rho^abs(i-j);

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

		Q.v <- Q.v/2

		n.x <- dim(X)[2];
		W0 <- array(0, dim=c(k,k));
		W1 <- array(0, dim=c(k,n.x));
		W2 <- array(0, dim=c(n.x,n.x));
		W3 <- array(0, dim=c(n.x,k));
		for(i in 1:n)
		{
			m.i <- M.j_x[i] 
			kr.Z<- kronecker( Z[i,,drop=F], array(1, dim=c(m.i,1)) )
			kr.X<- kronecker( X[i,,drop=F], array(1, dim=c(m.i,1)) )
			W0 <- W0 + t(kr.Z) %*% V.j_x[[i]] %*% kr.Z;
			W1 <- W1 + t(kr.Z) %*% V.j_x[[i]] %*% kr.X;
			W2 <- W2 + t(kr.X) %*% V.j_x[[i]] %*% kr.X;
			W3 <- W3 + t(kr.X) %*% V.j_x[[i]] %*% kr.Z;
		}
		
		Q.w <- W0 - W1 %*% solve(W2) %*% W3;
	}
	
	return(list(v=Q.v, w=Q.w/2));
}

est_gen_Q<-function(Y.delt, maf, Z, X, par_null, run.cpp=T)
{
	r0 <- NA;

	t0<-proc.time()

	if(!run.cpp)
		r0<- est_gen_Q.R(Y.delt, maf, Z, X, par_null)
	else	
		r0<-.Call("est_gen_Q_C", Y.delt, Z,  X, maf, par_null );

	t1<-proc.time()
	
#if( r1$v!=r0$v ) cat("V1<>V0")
#if( all(r1$w!=r0$w) ) cat("W1<>W0")
#cat("EST_GEN_Q Runtime:", (t1-t0)[3], "\n");

	return(r0);		
}

longskat_gene_run<-function(Y.delt, X, snp.mat, par_null, weights.rare, weights.common, run.cpp=T, debug=F ) 
{
	Z.scale <- SKAT_Scale_Genotypes(X, t(snp.mat), weights.rare=weights.rare, weights.common=weights.common )

	Q <- est_gen_Q( Y.delt, Z.scale$maf, Z.scale$new, X, par_null, run.cpp=run.cpp );

	P <- get_Q_pvale(Q$v, Q$w);

	if(debug) cat( " MAF=", (Z.scale$maf), "\n");
	
	return( list(snp.total=length(Z.scale$maf), snp.rare=Z.scale$rare, qv=Q$v, pv=P$p.value) )
}

longskat_gene_task<-function(gene.range, PF, Y.delt, X, par_null, gen.lst, weights.common, weights.rare, run.cpp, debug) 
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

		if(is.na(gen.mat)) next;
		if (length(which(gen.mat$maf>0.5))>0) browser();

		#snp.filter <- which( gen.mat$maf > 0.025 | gen.mat$maf < 0.01 );
		#if(length(snp.filter)>0)
		#{
		#	gen.mat$maf <- gen.mat$maf[ -snp.filter ];
		#	gen.mat$snp <- gen.mat$snp[ -snp.filter, ,drop=F];
		#	cat("Remove", length(snp.filter), "SNP due to 0.025-0.01\n");
		#}	


		if (length(gen.mat$maf)==0) next;

		ls <- longskat_gene_run(Y.delt, X, gen.mat$snp, par_null, 
				    weights.rare=weights.rare, 
				    weights.common=weights.common, 
				    run.cpp=run.cpp,
				    debug = debug);


		rs.lst.idx <- rs.lst.idx + 1;
		rs.lst[[rs.lst.idx]] <- c(i, min(gen.mat$info[,2]), min(gen.mat$info[,3]), ls$snp.total, ls$snp.rare, ls$qv, ls$pv );
		rs.name <- c(rs.name, as.character(gen$name));

		if(debug) cat(" [", i, "]", as.character(gen$name), ls$snp.total,ls$snp.rare, ls$qv, ls$pv, "\n");
	}	

	rs <- do.call( "rbind", rs.lst );
	ret <- data.frame(id=rs[,1], name=rs.name, chr=rs[,2], min.pos=rs[,3], snp=rs[,4], rare=rs[,5], Q=rs[,6], pv = rs[,7] );

	return(ret);
}

longskat_gene_plink<-function( file.plink.bed, file.plink.bim, file.plink.fam, file.phe.long, file.phe.cov, file.gene.set, gene.range = NA,  
	options=list(g.cov= 2, g.btime= 0, g.maxiter = 20, w.common=c(0.5,0.5), w.rare=c(1,25), run.cpp=T, debug=F, n.cpu=1) )
{
	cat( "[ LONGSKAT_GENE_PLINK ] Procedure.\n");
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
	
	PF <- read_gen_phe_cov ( file.plink.bed, file.plink.bim, file.plink.fam, file.phe.long, file.phe.cov ); 
	PF.pars <- list(file.plink.bed = file.plink.bed, 
			file.plink.bim = file.plink.bim, 
			file.plink.fam = file.plink.fam, 
			file.phe.long  = file.phe.long, 
			file.phe.cov   = file.phe.cov,
			file.gene.set  =  file.gene.set);

	y.ncol <- dim(PF$phe.log)[2];

	y <- matrix( as.integer( as.matrix( PF$phe.log[,c(2:y.ncol)] )), ncol=y.ncol-1)

	cat( "Starting to estimate the SIGMA_A, SIGMA_B, SIGMA_E and other parameters......\n");

	h0 <- longskat_est_model(y, PF$phe.cov[,c(3:(2+options$g.cov))], bTime=options$g.btime, g.maxiter=options$g.maxiter, debug=options$debug);
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
	X <- as.matrix(cbind(1, PF$phe.cov[,c(3:(2+options$g.cov))] ));
	par_null <- c(h0$sig_a, h0$sig_b, h0$sig_e, h0$rho);
	
	gen.lst <- read_gen_dataset(file.gene.set);
	if (is.na(gene.range)) gene.range<- c(1:gen.lst$len)
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
		
		PF <- read_gen_phe_cov ( PF.pars$file.plink.bed, PF.pars$file.plink.bim, PF.pars$file.plink.fam, PF.pars$file.phe.long, PF.pars$file.phe.cov );
		
		res.cluster <- longskat_gene_task( g.range0, PF, Y.delt, X, par_null, gen.lst, options$w.common, options$w.rare, options$run.cpp, options$debug );
		
		return(res.cluster);
	}

	gen.ret<-c();
	if( options$n.cpu>1 && require(snowfall) )
	{
		cat("Starting parallel computing, snowfall/snow......\n"); 
		sfInit(parallel=TRUE, cpus=options$n.cpu, type="SOCK")

		n.percpu <- ceiling( length(gene.range)/options$n.cpu );
		sfExport("n.percpu", "gene.range", "Y.delt", "X", "par_null", "gen.lst", "options", "PF.pars" );
		
		res.cluster <- sfClusterApplyLB( 1:options$n.cpu, cpu.fun);
		sfStop();

		cat("Stopping parallel computing......\n"); 
		for(i in 1:length(res.cluster))
			gen.ret <- rbind(gen.ret, res.cluster[[i]]);
	}
	else
	{
		cat("Starting the LSKAT estimate for each gene......\n");
		gen.ret <- longskat_gene_task(gene.range, PF, Y.delt, X, par_null, gen.lst, options$w.common, options$w.rare, options$run.cpp, options$debug );
	}
	
	tm <- proc.time() - tm.start;
	cat( "* RUNTIME =", tm[3], "seconds \n");

	ret = list(gen.ret=gen.ret, pars=PF.pars);
	class(ret) <- "LSKAT.gen.ret";

	return(ret);
}

summary.LSKAT.gen.ret<-function(r.gen)
{
	df.gene <- r.gen$gen.ret[, c(2,3,4,5,6,8),drop=F];
	df.gene <- df.gene[ order(df.gene[,6]),,drop=F ]; 
	df.gene[,6] <- df.gene[,6]*NROW(df.gene);
	
	n.sig <- c();
	n.sig.i <- c();
	for( i in 20:4 )
	{
		n <- length( which(df.gene[,6]<10^(-i)) );
		if(n>0)
		{
			n.sig   <- c(n.sig, n);
			n.sig.i <- c(n.sig.i, i);
		}	
	}

	cat("  Summary of LSKAT result:\n");	
	cat("* PHE.LOG.FILE = ", r.gen$pars$file.phe.long, "\n")	
	cat("* PHE.COV.FILE = ", r.gen$pars$file.phe.cov, "\n")	
	if (length(n.sig.i)>0)
	{
		cat("* Significant Genes (Bonferroni Correct):\n")	
		for(i in 1:length(n.sig.i))
			cat("* <= 10^(-", n.sig.i[i], ") Levels:", n.sig[i], "\n")	
	}
	
	cat("* Top 20 Genes (Bonferroni Correct):\n")	
	cat("* [", "X", "]", "Gene Name", "\t", "Chromosome", "\t", "Position", "\t", "SNP Counts", "\t", "p-Value","\n");	
	for(i in 1:20)
	{
		if( i > NROW(df.gene) ) next;
		cat("* [", i, "]", as.character(df.gene[i,1]), "\t", df.gene[i,2], "\t", df.gene[i,3], "\t", df.gene[i,4], "\t", df.gene[i,6],"\n");	
	
	}
}

plot.LSKAT.gen.ret<-function( r.gen, pdf.file=NA, title="",  y.max=NA )
{
	if(is.na(pdf.file)) pdf.file<-paste(r.gen$pars$file.phe.long, ".pdf", sep="");

	#CHR, POS, NAME, PV
	df.gene <- r.gen$gen.ret[,c(3,4,2,8)];
	df.gene <- df.gene[ order(df.gene[,1], df.gene[,2]), ]; 
	df.gene[,4] <- df.gene[,4]*NROW(df.gene);

	pdf( pdf.file, width=6, height=3.5)

	par(mar=c(3.2, 3.5,  2.5, 1))
	par(mgp=c(1.4, 0.3, 0), tck=-0.03)
	draw_manhattan( df.gene[,c(1,3,4)], map.title=title, 0.0001/NROW(df.gene), 0.7, y.max= y.max );

	dev.off();

	cat("* Manhattan plot is save into ", pdf.file, "\n");	
}
