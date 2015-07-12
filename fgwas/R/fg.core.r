#------------------------------------------------------------------------------------------
#
# Functional Mapping Method
#
# !!!! fy[which(fy<=.Machine$double.eps)] <- .Machine$double.eps;
#
#------------------------------------------------------------------------------------------

library(mvtnorm);

proc_h0_mu<-function( parin, f.curve, d.phe, d.times )
{
	gen_mu  <- f.curve$func( parin, d.times )
	mu.delt <- t(t(d.phe) - gen_mu);

	A <- sum(mu.delt^2);
	return(A);
}

proc_h0_cov<-function( parin, f.covar, gen_mu, d.phe, d.times )
{
	cov.mat <- f.covar$func( parin, d.times )
	if (any(is.na(cov.mat)))
		return (NaN);

	pv <- NA;
	try(pv <- dmvnorm( d.phe, gen_mu, cov.mat, log=T), fg.sys$try.slient );
	if ( class(pv)=="try-error" || is.na(pv) ) 
		return (NaN);

	A <- sum(pv);
	return(-A);
}

proc_h0_cov_bfgs<-function( parin, f.covar, gen_mu, d.phe, d.times )
{
	cov.mat <- f.covar$func( parin, d.times )
	if (any(is.na(cov.mat)))
		return (NaN);

	pv <- dmvnorm( d.phe, gen_mu, cov.mat, log=T);
	A <- sum(pv);

	return(-A);
}

proc_mle_mu <- function( f.curve, d.phe, d.umiss, d.times, n.loop=5)
{
	mu <- proc_cache_find( d.umiss );
	if ( !is.null( mu) ) 
	{
		cat( "MU(C0)", mu, "\n");
		return(mu)
	}
	
	loop <- 0;
	h0.better <- NULL;
	h2 <- NULL;
	unfit <- 0;

	while(loop<n.loop)
	{
		h0 <- NA;
		parin <- f.curve$par * runif( f.curve$pnum, 0.95, 1.05 );
		try( h0 <- optim( parin, proc_h0_mu, f.curve=f.curve, d.phe=d.phe, d.times=d.times, method = "L-BFGS-B", lower=f.curve$lower, upper=f.curve$upper ),  fg.sys$try.slient );
		if ( class(h0)=="try-error" || is.na(h0) || h0$value<0 )
			next;

		uc.ret <- compare_unfit( unfit, d.phe, h0, h0.better, 1:f.curve$pnum, f.curve$func, d.times )
		if ( !uc.ret$fit )
		{
			h0.better <- uc.ret$h2.better;
			unfit <- unfit+1;

			reset_seed();
			next;
		}

if( fg.sys$log >= LOG_DEBUG ) cat( "MU(L", loop, ")", h0$value, h0$par, "\n");

		if (is.null(h2)) h2 <- h0; ;
		if (h2$value > uc.ret$h2.better$value ) h2 <- uc.ret$h2.better;

		h0.better <- NULL;
		loop <- loop+1;
	}

	proc_cache_save( d.umiss, h2$par, h2$value );
	
	return(h2$par);
}


proc_mle_h0 <- function( f.curve, f.covar, f.parX, pheY, pheX, pheT, debug, n.loop=5, b.permu=F)
{
	loop <- 0;
	h0   <-NA;
	loop.optim <- 0;
	
	parin <- c(f.parX$est_par, f.curve$est_par, f.covar$est_par);
	while ( loop < n.loop )
	{
		h2 <- NA;
		parinx<- parin*runif( length(parin), 0.8, 1.2);
		
		loop.optim <- loop.optim + 1;
		if ( loop.optim%%2 == 1 )
			try( h2 <- optim( parinx, proc_h0_cov, f.covar=f.covar, gen_mu=gen_mu, d.phe=d.phe, d.times=d.times, method = "Nelder-Mead" ), fg.sys$try.slient )
		else
			try( h2 <- optim( parinx, proc_h0_cov_bfgs, f.covar=f.covar, gen_mu=gen_mu, d.phe=d.phe, d.times=d.times, method = "L-BFGS-B", lower=f.covar$lower, upper=f.covar$upper ), fg.sys$try.slient );

		if (class(h2)=="try-error" || any(is.na(h2)) || h2$value<0 )
		{
			reset_seed();
			next;
		}
		
		if( h2$convergence!=0 )
			cat(" Optim Failed ", h2$convergence, "\n")
		else 
		{
if( fg.sys$log >= LOG_DEBUG && !b.permu ) cat( "H0(L", loop, ")", h2$value, h2$par, "\n");
			if (any(is.na(h0))) h0 <-h2;
			if (h2$value < h0$value ) h0 <- h2;
			loop <- loop+1;
		}
	}


if( fg.sys$log >= LOG_INFO && !b.permu ) cat( "H0(F", loop, ")", h0$value, h0$par, "\n");

	h0$par <- c(h0$par, parin.mu);
	
	return(h0);
}


proc_top_select<-function( mat.scan, n.top)
{
	top.snp <- NA;
	sort.lr2 <- sort.int( mat.scan[,4], decreasing=T, index.return=T );
	if(n.top<=dim(mat.scan)[1])
		top.snp <- mat.scan[sort.lr2$ix[1:n.top], ];
		
	return(top.snp);
}

proc_sig_select<-function(mat.scan, fg.perm)
{
	re0 <- mat.scan[,c(1,2,3,4)]
	
	sig.05 <- NA;
	if( is.finite(fg.perm$pcut.05) )
		sig.05 <- mat.scan[which(re0[,4]>=fg.perm$pcut.05), ];
		
	sig.01 <- NA;
	if( is.finite(fg.perm$pcut.01) )
		sig.01 <- mat.scan[which(re0[,4]>=fg.perm$pcut.01), ];
	
	return(list(sig.05=sig.05, sig.01=sig.01));
}

compare_unfit<-function( unfit, d.phe, h2, h2.better, parset, f.curveunc, d.times )
{
	phe.mat <- as.matrix(d.phe);
    	py.mean <- colMeans(phe.mat);
	py.sd   <- apply(phe.mat, 2, function(x) sd(x) )

	mu <- f.curveunc ( h2$par[parset], d.times );
	mu.r0 <- range(abs(phe.mat - mu))
	py.prob  <- pnorm(abs( mu - py.mean), mean=0, sd=py.sd, lower.tail=F);

	py.05 <- length(which(py.prob<0.05));
	
	fit <- FALSE;
	if (py.05 <= 6)
		fit <- TRUE
	else if (unfit>5 && py.05 <= 9)
		fit <- TRUE
	else if (unfit>10 && py.05 <= 12)
		fit <- TRUE
	else if (unfit>100)
		fit <- TRUE

if(fg.sys$log>=LOG_DEBUG) cat("UNFIT0:", py.05, "/", unfit, "\n");
		
	if (is.null(h2.better))
	{
		h2.better  <- h2
		h2.better$py05 <- py.05;
	}

	if (h2.better$py05>py.05)
	{
		h2.better  <- h2
		h2.better$py05 <- py.05;
	}
	
	return(list(fit=fit, h2.better = h2.better));
}

proc_mle<-function( snp.idx, snp, pheY, pheX, pheT, phe.est, debug, n.loop=5, b.permu=F )
{
	d.scaf<- snp[1]; 
	d.loci<- snp[2];
	d.gen <- snp[-c(1,2)];
	gen.par <- max(d.gen0) - min(d.gen0) +1;

	h0 <- proc_mle_h0( phe.est$curve, phe.est$covar, phe.est$parX, pheY, pheX, pheT, debug, n.loop, b.permu );
	h1 <- proc_mle_h1( phe.est$curve, phe.est$covar, phe.est$parX, pheY, pheX, pheT, snp, debug, n.loop, b.permu );

	r.val <- ( h0$value - h1$value )*2
	r.pv  <- pchisq(r.val, df=(gen.par-1)*fg.dat$f.curve$pnum, lower.tail=FALSE);

	re <- c( unlist(d.scaf), unlist(d.loci), gen.par, r.val, r.pv, h0$par, h1$par );

	if( fg.sys$log >= LOG_INFO) 
		if (!b.permu)
			cat("snp[", snp.idx, "]" , dim(d$phe0)[1], "loci=", re[1], "/", re[2], re[3], "LR2=", re[4], re[5], "Parm=", re[6:12], "\n");
		
	return(re); 
}

proc_optim_loop <- function( parin, fn_mle, extra_mel_par, debug, n.loop=1, b.permu=F)
{
	loop       <- 0;
	h.best     <- NULL;
	loop.optim <- 0;
	
	parinx     <- parin;
	while ( loop < n.loop )
	{
		loop.optim <- loop.optim + 1;
		if ( loop.optim%%2 == 1 )
			try( h2 <- optim( parinx, fn_mle, extra_par=extra_mel_par, method = "Nelder-Mead" ), fg.sys$try.slient )
		else
			try( h2 <- optim( parinx, fn_mle, extra_par=extra_mel_par, method = "BFGS"), fg.sys$try.slient );

		if (class(h2)=="try-error" || any(is.na(h2)) || h2$value<0 )
		{
			parinx<- parin * runif( length(parin), 0.8, 1.2);
			reset_seed();
			next;
		}
		
		if( h2$convergence!=0 )
			cat(" Optim Failed ", h2$convergence, "\n")
		else 
		{
if( fg.sys$log >= LOG_DEBUG && !b.permu ) cat( "H0(L", loop, ")", h2$value, h2$par, "\n");

			if ( is.null(h.best) ) h.best <-h2;
			if ( h2$value < h.best$value ) h.best <- h2;
			loop <- loop+1;
		}

		parinx<- parin * runif( length(parin), 0.9, 1.1);
	}


if( fg.sys$log >= LOG_INFO && !b.permu ) cat( "H0(F", loop, ")", h.best$value, h.best$par, "\n");

	return(h.best);
}

proc_dmvnorm<-function(y_delt, mu_gen, cov_mat)
{
	pv <- c();
	if(is.vector( mu_gen ))
		try( pv <- dmvnorm( y_delt, rep(0, NCOL(mu_gen)), cov_mat, log=T), fg.sys$try.slient )
	else
	{
		for(i in 1:NROW(y.delt))
			pv <- c( pv, dmvnorm( y_delt[i,,drop=F], rep(0, NCOL(mu_gen)), cov_mat[[i]], log=T) );
	}
	
	return(pv);
}

proc_mle_h0 <-function ( f.curve, f.covar, f.parX, pheY, pheX, pheT, debug, n.loop=1, b.permu=F )
{
	fn_mle_h0 <- function( parin, extra_par)
	{
		f.covar  = extra_par$f.covar; 
		f.curve  = extra_par$f.curve; 
		f.parX   = extra_par$f.parX; 
		pheY     = extra_par$pheY; 
		pheX     = extra_par$pheX; 
		pheT     = extra_par$pheT;
		
		par_X    = parin[1:(NCOL(pheX)+1)];
		parin    = parin[-c(1:(NCOL(pheX)+1))];
		par_curve= parin[1:f.curve$pnum];
		parin    = parin[-c(1:f.curve$pnum)];
		par_covar= parin;
		
		cov_mat <- f.covar$func( par_covar, pheT )
		if (any(is.na(cov_mat)))
			return (NaN);

		mu_gen <- f.curve$func( par_curve, pheT )
		if( is.vector( mu_gen ) )
			y_delt <- t(t(pheY) - mu_gen) - rowSums( cbind(1, pheX) * par_X )
		else
			y_delt <- pheY - mu_gen  - rowSums( cbind(1, pheX) * par_X );

		pv <- proc_dmvnorm (y_delt, mu_gen, cov_mat);
		if ( class(pv) == "try-error" || any( is.na( pv ) ) )
			return (NaN);

		A <- sum(pv);
		return( -A );
	}
		
	parin <- c(f.parX$est_par, f.curve$est_par, f.covar$est_par);
	extra_par <- list(f.curve= f.curve, f.covar=f.covar, f.parX=f.parX, pheY=pheY, pheX=pheX, pheT=pheT);

	h0 <- proc_optim_loop( parin, fn_mle_h0, extra_par, debug, n.loop, b.permu);
	return(h0);
}

proc_mle_h1<-function( f.curve, f.covar, f.parX, pheY, pheX, pheT, snp, debug, n.loop=1, b.permu=F )
{
	par.curve <- c();

	if (length(which(snp==0))>0)
	{
		id.sel <- which(snp==0);
		par.mu <- proc_mle_mu( f.curve, f.parX, pheY[id.sel,,drop=F], pheX[id.sel,,drop=F], pheT[id.sel,,drop=F], n.loop );
		par.curve <- rbind( par.curve, par.mu); 		
	}
	else
		par.curve <- rbind( par.curve, rep(NA, f.curve$pnum) ); 

	if (length(which(snp==1))>0)
	{
		id.sel <- which(snp==1);
		par.mu <- proc_mle_mu( f.curve, f.parX, pheY[id.sel,,drop=F], pheX[id.sel,,drop=F], pheT[id.sel,,drop=F], n.loop );
		par.curve <- rbind( par.curve, par.mu); 		
	}
	else
		par.curve <- rbind( par.curve, rep(NA, f.curve$pnum) ); 

	if (length(which(snp==2))>0)
	{
		id.sel <- which(snp==2);
		par.mu <- proc_mle_mu( f.curve, f.parX, pheY[id.sel,,drop=F], pheX[id.sel,,drop=F], pheT[id.sel,,drop=F], n.loop );
		par.curve <- rbind( par.curve, par.mu); 		
	}
	else
		par.curve <- rbind( par.curve, rep(NA, f.curve$pnum) ); 
	
	
	gen_mu <- t(array(0, dim=dim(d.phe)));
	qq.idx <- which(d.gen==0)
	if (length(qq.idx)>0)  gen_mu[,qq.idx]<-f.curve$func( par.fin[1,], d.times  );
	Qq.idx <- which(d.gen==1)
	if (length(Qq.idx)>0)  gen_mu[,Qq.idx]<-f.curve$func( par.fin[2,], d.times  );
	QQ.idx <- which(d.gen==2)
	if (length(QQ.idx)>0)  gen_mu[,QQ.idx]<-f.curve$func( par.fin[3,], d.times  );
	
	phe.delt <- d.phe - t(gen_mu);
	
	gen_mu <-  rep(0, length(d.times) )
	loop.optim <- 0;
	loop <-0;
	h0<-NA;
	
	if( fg.sys$log >= LOG_DEBUG  && !b.permu ) show(par.fin);
	
	while ( loop < n.loop )
	{
		h2 <- NA;
		parinx<- f.covar$par * runif( f.covar$pnum, 0.8, 1.2 );

		loop.optim <- loop.optim + 1;
		if (loop.optim%%2==1)
			try( h2 <- optim( parinx, proc_h0_cov, f.covar=f.covar, gen_mu=gen_mu, d.phe=phe.delt, d.times=d.times, method = "Nelder-Mead" ), fg.sys$try.slient )
		else
			try( h2 <- optim( parinx, proc_h0_cov_bfgs, f.covar=f.covar, gen_mu=gen_mu, d.phe=phe.delt, d.times=d.times, method = "L-BFGS-B", lower=f.covar$lower, upper=f.covar$upper ), fg.sys$try.slient );
		
		if (class(h2)=="try-error" || any(is.na(h2)) || h2$value<0 )
		{
			reset_seed();
			next;
		}
		else
		{

		if( fg.sys$log >= LOG_DEBUG  && !b.permu ) cat( "H1(L", loop, ")", h2$value, h2$par, "\n");

			if (any(is.na(h0))) h0 <-h2;
			if (h2$value < h0$value ) h0 <- h2;
		}

		loop <- loop+1;
	}

	h0$par <- c(h0$par, par.fin[1,], par.fin[2,], par.fin[3,]);

if( fg.sys$log >= LOG_INFO  && !b.permu ) cat( "H1(F", loop, ")", h0$value, h0$par, "\n");

	return(h0);

}
