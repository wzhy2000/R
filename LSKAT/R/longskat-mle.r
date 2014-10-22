library(mvtnorm)

longskat_est_model<-function( y.log, y.cov, bTime=0, g.maxiter=20, debug=F )
{
	nrow <- dim(y.log)[1];
	ncol <- dim(y.log)[2];
	nCov <- dim(y.cov)[2];

	get_par<-function(par, y, y.cov)
	{	
		par_rho<-par[1]
		if ( par_rho<0 || par_rho>=0.99 )
			return(NaN);

		sig_a<-par[2]
		sig_b<-par[3]
		sig_e<-par[4]
		par_u<-par[5]

		par_cov <- par[c(6:(5+nCov)) ];
		par_t  <- 0;
		par_t2 <- 0;
		if( bTime!=0 )
		{
			par_t  <- par[6+nCov];
			par_t2 <- par[7+nCov];
		}
		
		y.delt <- y - par_u;
		for(i in 1:nCov)
			y.delt <- y.delt - y.cov[,i]*par_cov[i];

		y.delt <- y.delt - par_t * (c(1:ncol)-1) - par_t2 * (c(1:ncol)-1)^2;		
		y.temp <- c(y.delt)
		if ( length(which(is.na(y.temp)))>0) y.temp[ is.na(y.temp) ]<-0;
		y.delt <- array(y.temp, dim=c(nrow, ncol));
		
		AR.1 <- array(0,dim=c(ncol,ncol));
		for(i in 1:ncol)
		for(j in 1:ncol)
			AR.1[i,j] <- par_rho^abs(i-j);
		
		sigma <- diag(sig_e^2, ncol) + sig_a^2 + sig_b^2*AR.1;

		A <- -sum( log( dmvnorm( y.delt,  rep(0, ncol), sigma ) ) )
		return( A );
	}

	par.init <- c( 0.5, sd(y.log, na.rm=T), sd(y.log, na.rm=T), sd(y.log, na.rm=T), mean(y.log, na.rm=T) );
	for(i in 1:nCov)
		par.init <- c(par.init, 0.1)
		
	if (bTime!=0)
		par.init <- c( par.init, 0.1, 0.1 )
	
	#par.init <- c(rho=0.75, sig_a=0.2, sig_b=0.3, sig_e=0.1, u=1, a=0.5, b=0.5); 
	loop.n <- 0;
	min.val <- Inf;
	min.par <- par.init;
	
	while( loop.n <= g.maxiter )
	{
		r0 <- try( optim( par.init, get_par, y = y.log, y.cov=y.cov, method = "BFGS", control=list(maxit=500) ), silent = TRUE );
		if (class(r0)=="try-error")
		{
			if (min.val >= 1e8)
				loop.n <- loop.n + 0.2;
			par.init <- min.par*runif( length(min.par) );
			next;
		}

		loop.n <- loop.n+1;
		if (r0$val < min.val)
		{
			if(debug) cat("  LOOP =", loop.n, "/", g.maxiter, " val=", r0$val,  "par=", r0$par, "\n");
			min.val <- r0$value;
			min.par <- r0$par;
			par.init <- min.par;
		}

		par.init <- par.init*runif( length(par.init), 0.8, 1.2 );

	}
	
	if(debug) cat("  LOOP =", loop.n, " min.val=", min.val, "par=", min.par, "\n");
	
	if( is.infinite( min.val) || is.na(min.val) || any(is.na(min.par))  )
		return( list(bSuccess=F) );

	par_cov <- min.par[c(6:(5+nCov))];
	par_t  <- 0;
	par_t2 <- 0;
	if( bTime!=0 )
	{
		par_t  <- min.par[6+nCov];
		par_t2 <- min.par[7+nCov];
	}
			
	y.delt <- y.log - min.par[5];
	for(i in 1:nCov)
		y.delt <- y.delt - y.cov[,i]*par_cov[i];
	y.delt <- y.delt - par_t * (c(1:ncol)-1) - par_t2 * (c(1:ncol)-1)^2;	
	
	r <- list(bSuccess=T, u=min.par[5], sig_a=abs(min.par[2]), sig_b=abs(min.par[3]), sig_e=abs(min.par[4]), 
		  cov=min.par[6:(5+nCov)], rho=min.par[1], val=min.val, y.delt=y.delt );

	return(r);

}
