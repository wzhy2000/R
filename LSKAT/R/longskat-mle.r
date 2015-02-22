library(mvtnorm)

longskat_est_model<-function( y.log, y.time = NULL, y.cov, nTime=0, g.maxiter=20, debug=F )
{
	nrow <- dim(y.log)[1];
	ncol <- dim(y.log)[2];
	nCov <- dim(y.cov)[2];

	get_par<-function(par, y, y.time, y.cov)
	{	
		par_rho<-par[1]
		if ( par_rho<0 || par_rho>=0.99 )
			return(NaN);

		sig_a<-par[2]
		sig_b<-par[3]
		sig_e<-par[4]
		par_u<-par[5]
		par_cov <- par[ c(6:(5+nCov)) ];

		AR.1 <- array(0,dim=c(ncol,ncol));
		for(i in 1:ncol)
		for(j in 1:ncol)
			AR.1[i,j] <- par_rho^abs(i-j);
		
		sigma <- diag(sig_e^2, ncol) + sig_a^2 + sig_b^2*AR.1;

		y.delt <- y - par_u;
		for(i in 1:nCov)
			y.delt <- y.delt - y.cov[,i]*par_cov[i];
			
		par_t   <- c();
		if(nTime >0 )
		{
			par_t  <-  par[5 + nCov + c(1:nTime)];
			for(i in 1:nTime)
				y.delt <- y.delt - (y.time^i) * par_t[i];
		}
		
		# y.temp <- c(y.delt)
		# if ( length(which(is.na(y.temp)))>0) y.temp[ is.na(y.temp) ]<-0;
		# y.delt <- array(y.temp, dim=c(nrow, ncol));
		# A <- -sum( dmvnorm( y.delt,  rep(0, ncol), sigma, log=T ) ) 
		
		#y.list <- get_ylog_list(y.delt);
		#A <- 0;
		#for(i in 1:length(y.list) )
		#{
		#	if(!is.null(y.list[[i]]) && NROW(y.list[[i]])>0)
		#		A <- A - sum( dmvnorm( y.list[[i]], rep(0, i), sigma[1:i, 1:i, drop=F], log=T ) );
		#}


		#y.list <- get_ylog_list(y.delt);
		A <- 0;
		for(i in 1:NROW(y.delt) )
		{
			t.sel <- which( !is.na(c(y.time[i,]) ) );
			sig <- sigma[t.sel, t.sel, drop=F];
			A <- A - sum( dmvnorm( y.delt[i,t.sel,drop=F], rep(0, length(t.sel)), sig, log=T ) );
		}

			
		return( A );
	}

	#par.init <- c(rho=0.75, sig_a=0.2, sig_b=0.3, sig_e=0.1, u=1, a=0.5, b=0.5); 
	par.init <- c( 0.5, sd(y.log, na.rm=T), sd(y.log, na.rm=T), sd(y.log, na.rm=T), mean(y.log, na.rm=T) );
	for(i in 1:nCov)
		par.init <- c(par.init, 1/( mean(y.cov[,i], na.rm=T)^2 + 1) )
	
	if (nTime>0)
	{
		par.init <- c( par.init, rep( 1/( mean(y.time, na.rm=T)^2 +1) , nTime ) )
		if(is.null(y.time))
			y.time <- t( t(ifelse(is.na(y.log),NA, 1))*(rep(1:NROW(y.log))))
	}
	
	loop.n <- 0;
	min.val <- Inf;
	min.par <- par.init;
	
	while( loop.n <= g.maxiter )
	{
		r0 <- try( optim( par.init, get_par, y = y.log, y.time=y.time, y.cov=y.cov, method = "BFGS", control=list(maxit=500) ), silent = F );
		if (class(r0)=="try-error")
		{
			if (min.val >= 1e8)
				loop.n <- loop.n + 0.2;
			par.init <- min.par*runif( length(min.par) );
			next;
		}

		loop.n <- loop.n+1;
		if ( r0$convergence==0 && r0$val<min.val)
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
	
	y.delt <- y.log - min.par[5];
	for(i in 1:nCov)
		y.delt <- y.delt - y.cov[,i]*par_cov[i];
		
	par_t <- c();	
	if( nTime>0 )
	{
		par_t  <- min.par[ 5 + nCov +c(1:nTime)];
		for( i in 1:length(par_t))
			y.delt <- y.delt - (y.time^i) * par_t[i];	
	}
	
	r <- list(bSuccess=T, u=min.par[5], rho=min.par[1], sig_a=abs(min.par[2]), sig_b=abs(min.par[3]), sig_e=abs(min.par[4]), 
		  par_cov=par_cov, par_t=par_t, val=min.val, y.delt=y.delt, y.t = y.time );

	return(r);

}

get_ylog_list<-function(y.long)
{
	y.long.list <-list(); 
	y.una <-apply(y.long, 1, function(y.i){ length( which(!is.na(y.i))); } )

	for(i in 1:NCOL(y.long))
		y.long.list[[i]] <- y.long[which(y.una==i),c(1:i),drop=F]
	
	return(y.long.list)
}
