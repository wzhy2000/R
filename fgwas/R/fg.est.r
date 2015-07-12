library(mvtnorm);

fg_dat_est<-function( pheY, pheX, pheT, no.curve=NULL, no.covar=NULL )
{
	if( is.null(no.curve) )
		no.curve <- fg_fit_curve( pheY, pheX, pheT );

	if( is.null(no.covar) )
		no.covar <- fg_fit_covar( pheY, pheX, pheT );

	f.curve <- fg_get_curve( no.curve );
	r <- proc_est_curve( pheY, pheX, pheT, f.curve )
	if( r$error )
		return(list(error=T, err.info="Can not estimate the parameter of mean vector according to the curve function" ) )

	parX.est       <- r$par[1:(1+NCOL(pheX))];
	par.curve      <- r$par[-(1:(1+NCOL(pheX)))];
		
	range <- proc_est_curve_range(pheY, pheX, pheT, f.curve, par.init=r$par);
		
	parX.lower     <- range$lower[1:(1+NCOL(pheX))];
	parX.upper     <- range$upper[1:(1+NCOL(pheX))];

	f.curve$est_par <- par.curve;
	f.curve$lower   <- range$lower[-(1:(1+NCOL(pheX)))];
	f.curve$upper   <- range$upper[-(1:(1+NCOL(pheX)))];

	f.covar <- fg_get_covar( no.covar );
	r.est <- fg_est_covar( r$y.delt, NULL, pheT, f.curve, f.covar ); 
	if ( r.est$error ) 
		return(list(error=T, err.info="Can not estimate the parameter of mean vector according to the curve function" ) );

	f.covar$est_par <- r.est$par;
	
	return(list(error=F, 
		curve = f.curve, 
		covar = f.covar, 
		parX  = list( est_par = parX.est, lower = parX.lower, upper = parX.upper)));
}

fn_get_delt<-function( pheY, pheX, pheT, f.curve, parin )
{
	par_X <- c( parin[1] );
	if ( !is.null(pheX) )
		par_X <- c( par_X, parin[ 2:(NCOL(pheX)+1)]);

	par_c <- parin[ -(1:(NCOL(pheX)+1)) ];

	mu_gen <- f.curve$func( par_c, pheT )
	if(any(is.na( mu_gen ))) 
		return(NaN);

	if( is.vector( mu_gen ) )
		y_delt <- t(t(pheY) - mu_gen) - rowSums( cbind(1, pheX) * par_X )
	else
		y_delt <- pheY - mu_gen  - rowSums( cbind(1, pheX) * par_X );

	return( y_delt );	
}
	
proc_est_curve<-function(  pheY, pheX, pheT, f.curve, par.init=NULL, n.loop=10 )
{
	get_init_curve_par<-function( pheY, pheX, pheT, f.curve )
	{
		if ( !is.null( f.curve$fn_est ) )
			return( f.curve$fn_est(pheY, pheX, pheT) );
		
		r.max <- max( pheY, na.rm=T );
		r.min <- min( pheY, na.rm=T );

		return(runif( f.curve$pnum, r.min, r.max) );
	}

	get_rand_curve_par<-function( pheY, pheX, pheT, f.curve, parin )
	{
		uc <- check_fittness( pheY, pheX, pheT, f.curve, parin )

		if(uc$fit)
		{
			return( parin*runif(f.curve$pnum, 0.9, 1.1));
		}
		else
		{
			if ( uc$over0.05 < NCOL(pheY)/4 )	
				return( parin*runif(f.curve$pnum, 0.5, 1.5))
			else
			{
				par <- get_curve_init( pheY, pheX, pheT, f.curve );
				return( par * runif( f.curve$pnum, 0.5, 1.5 ) );
			}	
		}
	}	

	
	# parin : par.X[1..n], par.curve.
	fn_mle_est<-function( parin, extra_par )
	{
		pheY <- extra_par$pheY; 
		pheX <- extra_par$pheX;
		pheT <- extra_par$pheT; 
		fn_curve  <- extra_par$fn_curve; 
		
		y_delt <- fn_get_delt( pheY, pheX, pheT, f.curve, parin );
		
		A <- sum( y_delt^2, na.rm=T );

		return(A);
	}	
	

	if( is.null( par.init) )
	{
		par.init <-  c( mean(pheY, na.rm=T) );
		if( !is.null(pheX)) par.init <- c( par.init, rep(mean(pheY, na.rm=T), NCOL(pheX) ) );
		par.init <- c( par.init, get_init_curve_par(  pheY, pheX, pheT, f.curve ) );
	}

	h0 <- proc_mle_loop(  pheY, pheX, pheT, 
			f.curve, 
			fn_mle_est, 
			mle_extra_par = list( pheY=pheY, pheX=pheX, pheT=pheT, fn_curve=f.curve$func ),
			parin = par.init, 
			fn_init_par = get_init_curve_par, 
			fn_rand_par = get_rand_curve_par, 
			n.loop = 10 )

	if(fg.sys$log>=LOG_INFO)  
		cat( "MU[F]", h0$value, h0$par, "\n");

	r.check <- check_fittness( pheY, pheX, pheT, f.curve, h0$par );
	if( r.check$fit && is.finite( h0$value ) )
	{
		y.delt <- fn_get_delt( pheY, pheX, pheT, f.curve, h0$par );
		return(list(error=F, par=h0$par, value=h0$val, y.delt=y.delt))
	}
	else
		return(list(error=T, par=NA, val=NA));
		
}

proc_est_curve_range<-function( pheY, pheX, pheT, f.curve, par.init )
{
	n.obs <- NROW( pheY );
	mu.pars <- c();
	
	loop <- 0;
	while( loop < fg.sys$mu.range.loop )
	{
		y.sub <- sample(n.obs)[1:round( n.obs * runif(1,0.5,0.9) )]
		
		pheY0  <- pheY[y.sub,,drop=F];
		if( !is.null(pheX) ) pheX0 <- pheX[y.sub,,drop=F];
		if( !is.null(dim(pheT)) )  pheT0 <- pheT[y.sub,,drop=F] else pheT0 <- pheT;
		
		r <- proc_est_curve( pheY0, pheX0, pheT0, f.curve, par.init, n.loop = 2 );
		if (r$error) next;
					
		mu.pars <- rbind(mu.pars, r$par);
		loop <- loop + 1;
	}

	mu.lower <-c();	
	mu.upper <-c();	
	
	for(i in 1:NCOL(mu.pars) )
	{
		mu.lower <- c(mu.lower, min(mu.pars[,i]) )
		mu.upper <- c(mu.upper, max(mu.pars[,i]) )
	}

	return(list(lower=mu.lower, upper=mu.upper))
}

fg_est_covar<-function( Y.delt, pheX, pheT, f.curve, f.covar, par.init=NULL, n.loop=10 )
{
	library(mvtnorm);
	
	get_init_covar_par<-function( Y.delt, pheX, pheT, f.covar)
	{
		if ( !is.null(f.covar$fn_est) )
			return( f.covar$fn_est(Y.delt, pheX, pheT) );
		
		par <- c();
		for(i in 1:f.covar$pnum)
		{
			r.min <- ifelse( is.finite( f.covar$lower[i]), f.covar$lower[i], -1*sd(Y.delt, na.rm=T) );
			r.max <- ifelse( is.finite( f.covar$upper[i]), f.covar$upper[i], sd(Y.delt, na.rm=T) )

			par <- c(par, runif(1, r.min, r.max) );
		}

		return( par );
	}
	
	get_rand_covar_par<-function( Y.delt, pheX, pheT, f.covar, parin )
	{
		par <- parin* runif(f.covar$pnum, 0.9, 1.1);
		return(par);
	}

	#parin:
	# phi1, s1, phi2, s2, a1, b1, r12, a2, b2, r21
	fn_mle_est<-function( parin, extra_par)
	{
		y.delt  <- extra_par$Y.delt; 
		pheX    <- extra_par$pheX;
		pheT    <- extra_par$pheT; 
		f.covar <- extra_par$f.covar;
		
		cov.mat  <- f.covar$func( parin, pheT);
		if(any(is.na(cov.mat))) return(NaN);

		pv <- dmvnorm( y.delt, rep(0, NCOL(y.delt)), cov.mat, log=T);

		if(any(is.infinite(pv)))
			return(NaN);

		A <- sum( pv );
		return( -A );
	}

	if( is.null( par.init) )
		par.init <- get_init_covar_par(  Y.delt, pheX, pheT, f.covar );

	h0 <- proc_mle_loop( Y.delt, pheX, pheT, 
			f.covar, 
			fn_mle_est, 
			mle_extra_par=list(Y.delt=Y.delt, pheX=pheX, pheT=pheT, f.covar=f.covar),
			parin = par.init, 
			fn_init_par = get_init_covar_par, 
			fn_rand_par = get_rand_covar_par, 
			n.loop=10 )

	if( fg.sys$log >= LOG_INFO ) cat( "COV(F)", h0$value, h0$par, "\n");

	if( is.finite( h0$value ) )
		return(list(error=F, par=h0$par))
	else
		return(list(error=T, par=NA));
}

proc_mle_loop<-function(  pheY, pheX, pheT, f.obj, fn_mle, mle_extra_par, parin, fn_init_par, fn_rand_par, n.loop=10 )
{
	h0<-list( value=Inf, par=parin ); 
	
	while( is.infinite(h0$value) )
	{
		try( h0 <- optim( parin, fn_mle, extra_par=mle_extra_par, method = "BFGS" ) );
		if ( class(h0)=="try-error" || is.na(h0) || is.infinite(h0$value) )
		{
			parin <- fn_init_par( pheY, pheX, pheT, f.obj );
			reset_seed();
			next;
		}
	}

	if( fg.sys$log>=LOG_INFO ) cat( "X[0]", h0$value, h0$par, "\n");

	parin0 <- h0$par;
	loop <- 1;
	unfit <- 0;
	h2.better <-NULL;
	loop.optim <-1;
	
	min.fit<-Inf;
	
	while ( loop < n.loop )
	{
		h2 <- NA;
		parinx<- fn_rand_par( pheY, pheX, pheT, f.obj, parin0 );
		loop.optim <- loop.optim+1;
		try( h2 <- optim( parinx, fn_mle, extra_par=mle_extra_par,  method = ifelse(loop.optim%%2==1, "Nelder-Mead", "BFGS") ), TRUE );

		if (class(h2)=="try-error" || any(is.na(h2)) || h2$value<0  )
		{
			if(loop.optim>100) 
			{
				loop.optim <- 0;
				loop <- loop+1;

				if (is.infinite(h0$value) && n.loop - loop< 2 )
					n.loop <- n.loop + 1;
			}
			
			reset_seed();
			next;
		}
		else
		{
			uc <- check_fittness( pheY, pheX, pheT, f.obj, h2$par )
			if ( !uc$fit )
			{	
				if (uc$over0.05<min.fit)
				{
					parin0    <- h2$par;
					min.fit <- uc$over0.05;
				}

				reset_seed();
				next;
			}

			if(fg.sys$log>=LOG_DEBUG) cat( "MU(L", loop, ")", uc$over0.05, "/", h2$value, h2$par, "\n");

			if ( h2$value < h0$value && h2$value>0 ) h0 <- h2;
			if ( h2$value >0 && h0$value<0 ) h0 <- h2;

			parin0 <- h0$par;
			min.fit <- Inf;
			
			if (is.infinite(h0$value) && n.loop - loop< 2 )
				n.loop2 <- n.loop2 + 1;
		}

		loop.optim <- 0;
		loop <- loop+1;
	}

	if( fg.sys$log>=LOG_INFO )  cat( "X[F]", h0$value, h0$par, "\n");

	uc <- check_fittness( pheY, pheX, pheT, f.obj, h0$par );
	if( uc$fit && is.finite(h0$value) )
		return(list(error=F, par=h0$par, value=h0$value))
	else
		return(list(error=T, par=NA, value=NA));
}

check_fittness<-function( pheY, pheX, pheT, f.obj, parin )
{
	if( !is.curve(f.obj) ) return(list(fit=T));
		
	y_delt <- fn_get_delt( pheY, pheX, pheT, f.obj, parin );
	y_sd   <- apply( pheY, 2, function(x) sd(x) )

	py.prob  <- pnorm( colMeans(y_delt), mean=0, sd=y_sd, lower.tail=F );

	py.05 <- length( which(py.prob<0.05) );
	
	return(list(fit=ifelse(py.05>1, F, T), over0.05=py.05) );
}

fg_fit_curve<-function( pheY, pheX, pheT )
{
#todo;
}

fg_fit_covar<-function( pheY, pheX, pheT )
{
#todo;
}