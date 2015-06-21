library(mvtnorm)

# y.log format:  shareid, trait1, ..., traitN
# y.time format: shareid, time1, ..., timeN
# y.cov format:  shareid, cov1, ..., covM

#public
longskat_est_model<-function( y.log0, y.cov0, y.time0 = NULL, y.cov.time=0, g.maxiter=20, debug=F )
{
	nrow <- NROW(y.log0);
	ncol <- NCOL(y.log0)-1;
	nCov <- NCOL(y.cov0)-1;
	
	#check id matched!
	
	if(is.data.frame(y.log0)) y.log <- data.matrix(y.log0[,-1,drop=F])
	if(is.data.frame(y.cov0)) y.cov <- data.matrix(y.cov0[,-1,drop=F])
	if(is.data.frame(y.time0)) y.time<- data.matrix(y.time0[,-1,drop=F])

	if(is.matrix(y.log0))  y.log <- y.log0[,-1,drop=F];
	if(is.matrix(y.cov0))  y.cov <- y.cov0[,-1,drop=F];
	if(is.matrix(y.time0)) y.time<- y.time0[,-1,drop=F];
	
	if( is.null(y.time0))
		y.time <- t(matrix(rep(c(1:ncol),nrow), nrow=ncol)) ;

show(head(y.log));
show(head(y.cov));
show(head(y.time));

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
		if( y.cov.time >0 )
		{
			par_t  <-  par[5 + nCov + c(1:y.cov.time)];
			for(i in 1:y.cov.time)
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
	
	est_par_cov<-function()
	{
		par.init.cov <- c(mean(y.log, na.rm=T));
		for(i in 1:nCov)
			par.init.cov <- c( par.init.cov, 1/( mean(y.cov[,i], na.rm=T)^2 + 1) );

		return(par.init.cov);			
	}

	# par[1] = rho
	# par[2] = sig_a
	# par[3] = sig_b
	# par[4] = sig_e
	# par[5] = par_u
	# par[6, 5+nCov] = par_cov
	# par[6+nCov, (6:7+nCov) ] = par_cov_time
	# e.g.
	# par.init <- c(rho=0.75, sig_a=0.2, sig_b=0.3, sig_e=0.1, u=1, a=0.5, b=0.5); 

	par.init.cov <- est_par_cov();
	par.init <- c( 0.5, sd(y.log, na.rm=T), sd(y.log, na.rm=T), sd(y.log, na.rm=T), par.init.cov );

	if (y.cov.time>0)
	{
		par.init <- c( par.init, rep( 1/( mean(y.time, na.rm=T)^2 +1) , y.cov.time ) )
		if(is.null(y.time))
			y.time <- t( t(ifelse(is.na(y.log),NA, 1))*(rep(1:NROW(y.log))))
	}
	
	tolerance <- 1;
	loop.n <- 0;
	min.val <- Inf;
	min.par <- par.init;
	
	while( loop.n <= g.maxiter && tolerance > 1e-5 )
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
			
			tolerance <- max( c( min.val, par.init) - c(r0$value, r0$par) );
			
			min.val <- r0$value;
			min.par <- r0$par;
			par.init<- min.par;
		}

		par.init <- par.init*runif( length(par.init), 0.8, 1.2 );
	}
	
	if(debug) cat("  Final =", loop.n, " min.val=", min.val, "par=", min.par, "\n");
	
	if( is.infinite( min.val) || is.na(min.val) || any(is.na(min.par))  )
		return( list(bSuccess=F) );

	par_cov <- min.par[c(6:(5+nCov))];
	
	y.delt <- y.log - min.par[5];
	for(i in 1:nCov)
		y.delt <- y.delt - y.cov[,i]*par_cov[i];
		
	par_t <- c();	
	if( y.cov.time>0 )
	{
		par_t  <- min.par[ 5 + nCov +c(1:y.cov.time)];
		for( i in 1:length(par_t))
			y.delt <- y.delt - (y.time^i) * par_t[i];	
	}
	
	pars <- list( mu     = min.par[5], 
			  	  rho    = min.par[1], 
			  	  sig_a  = abs(min.par[2]), 
			  	  sig_b  = abs(min.par[3]), 
			  	  sig_e  = abs(min.par[4]), 
		  	      par_cov= par_cov, 
		  	      par_t  = par_t );
	r.model <- list(par = pars, likelihood = min.val, y.cov.time=y.cov.time, y.delt=y.delt, y.time = y.time, y.cov = y.cov );
	
	class(r.model) <- "LSKAT.null.model";

	return(r.model);
}

#public
print.LSKAT.null.model <- function(r.model, useS4=FALSE)
{
	cat("  LSKAT/Gene Summary\n");	
	cat("[1] MLE Results:\n");	
	cat("* SIGMA_A =", 		   r.model$par$sig_a, "\n");
	cat("* SIGMA_B =",         r.model$par$sig_b, "\n");
	cat("* SIGMA_E =",         r.model$par$sig_e, "\n");
	cat("* RHO =",             r.model$par$rho, "\n");
	cat("* MU =",              r.model$par$mu, "\n");
	cat("* Beta(Cov)=",        r.model$par$par_cov, "\n");
	cat("* Beta(Time)=",       r.model$par$par_t, "\n");
	cat("* L(min) =",          r.model$likelihood, "\n");
}


#private
get_ylog_list<-function(y.long)
{
	y.long.list <-list(); 
	y.una <-apply(y.long, 1, function(y.i){ length( which(!is.na(y.i))); } )

	for(i in 1:NCOL(y.long))
		y.long.list[[i]] <- y.long[which(y.una==i),c(1:i),drop=F]
	
	return(y.long.list)
}
