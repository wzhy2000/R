#--------------------------------------------------------------
# curve.curve_mle
#
#  Used by permutaion.
#
#	prob: the probility of Qq
#     y : the phenotype * times data
#    par: initial parameter for likelihood
#         including a1,b1,r1,a0,b0,r0,rho,sigma2
#--------------------------------------------------------------
curve.mlefunc<-function( par, y, time.std )
{
	len.cov <- FM2.covar$get_par_num(y);
	par.covar <- par[1:len.cov];
	if (!FM2.covar$is_valid(par.covar))
		return(NaN);
	
	sig.inv <- FM2.covar$get_inv_mat(par.covar, time.std, FM2.curve$trait_num);
	sig.det <- FM2.covar$get_mat_det(par.covar, time.std, FM2.curve$trait_num);

	m  <- length(y[1,]);
	n  <- length(y[,1]);

	gen_par<- par[(len.cov+1):(len.cov+ FM2.curve$par_num )];
	mu0 <- FM2.curve$get_mu( gen_par, time.std, FM2.curve$trait_num );
	yy0 <- y - matrix( rep(1,n),byrow=n )%*%matrix( mu0, nrow=1 ); 

	for ( i in 1:length(y[1,]) )
	{
		y.miss <- which( is.na(y[,i]) );
		if(!is.null(yy0)) yy0[y.miss, i]<-0;
	}
	
	fy0 <- c();
	for(i in 1:n)
		fy0<- c( fy0, (2*pi)^(-m/2)*(sig.det)^(-1/2)*exp(t(yy0[i,])%*%sig.inv%*%yy0[i,]/2*(-1) ) );
	A <- -sum(log(fy0));

	return (A);
}

curve.get_est_param<-function(dat)
{
	phenos <- as.matrix( dat$phenos_table );
	par.num <- FM2.curve$par_num;

	par <- c();
	lr2 <- NA;
	loop <-0;

	while( loop <= 5 )
	{
		parin <- FM2.covar$get_init_rand(dat);
		if (length(par)>2)
		{
			for (i in 1:par.num)
				parin <- c(parin, runif(1)* par[i+2]  );
		}
		else
			for (i in 1:par.num)
				parin <- c(parin, runif(1)*mean(phenos, na.rm=TRUE)  );

		r0 <- NA;
		try( r0<- optim( parin, curve.mlefunc, y = phenos, time.std=dat$sample_times, method ="BFGS" ), TRUE);
		if (class(r0)=="try-error")
			next;
		
		if(any(is.na(r0)))
			next;
		
		if (is.na(lr2))
		{
			par <- r0$par;
			lr2 <- r0$value;
		}

		if (r0$value < lr2)
		{
			par <- r0$par;
			lr2 <- r0$value;
		}

		cat("TRY:", r0$value, r0$par, "\n");
		loop <- loop + 1;
	}
	
	cat("Estimated parameters:", "log(L)=", lr2, "PAR=", par, "\n");
	return(par);
}

class_curve<-list(
	name 	 	= "Base Curve",
	get_est_param 	= curve.get_est_param,
	get_mlefunc     = curve.mlefunc);

