#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
#  Non-stationary structured antedependence model (2)
#  
#  Variance : SAD1(2)
#  
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

SAD2.get_inv_mat <- function(par0, times, traits=1, options=list())
{
	par <- par0;
	if (class(par)!="list")
	{
		par<-list(rho=par0[1], 
			phi1  = par0[2],
			phi2  = par0[3],
			psi1  = par0[4],
			psi2  = par0[5],
			gama1 = par0[6],
			gama2 = par0[7] )
	}


	T <- length(times)*traits;

	sigepsi<-matrix(0,T,T);
  	L<-matrix(0,T,T)
 	diag(L) <- 1

 	#first get L
 	for (i in 1:T)
 	{
  	 	for (j in 1:T)
  	 	{
		 	if (((i-j)==1) & ((i %% 2)==0))
		 	{
			 	L[i,j]<-0
		 	}
		 	if (((i-j)==1) & ((i %% 2)==1))
		 	{
				L[i,j]<- -par$psi1
		 	}
		 	if (((i-j)==2) & ((i %% 2)==1))
		 	{
			 	L[i,j]<- -par$phi1
		 	}
		 	if (((i-j)==2) & ((i %% 2)==0))
		 	{
			 	L[i,j]<- -par$phi2
		 	}
		 	if (((i-j)==3) & ((i %% 2)==0))
		 	{
			 	L[i,j]<- -par$psi2
		 	}
		 	if (((i-j)==3) & ((i %% 2)==1))
		 	{
			 	L[i,j]<- 0
		 	}
	   	}
	}

	#get sigepsi
	for (i in 1:T)
	{
		for (j in 1:T){
			if ((i==j) & ((i %% 2)==1))
			{
				sigepsi[i,j]<- par$gama1^2
			}

			if ((i==j) & ((i %% 2)==0))
			{
				sigepsi[i,j]<- par$gama2^2
			}

			if (((i-j)==1) & ((i %% 2)==0))
			{
				sigepsi[i,j]<- par$rho * par$gama1 * par$gama2
			}

			if (((i-j)==1) & ((i %% 2)==1))
			{
				sigepsi[i,j]<-0
			}

			if (((j-i)==1) & ((j %% 2)==0))
			{
				sigepsi[i,j]<- par$rho * par$gama1 * par$gama2
			}

			if (((j-i)==1) & ((j %% 2)==1))
			{
				sigepsi[i,j]<-0
			}
		}
	}

	sigma.inv <- t(L)%*%solve(sigepsi)%*%L;
	sigma.inv <- sigma.inv[c(seq(1, T, 2), seq(2, T, 2)) ,c(seq(1, T, 2), seq(2, T, 2))];

	return(sigma.inv);
}

SAD2.get_mat <- function(par, times, traits=1, options=list())
{
	return(solve(SAD2.get_inv_mat(par, times, traits, options)));
}

SAD2.get_mat_det <- function(par, times, traits=1, options=list())
{
	T <- length(times)*traits;
	par0 <- par;
	if (class(par)!="list")
	{
		par0<-list(rho=par[1], 
			phi1  = par[2],
			phi2  = par[3],
			psi1  = par[4],
			psi2  = par[5],
			gama1 = par[6],
			gama2 = par[7] )
	}
	
	return( ((1-par0$rho^2)*par0$gama1^2*par0$gama2^2)^(T/2) );
}

SAD2.get_init_rand<-function(dat)
{
	par<-runif(7);
	return(list(rho=par[1], 
			phi1  = par[2],
			phi2  = par[3],
			psi1  = par[4],
			psi2  = par[5],
			gama1 = par[6],
			gama2 = par[7]));

	return(list(rho=0.2*runif(1, min=0.9,max=1.1), 
			phi1  = 0.1*runif(1, min=0.9,max=1.1),
			phi2  = 0.2*runif(1, min=0.9,max=1.1),
			psi1  = 0.2*runif(1, min=0.9,max=1.1),
			psi2  = 0.4*runif(1, min=0.9,max=1.1),
			gama1 = 1.41*runif(1, min=0.9,max=1.1),
			gama2 = 0.01*runif(1, min=0.9,max=1.1)));

}

SAD2.is_valid<-function(par)
{
	if (par[1]>0.999 || par[1]<0 )
		return(FALSE);

	return(TRUE);
}

SAD2.get_sum_par <- function(par, fmt)
{
	return("");
}


SAD2.get_est_param<-function( dat, par.cross, par.covar, par.model)
{
	if (!SAD2.is_valid(par.covar))
		return(NaN);

	qtl.prob <- SM.cross$get_qtl_prob( dat$genos_table[,1], NULL, 0,0, par.cross );
	time.std <- dat$sample_times;
	y <- as.matrix( dat$phenos_table );

	mle<-function(par)
	{
		if (!SAD2.is_valid(par))
			return(NaN);
		
		sig.inv <- SAD2.get_inv_mat(par, dat$sample_times, FM.curve$trait_num);
		sig.det <- SAD2.get_mat_det(par, dat$sample_times, FM.curve$trait_num);

		A <- FM.curve$mlefunc(y, par.model, sig.inv, sig.det, qtl.prob, time.std)
		return(A$val);
	}

	h1<-NaN;
	h1 <- try( optim( unlist(par.covar), mle, method = "Nelder-Mead", control=list(maxit=200) ) );
	if ( class(h1)!="try-error" && !is.nan(h1$value) )
		return(list(rho=h1$par[1], 
					phi1=h1$par[2], 
					phi2=h1$par[3], 
					psi1=h1$par[4], 
					psi2=h1$par[5], 
					gama1=h1$par[6], 
					gama2=h1$par[7]))
	else
		return(NaN);
}

SAD2.get_par_num<-function(dat)
{
	par_num = 7;
	return(par_num);
}

SAD2.guess_simu_sigma<-function( sim.mu, trait_num, prob, H2 )
{
	n.len <- dim(sim.mu)[2];
	n.trait <- dim(sim.mu)[2]%/%trait_num;
	
	phe.max <- array(0, dim=c(3, trait_num));
	for (i in 1:trait_num)
	{
		s1<- sim.mu[1, ((i-1)*n.trait+1):(i*n.trait) ] - sim.mu[2, ((i-1)*n.trait+1):(i*n.trait) ]; 
		s2<- sim.mu[2, ((i-1)*n.trait+1):(i*n.trait) ] - sim.mu[3, ((i-1)*n.trait+1):(i*n.trait) ];
		i.max <- which.max(s1*s1+s2*s2);
		phe.max[,i] <- sim.mu[, (i-1)*n.trait+i.max];
	}

	a <- phe.max[,1]
	b <- phe.max[,2]
	q<-prob;

	exQ  <- q^2*a[1]+2*q*(1-q)*a[2]+(1-q)^2*a[3]
	exQ2 <- q^2*a[1]^2+2*q*(1-q)*a[2]^2+(1-q)^2*a[3]^2
	varQ <- exQ2-exQ^2
	exC  <- q^2*b[1]+2*q*(1-q)*b[2]+(1-q)^2*b[3]
	exC2 <- q^2*b[1]^2+2*q*(1-q)*b[2]^2+(1-q)^2*b[3]^2
	varC <- exC2-exC^2

	phi1<-0.1  ###
	phi2<-0.2  ###

	veQH1<- (1-H2)/H2*varQ
	gama1.sq<- veQH1*(1-phi1^2)/(1-phi1^n.len)
	veCH1<- (1-H2)/H2*varC
	gama2.sq<-veCH1*(1-phi2^2)/(1-phi2^n.len)

	return(c(sqrt(gama1.sq), sqrt(gama2.sq)) );
}

covar_SAD2<-list(
	name 	= "SAD2",
	desc 	= "Structured antedependence model (2)",
	is_valid= SAD2.is_valid,
	get_par_num = SAD2.get_par_num, 
	get_mat = SAD2.get_mat,
	get_est_param = SAD2.get_est_param,
	get_init_rand= SAD2.get_init_rand,
	get_mat_det = SAD2.get_mat_det,
	get_inv_mat = SAD2.get_inv_mat,
	get_sum_par = SAD2.get_sum_par,
	guess_simu_sigma = SAD2.guess_simu_sigma);


#COVAR_SAD2 <<- FM2.reg_covar( covar_SAD2 );


