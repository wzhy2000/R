#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
#  Autoregressive model(1)
#  
#  Variance : AR1
#  
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

AR1.get_mat <- function(par0, times, traits=1, options=list())
{
	par<-par0;
	if (class(par0)=="list")
		par <- unlist(par0);

	t_len <- length( times );

	Ar.1 <- array(0, dim=c(t_len*traits,t_len*traits));
	for (i0 in 1:traits)
	for (i1 in 1:traits)
	{
		if (i0==i1)
			for (k0 in 1:t_len)
			for (k1 in 1:t_len)
			{
					Ar.1[(i0-1)*t_len+k0,(i1-1)*t_len+k1] <- par[i0*2]^2 * par[i0*2-1]^abs( k0 - k1 );
			}
	}

	return(Ar.1);
}

AR1.get_inv_mat <- function(par, times, traits=1, options=list())
{
	Ar.1 <- AR1.get_mat(par, times, traits, options)
	return(solve(Ar.1));
}

AR1.get_mat_det <- function(par, times, traits=1, options=list())
{
	Ar.1 <- AR1.get_mat(par, times, traits, options)
	return(det(Ar.1));
}

AR1.get_init_rand<-function(dat)
{
	t.last<-dat$trait_num*(length(dat$sample_times)-1)
	dat.t<-dat$phenos_table[,t.last+c(1:dat$trait_num)];

	if (dat$trait_num>=2)
		par.s2 <- apply(dat.t, 2, function(x){sd(x)})
	else
		par.s2 <- sd(dat.t);

	par <-c();
	for(i in 1:dat$trait_num)
		par<-c(par, runif(1), par.s2[i]);

	
	return(par);
	#return(par_STEM_AR1_NP1$simu_covar);
}

AR1.is_valid<-function(par)
{
	par.rho<-par[seq(1,length(par),2)];
	par.s2<-par[seq(2,length(par),2)];
	over.n1 <- which(par.rho<= 0.000001 | par.rho>0.99)

	return(length(c(over.n1))==0);
}

AR1.get_sum_par <- function(par, fmt)
{
	return("");
}


AR1.get_est_param<-function( dat, par.cross, par.covar, par.model)
{
	if (!AR1.is_valid(par.covar))
		return(NaN);

	qtl.prob <- SM.cross$get_qtl_prob( dat$genos_table[,1], NULL, 0,0, par.cross );
	time.std <- dat$sample_times;
	y <- as.matrix( dat$phenos_table );

	mle<-function(par)
	{
		if (!AR1.is_valid(par))
			return(NaN);
		
		sig.mat<-AR1.get_mat(par, dat$sample_times, FM2.model$trait_num)
		
		sig.inv <- solve(sig.mat);
		sig.det <- det(sig.mat);

		A <- FM.curve$mlefunc(y, par.model, sig.inv, sig.det, qtl.prob, time.std)
		return(A$val);
	}

	h1<-NaN;
	h1 <- try( optim( unlist(par.covar), mle, method = "Nelder-Mead", control=list(maxit=200) ) );
	if ( class(h1)!="try-error" && !is.nan(h1$value) )
		return(h1$par)
	else
		return(NaN);
}

AR1.get_par_num<-function(dat)
{
	if (class(dat)!="list")
		par_num<-2
	else		
		par_num = 2*dat$trait_num;
	return(par_num);
}


AR1.guess_simu_sigma<-function( sim.mu, trait_num, prob, H2 )
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

	mu <- apply(phe.max, 2, function(x){x[1]/2+x[3]/2} );
	a  <- phe.max[1,]-mu;
	d  <- phe.max[2,]-mu;

	q <- 1-prob;

	alpha <- a+(q-prob)*d;
	sig.a <- 2*prob*q*alpha;
	var_a <- 2*prob*q*alpha^2;
	var_d <- (2*prob*q*d)^2;
	var_g <- var_a + var_d;
	var_e <- (1/H2-1)*var_g;

	return(sqrt(var_e));
}

covar_AR1<-list(
	name 	= "AR1",
	desc 	= "parametric stationary autoregressive(AR) model",
	get_par_num  = AR1.get_par_num, 
	is_valid     = AR1.is_valid,
	get_mat      = AR1.get_mat,
	get_inv_mat  = AR1.get_inv_mat,
	get_mat_det  = AR1.get_mat_det,
	get_est_param= AR1.get_est_param,
	get_init_rand= AR1.get_init_rand,
	get_sum_par  = AR1.get_sum_par,
	guess_simu_sigma  = AR1.guess_simu_sigma);

ZZZ.regcovar<-function()
{
	COVAR_AR1  <<- FM2.reg_covar( covar_AR1 );
	COVAR_SAD2 <<- FM2.reg_covar( covar_SAD2 );	
}



