g.a     <- 0.5;
g.b     <- 0.5;
#g.a <- -0.26;
#g.b <- 10.35;
g.rho   <- 0.75;
g.sig_a <- 0.8;
g.sig_b <- 0.8;
g.sig_e <- 0.8;
g.beta.c1   <- 0.1;
g.snprange  <- c(29*1000, 30*1000)
g.phe.len   <- 6;
g.alevel    <- 10^(-6)

source("long-simu-run.r");

power.test<-function( nloop, nsample, phe.dist, phe.cov, par, power.rdata )
{
	f.simu <- f.mn.ar1;
	if (phe.dist=="mn" && phe.cov=="ar1") f.simu <- f.mn.ar1; 
	if (phe.dist=="mn" && phe.cov=="sad") f.simu <- f.mn.sad; 
	if (phe.dist=="mn" && phe.cov=="cm") f.simu  <- f.mn.cm; 
	if (phe.dist=="mt" ) f.simu <- f.mt.ar1; 
	if (phe.dist=="msn" ) f.simu <- f.sn.ar1; 
	if (phe.dist=="mmn" ) f.simu <- f.mmn.ar1; 
	
	par$beta.effect<- c( 0.1, 0.1, 0, 0.8 ); 
	par$c1<-c(0.2828, 0.2828, 0, 0);
	r.1 <-check.power( nloop, nsample, g.snprange, c(1,25), c(0.5,0.5), par, f.simu )
save(par, r.1, file=power.rdata);

	par$beta.effect<- c( 0.2, 0.2, 0, 0.6 );
	par$c1<-c(0.2, 0.2, 0, 0);
	r.2 <-check.power( nloop, nsample, g.snprange, c(1,25), c(0.5,0.5), par, f.simu )
save(par, r.1, r.2, file=power.rdata);

	par$beta.effect<- c( 0.3, 0.3, 0, 0.4 );
	par$c1<-c(0.1633, 0.1633, 0, 0);
	r.3 <-check.power( nloop, nsample, g.snprange, c(1,25), c(0.5,0.5), par, f.simu )
save(par, r.1, r.2, r.3, file=power.rdata);

	par$beta.effect<- c( 0.1, 0.1, 0.1, 0.1, 0.6 );
	par$c1<-c(0.2, 0.0, -0.2,  0, 0);
	r.4 <-check.power( nloop, nsample, g.snprange, c(1,25), c(0.5,0.5), par, f.simu )
save(par, r.1, r.2, r.3, r.4, file=power.rdata);

	par$beta.effect<- c( 0.2, 0.2, 0.0, 0.6 );
	par$c1<-c(0.24, 0.15, 0, 0);
	r.5 <-check.power( nloop, nsample, g.snprange, c(1,25), c(0.5,0.5), par, f.simu )
save(par, r.1, r.2, r.3, r.4, r.5, file=power.rdata);

	par$beta.effect<- c( 0.2, 0.0, 0.6 );
	par$c1<-c(NA, 0, 0);
	r.6 <-check.power( nloop, nsample, g.snprange, c(1,25), c(0.5,0.5), par, f.simu )

	save(par, r.1, r.2, r.3, r.4, r.5, r.6, file=power.rdata);

	a.level<-10^(-6);
	r.power<-c();
	r.power<-rbind( r.power, c(nsample, 1, length(which(r.1[,2]<a.level)), length(which(r.1[,3]<a.level)), length(which(r.1[,4]<a.level)) ) );
	r.power<-rbind( r.power, c(nsample, 2, length(which(r.2[,2]<a.level)), length(which(r.2[,3]<a.level)), length(which(r.2[,4]<a.level)) ) );
	r.power<-rbind( r.power, c(nsample, 3, length(which(r.3[,2]<a.level)), length(which(r.3[,3]<a.level)), length(which(r.3[,4]<a.level)) ) );
	r.power<-rbind( r.power, c(nsample, 4, length(which(r.4[,2]<a.level)), length(which(r.4[,3]<a.level)), length(which(r.4[,4]<a.level)) ) );
	r.power<-rbind( r.power, c(nsample, 5, length(which(r.5[,2]<a.level)), length(which(r.5[,3]<a.level)), length(which(r.5[,4]<a.level)) ) );
	r.power<-rbind( r.power, c(nsample, 6, length(which(r.6[,2]<a.level)), length(which(r.6[,3]<a.level)), length(which(r.6[,4]<a.level)) ) );
	colnames(r.power) <- c("sample", "beta.plus", "long.skat", "skat.bl", "skat.mu");
	r.power6 <- r.power;
	
	a.level<-10^(-4);
	r.power<-c();
	r.power<-rbind( r.power, c(nsample, 1, length(which(r.1[,2]<a.level)), length(which(r.1[,3]<a.level)), length(which(r.1[,4]<a.level)) ) );
	r.power<-rbind( r.power, c(nsample, 2, length(which(r.2[,2]<a.level)), length(which(r.2[,3]<a.level)), length(which(r.2[,4]<a.level)) ) );
	r.power<-rbind( r.power, c(nsample, 3, length(which(r.3[,2]<a.level)), length(which(r.3[,3]<a.level)), length(which(r.3[,4]<a.level)) ) );
	r.power<-rbind( r.power, c(nsample, 4, length(which(r.4[,2]<a.level)), length(which(r.4[,3]<a.level)), length(which(r.4[,4]<a.level)) ) );
	r.power<-rbind( r.power, c(nsample, 5, length(which(r.5[,2]<a.level)), length(which(r.5[,3]<a.level)), length(which(r.5[,4]<a.level)) ) );
	r.power<-rbind( r.power, c(nsample, 6, length(which(r.6[,2]<a.level)), length(which(r.6[,3]<a.level)), length(which(r.6[,4]<a.level)) ) );
	colnames(r.power) <- c("sample", "beta.plus", "long.skat", "skat.bl", "skat.mu");	
	r.power4 <- r.power;

	a.level<-10^(-2);
	r.power<-c();
	r.power<-rbind( r.power, c(nsample, 1, length(which(r.1[,2]<a.level)), length(which(r.1[,3]<a.level)), length(which(r.1[,4]<a.level)) ) );
	r.power<-rbind( r.power, c(nsample, 2, length(which(r.2[,2]<a.level)), length(which(r.2[,3]<a.level)), length(which(r.2[,4]<a.level)) ) );
	r.power<-rbind( r.power, c(nsample, 3, length(which(r.3[,2]<a.level)), length(which(r.3[,3]<a.level)), length(which(r.3[,4]<a.level)) ) );
	r.power<-rbind( r.power, c(nsample, 4, length(which(r.4[,2]<a.level)), length(which(r.4[,3]<a.level)), length(which(r.4[,4]<a.level)) ) );
	r.power<-rbind( r.power, c(nsample, 5, length(which(r.5[,2]<a.level)), length(which(r.5[,3]<a.level)), length(which(r.5[,4]<a.level)) ) );
	r.power<-rbind( r.power, c(nsample, 6, length(which(r.6[,2]<a.level)), length(which(r.6[,3]<a.level)), length(which(r.6[,4]<a.level)) ) );
	colnames(r.power) <- c("sample", "beta.plus", "long.skat", "skat.bl", "skat.mu");	
	r.power2 <- r.power;
	
	save( r.power6, r.power4, r.power2, par, r.1, r.2, r.3, r.4, r.5, r.6, file=power.rdata);
}

type1.test<-function(nloop, nrep, nsample, phe.dist, phe.cov, par, type1.rdata)
{
	f.simu <- f.mn.ar1;
	if (phe.dist=="mn" && phe.cov=="ar1") f.simu <- f.mn.ar1; 
	if (phe.dist=="mn" && phe.cov=="sad") f.simu <- f.mn.sad; 
	if (phe.dist=="mn" && phe.cov=="cm") f.simu  <- f.mn.cm; 
	if (phe.dist=="mt" ) f.simu <- f.mt.ar1; 
	if (phe.dist=="msn" ) f.simu <- f.sn.ar1; 
	if (phe.dist=="mmn" ) f.simu <- f.mmn.ar1; 

	cat(" [Type1 Test]\n")
	cat(" * loop=", nloop, "\n")
	cat(" * rep=", nrep, "\n")
	cat(" * nsample=", nsample, "\n")
	cat(" * par=",  unlist(par), "\n")
	cat(" * phe.dist=", phe.dist, "\n")
	cat(" * phe.cov=", phe.cov, "\n")
	cat(" * ret.rdata=", type1.rdata, "\n")
	
	r.1  <- check.type1err( nloop, nrep, nsample, g.snprange, c(1,25), c(0.5,0.5), par, f.simu )
	save(par, r.1, file=type1.rdata);

	a5.level<-10^(-5);
	a4.level<-10^(-4);
	a3.level<-10^(-3);
	a2.level<-10^(-2);

	r.power<-c(nloop, nsample, 
	          length(which( r.1[,3]<a2.level)), 
	          length(which( r.1[,3]<a3.level)), 
	          length(which( r.1[,3]<a4.level)), 
	          length(which( r.1[,3]<a5.level)) );

	names(r.power) <- c("repi", "sample", "a2", "a3", "a4", "a5");
	save(r.power, par, r.1, file=type1.rdata);
}


test.main<-function()
{
	power.test <- as.integer(get_con_param("power"))==1
	nloop      <- as.integer(get_con_param("nloop"))
	nrep       <- as.integer(get_con_param("nrep"))
	nsample    <- as.integer(get_con_param("phe.sample"))
	phe.dist   <- as.character(get_con_param("phe.dist"))
	phe.cov    <- as.character(get_con_param("phe.cov"))
	phe.par    <- as.character(get_con_param("phe.par"))
	ret.rdata  <- as.character(get_con_param("ret.rdata"))

	#power.test <- T;
	#nloop      <- 2;
	#nrep       <- 1000;
	#phe.sample <- 500;
	#phe.dist   <- "msn"
	#phe.cov    <- "ar1"
	#phe.par    <- "0.7,0.5"
	#ret.rdata  <- "test-xxx.rdata"

	pars <- strsplit(phe.par, ",");
	par1 <- as.numeric(pars[[1]][1]);
	par2 <- as.numeric(pars[[1]][2]);
	par3 <- as.numeric(pars[[1]][3]);
	par4 <- as.numeric(pars[[1]][4]);
	par5 <- as.numeric(pars[[1]][5]);
	par6 <- as.numeric(pars[[1]][6]);

	par <- list(a=g.a, b=g.b, sig_a=g.sig_a, sig_b=g.sig_b, sig_e=g.sig_e, 
		 c1=g.beta.c1, times=g.phe.len, a.level=g.alevel, 
		 par1=par1, par2=par2, par3=par3, par4=par4, par5=par5, par6=par6 );

	sig_a <- as.numeric(get_con_param("sig.a"))
	sig_b <- as.numeric(get_con_param("sig.b"))
	sig_e <- as.numeric(get_con_param("sig.e"))	
	a     <- as.numeric(get_con_param("a"))
	b     <- as.numeric(get_con_param("b"))
	
	if(!is.na(sig_a)) par$sig_a <- sig_a;
	if(!is.na(sig_b)) par$sig_b <- sig_b;
	if(!is.na(sig_e)) par$sig_e <- sig_e;
	if(!is.na(a))     par$a     <- a;
	if(!is.na(b))     par$b     <- b;

	if ( power.test )
		power.test( nloop, nsample, phe.dist, phe.cov, par, ret.rdata )
	else
		type1.test( nloop, nrep, nsample, phe.dist, phe.cov, par, ret.rdata);
}



# $R  --vanilla --quite power=1 nloop=1000 phe.sample=500  phe.dist=mn phe.cov=ar1 phe.par=0.7,0.8 ret.rdata=power.L1.mn.ar1.500.rdata sig.a=0.8 sig.e=0.8   < long-simu-testplan.r > power-L1-mn-ar1-500.out
# $R  --vanilla --quite power=1 nloop=1000 phe.sample=500  phe.dist=mn phe.cov=ar1 phe.par=0.7,0.4 ret.rdata=power.L2.mn.ar1.500.rdata sig.a=0.4 sig.e=0.4   < long-simu-testplan.r > power-L2-mn-ar1-500.out
# $R  --vanilla --quite power=1 nloop=1000 phe.sample=500  phe.dist=mn phe.cov=ar1 phe.par=0.7,0.2 ret.rdata=power.L3.mn.ar1.500.rdata sig.a=0.2 sig.e=0.2   < long-simu-testplan.r > power-L3-mn-ar1-500.out
# $R  --vanilla --quite power=1 nloop=1000 phe.sample=500  phe.dist=mn phe.cov=ar1 phe.par=0.7,0.8 ret.rdata=power.L4.mn.ar1.500.rdata sig.a=0.2 sig.e=0.2   < long-simu-testplan.r > power-L4-mn-ar1-500.out
# $R  --vanilla --quite power=1 nloop=1000 phe.sample=500  phe.dist=mn phe.cov=ar1 phe.par=0.7,0.2 ret.rdata=power.L5.mn.ar1.500.rdata sig.a=0.8 sig.e=0.8   < long-simu-testplan.r > power-L5-mn-ar1-500.out
# $R  --vanilla --quite power=1 nloop=1000 phe.sample=500  phe.dist=mt phe.cov=ar1 phe.par=0.7,0.8 ret.rdata=power.L6.mt.ar1.500.rdata sig.a=0.8 sig.e=0.8  < long-simu-testplan.r > power-L6-mt-ar1-500.out
# $R  --vanilla --quite power=1 nloop=1000 phe.sample=500  phe.dist=mt phe.cov=ar1 phe.par=0.7,0.4 ret.rdata=power.L7.mt.ar1.500.rdata sig.a=0.4 sig.e=0.4  < long-simu-testplan.r > power-L7-mt-ar1-500.out
# $R  --vanilla --quite power=1 nloop=1000 phe.sample=500  phe.dist=msn phe.cov=ar1 phe.par=0.7,0.8 ret.rdata=power.L8.msn.ar1.500.rdata sig.a=0.8 sig.e=0.8  < long-simu-testplan.r > power-L8-msn-ar1-500.out
# $R  --vanilla --quite power=1 nloop=1000 phe.sample=500  phe.dist=msn phe.cov=ar1 phe.par=0.7,0.4 ret.rdata=power.L9.msn.ar1.500.rdata sig.a=0.4 sig.e=0.4  < long-simu-testplan.r > power-L9-msn-ar1-500.out
# $R  --vanilla --quite power=1 nloop=1000 phe.sample=500  phe.dist=mmn phe.cov=ar1 phe.par=0.5,0.8,0.5,0.8,0.6 ret.rdata=power.L10.mmn.ar1.500.rdata sig.a=0.8 sig.e=0.8  < long-simu-testplan.r > power-L10-mmn-ar1-500.out
# $R  --vanilla --quite power=1 nloop=1000 phe.sample=500  phe.dist=mmn phe.cov=ar1 phe.par=0.5,0.8,0.3,0.4,0.7 ret.rdata=power.L11.mmn.ar1.500.rdata sig.a=0.4 sig.e=0.4  < long-simu-testplan.r > power-L11-mmn-ar1-500.out
# $R  --vanilla --quite power=1 nloop=1000 phe.sample=500  phe.dist=mn phe.cov=sad phe.par=0.7,0.8 ret.rdata=power.L12.mn.sad.500.rdata sig.a=0.8 sig.e=0.8  < long-simu-testplan.r > power-L12-mn-sad-500.out
# $R  --vanilla --quite power=1 nloop=1000 phe.sample=500  phe.dist=mn phe.cov=sad phe.par=0.7,0.4 ret.rdata=power.L13.mn.sad.500.rdata sig.a=0.4 sig.e=0.4  < long-simu-testplan.r > power-L13-mn-sad-500.out
# $R  --vanilla --quite power=1 nloop=1000 phe.sample=500  phe.dist=mn phe.cov=cm phe.par=0.7,0.8 ret.rdata=power.L14.mn.cm.500.rdata sig.a=0.8 sig.e=0.8  < long-simu-testplan.r > power-L14-mn-cm-500.out
# $R  --vanilla --quite power=1 nloop=1000 phe.sample=500  phe.dist=mn phe.cov=cm phe.par=0.7,0.4 ret.rdata=power.L15.mn.cm.500.rdata sig.a=0.4 sig.e=0.4  < long-simu-testplan.r > power-L15-mn-cm-500.out
# $R  --vanilla --quite power=1 nloop=1000 phe.sample=500  phe.dist=mn phe.cov=ar1 phe.par=0.87,1.48 ret.rdata=power.L16.mn.ar1.500.rdata sig.a=0.001 sig.e=1.67   < long-simu-testplan.r > power-L16-mn-ar1-500.out
# 
# 
# $R  --vanilla --quite power=0 nloop=10 nrep=1000 phe.sample=500 phe.dist=mn phe.cov=ar1 phe.par=0.7,0.8 ret.rdata=type1.L1.mn.ar1.500.1.rdata sig.a=0.8 sig.e=0.8  < long-simu-testplan.r > type1-L1-mn-ar1-500-1.out


test.main();
