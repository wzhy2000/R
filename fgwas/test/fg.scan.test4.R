#1
if(1)
{
	mu.f <- fg_get_curve( CURVE_LOG );

	mu.f$par <- c(21.00, 0.9, 0.14);
	mu.f$QQ.par <- c(20.00, 0.8, 0.3);
	mu.f$Qq.par <- c(21.00, 0.9, 0.2);
	mu.f$qq.par <- c(24.00, 0.7, 0.4);

	cov.f <- fg_get_covariance( COV_AR1 );
	cov.f$par <- c(0.8, 2);

	file.simu <- "demo2.simu";
	rdata.scan<- "demo2.scan.rdata";
	rdata.perm<- "demo2.perm1.rdata";
	dat <- fg_simulate_data( mu.f, cov.f, 200, 30, 1:10, file.simu )
	r.perm <- fg_permutation( dat, 100, rdata.perm, n.loop=2 );
	r.perm2 <- fg_perm_merge(r.perm, rdata.perm);

	save(dat,r.perm,r.perm2,file="F1-4.rdata");
	str(r.perm)
	str(r.perm2)
}

#2
if(1)
{
	rdata.perm<- "demo2.perm1.rdata";
	r.perm <- matrix(1:6,nrow=2);

	r.perm2<-NA;
	try(r.perm2 <- fg_perm_merge(r.perm, rdata.perm));
	str(r.perm2)

}

#3
if(1)
{
	mu.f <- fg_get_curve( CURVE_LOG );

	mu.f$par <- c(21.00, 0.9, 0.14);
	mu.f$QQ.par <- c(20.00, 0.8, 0.3);
	mu.f$Qq.par <- c(21.00, 0.9, 0.2);
	mu.f$qq.par <- c(24.00, 0.7, 0.4);

	cov.f <- fg_get_covariance( COV_AR1 );
	cov.f$par <- c(0.8, 2);

	file.simu <- "demo2.simu";
	rdata.scan<- "demo2.scan.rdata";
	rdata.perm<- "demo2.perm3.mm";
	rdata1.perm<-"demo2.permerror.mm";
	dat <- fg_simulate_data( mu.f, cov.f, 200, 30, 1:10, file.simu )
	r.perm <- fg_permutation( dat, 10, rdata.perm, n.loop=2 );
	r.perm2<-NA;
	try(r.perm2 <- fg_perm_merge(r.perm, rdata1.perm));
	str(r.perm)
	str(r.perm2);
}

#4
if(1)
{
	mu.f <- fg_get_curve( CURVE_LOG );

	mu.f$par <- c(21.00, 0.9, 0.14);
	mu.f$QQ.par <- c(20.00, 0.8, 0.3);
	mu.f$Qq.par <- c(21.00, 0.9, 0.2);
	mu.f$qq.par <- c(24.00, 0.7, 0.4);

	cov.f <- fg_get_covariance( COV_AR1 );
	cov.f$par <- c(0.8, 2);

	file.simu <- "demo2.simu";
	rdata.scan<- "demo2.scan.rdata";
	rdata.perm<- c(1:100);
	dat <- fg_simulate_data( mu.f, cov.f, 200, 30, 1:10, file.simu )
	r.perm <- fg_permutation( dat, 10, rdata.perm, n.loop=2 );

	r.perm2<-NA;
	try(r.perm2 <- fg_perm_merge(r.perm, rdata.perm));
	str(r.perm)
	str(r.perm2)
}
