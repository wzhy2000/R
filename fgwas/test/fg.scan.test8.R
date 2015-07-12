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

	file.simu <- "demo2.simu(1)";
	rdata.scan<- "demo2.scan(1).rdata";
	rdata.perm<- "demo2.perm(1).rdata";
	pdf.simu  <- "demo2.simu(1).pdf";

	dat <- fg_simulate_data( mu.f, cov.f, 200, 30, 1:10, file.simu )
	r.scan <- fg_snpscan(d, 1, 30, 2, rdata.scan);
	r.perm <- fg_permutation( d, 10, rdata.perm, n.loop=2 );

	r.perm2 <- fg_perm_merge(r.perm, rdata.perm);
	r.perm2 <- fg_perm_merge(r.perm2, rdata.perm);
	r.perm2 <- fg_perm_merge(r.perm2, rdata.perm);
	r.perm2 <- fg_perm_merge(r.perm2, rdata.perm);
	r.perm2 <- fg_perm_merge(r.perm2, rdata.perm);
	r.perm2 <- fg_perm_merge(r.perm2, rdata.perm);
	r.perm2 <- fg_perm_merge(r.perm2, rdata.perm);
	r.perm2 <- fg_perm_merge(r.perm2, rdata.perm);
	r.perm2 <- fg_perm_merge(r.perm2, rdata.perm);
	r.perm2 <- fg_perm_merge(r.perm2, rdata.perm);
	r.perm2 <- fg_perm_merge(r.perm2, rdata.perm);
	r.perm2 <- fg_perm_merge(r.perm2, rdata.perm);
	r.perm2 <- fg_perm_merge(r.perm2, rdata.perm);
	r.perm2 <- fg_perm_merge(r.perm2, rdata.perm);
	r.perm2 <- fg_perm_merge(r.perm2, rdata.perm);

	r.report<-fg_report(dat, r.scan, r.perm2, pdf.simu);

	save(r.report,file="F1-8.rdata");
	str(r.report)
}

#2
if(1)
{
	mu.f <- fg_get_curve( CURVE_LOG );

	mu.f$par <- c(21.00, 0.9, 0.14);
	mu.f$QQ.par <- c(20.00, 0.8, 0.3);
	mu.f$Qq.par <- c(21.00, 0.9, 0.2);
	mu.f$qq.par <- c(24.00, 0.7, 0.4);

	cov.f <- fg_get_covariance( COV_AR1 );
	cov.f$par <- c(0.8, 2);

	file.simu <- "demo2.simu(2)";
	rdata.scan<- "demo2.scan(2).rdata";
	rdata.perm<- "demo2.perm(2).rdata";
	pdf.simu  <- "demo2.simu(2).pdf";

	dat <- c(2:200)
	r.scan <- fg_snpscan(d, 1, 30, 2, rdata.scan);
	r.perm <- fg_permutation( d, 10, rdata.perm, n.loop=2 );

	r.perm2 <- fg_perm_merge(r.perm, rdata.perm);
	r.perm2 <- fg_perm_merge(r.perm2, rdata.perm);
	r.perm2 <- fg_perm_merge(r.perm2, rdata.perm);
	r.perm2 <- fg_perm_merge(r.perm2, rdata.perm);
	r.perm2 <- fg_perm_merge(r.perm2, rdata.perm);
	r.perm2 <- fg_perm_merge(r.perm2, rdata.perm);
	r.perm2 <- fg_perm_merge(r.perm2, rdata.perm);
	r.perm2 <- fg_perm_merge(r.perm2, rdata.perm);
	r.perm2 <- fg_perm_merge(r.perm2, rdata.perm);
	r.perm2 <- fg_perm_merge(r.perm2, rdata.perm);
	r.perm2 <- fg_perm_merge(r.perm2, rdata.perm);
	r.perm2 <- fg_perm_merge(r.perm2, rdata.perm);
	r.perm2 <- fg_perm_merge(r.perm2, rdata.perm);
	r.perm2 <- fg_perm_merge(r.perm2, rdata.perm);
	r.perm2 <- fg_perm_merge(r.perm2, rdata.perm);

	r.report<-NA;
	try(r.report<-fg_report(dat, r.scan, r.perm2, pdf.simu));
	str(r.report)
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

	file.simu <- "demo2.simu(3)";
	rdata.scan<- "demo2.scan(3).rdata";
	rdata.perm<- "demo2.perm(3).rdata";
	pdf.simu  <- "demo2.simu(3).pdf";

	dat <- fg_simulate_data( mu.f, cov.f, 200, 30, 1:10, file.simu )
	r.scan <- matrix(1:6,nrow=2);
	r.perm <- fg_permutation( d, 10, rdata.perm, n.loop=2 );

	r.perm2 <- fg_perm_merge(r.perm, rdata.perm);
	r.perm2 <- fg_perm_merge(r.perm2, rdata.perm);
	r.perm2 <- fg_perm_merge(r.perm2, rdata.perm);
	r.perm2 <- fg_perm_merge(r.perm2, rdata.perm);
	r.perm2 <- fg_perm_merge(r.perm2, rdata.perm);
	r.perm2 <- fg_perm_merge(r.perm2, rdata.perm);
	r.perm2 <- fg_perm_merge(r.perm2, rdata.perm);
	r.perm2 <- fg_perm_merge(r.perm2, rdata.perm);
	r.perm2 <- fg_perm_merge(r.perm2, rdata.perm);
	r.perm2 <- fg_perm_merge(r.perm2, rdata.perm);
	r.perm2 <- fg_perm_merge(r.perm2, rdata.perm);
	r.perm2 <- fg_perm_merge(r.perm2, rdata.perm);
	r.perm2 <- fg_perm_merge(r.perm2, rdata.perm);
	r.perm2 <- fg_perm_merge(r.perm2, rdata.perm);
	r.perm2 <- fg_perm_merge(r.perm2, rdata.perm);

	r.report<-NA;
	try(r.report<-fg_report(dat, r.scan, r.perm2, pdf.simu));
	str(r.report)
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

	file.simu <- "demo2.simu(4)";
	rdata.scan<- "demo2.scan(4).rdata";
	rdata.perm<- "demo2.perm(4).rdata";
	pdf.simu  <- "demo2.simu(4).pdf";

	dat <- fg_simulate_data( mu.f, cov.f, 200, 30, 1:10, file.simu )
	r.scan <- fg_snpscan(d, 1, 30, 2, rdata.scan);
	r.perm <- fg_permutation( d, 10, rdata.perm, n.loop=2 );

	r.perm2 <- fg_perm_merge();

	r.report<-NA;
	try(r.report<-fg_report(dat, r.scan, r.perm2, pdf.simu));
	str(r.report)
}

#5
if(1)
{
	mu.f <- fg_get_curve( CURVE_LOG );

	mu.f$par <- c(21.00, 0.9, 0.14);
	mu.f$QQ.par <- c(20.00, 0.8, 0.3);
	mu.f$Qq.par <- c(21.00, 0.9, 0.2);
	mu.f$qq.par <- c(24.00, 0.7, 0.4);

	cov.f <- fg_get_covariance( COV_AR1 );
	cov.f$par <- c(0.8, 2);

	file.simu <- "demo2.simu(5)";
	rdata.scan<- "demo2.scan(5).rdata";
	rdata.perm<- "demo2.perm(5).rdata";
	pdf.simu  <- "demo2.simu(5).r";

	dat <- fg_simulate_data( mu.f, cov.f, 200, 30, 1:10, file.simu )
	r.scan <- fg_snpscan(d, 1, 30, 2, rdata.scan);
	r.perm <- fg_permutation( d, 10, rdata.perm, n.loop=2 );

	r.perm2 <- fg_perm_merge(r.perm, rdata.perm);
	r.perm2 <- fg_perm_merge(r.perm2, rdata.perm);
	r.perm2 <- fg_perm_merge(r.perm2, rdata.perm);
	r.perm2 <- fg_perm_merge(r.perm2, rdata.perm);
	r.perm2 <- fg_perm_merge(r.perm2, rdata.perm);
	r.perm2 <- fg_perm_merge(r.perm2, rdata.perm);
	r.perm2 <- fg_perm_merge(r.perm2, rdata.perm);
	r.perm2 <- fg_perm_merge(r.perm2, rdata.perm);
	r.perm2 <- fg_perm_merge(r.perm2, rdata.perm);
	r.perm2 <- fg_perm_merge(r.perm2, rdata.perm);
	r.perm2 <- fg_perm_merge(r.perm2, rdata.perm);
	r.perm2 <- fg_perm_merge(r.perm2, rdata.perm);
	r.perm2 <- fg_perm_merge(r.perm2, rdata.perm);
	r.perm2 <- fg_perm_merge(r.perm2, rdata.perm);
	r.perm2 <- fg_perm_merge(r.perm2, rdata.perm);

	r.report<-NA;
	try(r.report<-fg_report(dat, r.scan, r.perm2, pdf.simu));
	str(r.report)
}
