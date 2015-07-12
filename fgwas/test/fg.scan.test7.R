#1
if(1)
{
	mu.f <- fg_get_mufunc( CURVE_LOG );

	mu.f$par <- c(21.00, 0.9, 0.14);
	mu.f$QQ.par <- c(20.00, 0.8, 0.3);
	mu.f$Qq.par <- c(21.00, 0.9, 0.2);
	mu.f$qq.par <- c(24.00, 0.7, 0.4);

	cov.f <- fg_get_covariance( COV_AR1 );
	cov.f$par <- c(0.8, 2);

	rdata.scan <- "demo3.simu(1).rdata";
	file.simu <- "demo3(1).simu";

	r<-fg_simulate_test( mu.f, cov.f, 200, 30, 1:8, file.simu, rdata.scan, n.perm=100 )

	save(r,file="F1-7.rdata");
}

#2
if(1)
{
	mu.f <- fg_get_mufunc( CURVE_LOG2 );

	mu.f$par <- c(21.00, 0.9, 0.14);
	mu.f$QQ.par <- c(20.00, 0.8, 0.3);
	mu.f$Qq.par <- c(21.00, 0.9, 0.2);
	mu.f$qq.par <- c(24.00, 0.7, 0.4);

	cov.f <- fg_get_covariance( COV_AR1 );
	cov.f$par <- c(0.8, 2);

	rdata.scan <- "demo3.simu(2).rdata";
	file.simu <- "demo3(2).simu";

	r<-NA;
	try(r<-fg_simulate_test( mu.f, cov.f, 200, 30, 1:8, file.simu, rdata.scan, n.perm=100 ));
	str(r);
}

#3
if(1)
{
	mu.f <- fg_get_mufunc( CURVE_LOG );

	mu.f$par <- c(21.00, 0.9, 0.14);
	mu.f$QQ.par <- c(20.00, 0.8, 0.3);
	mu.f$Qq.par <- c(21.00, 0.9, 0.2);
	mu.f$qq.par <- c(24.00, 0.7, 0.4);

	cov.f <- fg_get_covariance( COV_SAD );
	cov.f$par <- c(0.8, 2);

	rdata.scan <- "demo3.simu(3).rdata";
	file.simu <- "demo3(3).simu";

	r<-NA;
	try(r<-fg_simulate_test( mu.f, cov.f, 200, 30, 1:8, file.simu, rdata.scan, n.perm=100 ));
	str(r);
}

#4
if(1)
{
	mu.f <- fg_get_mufunc( CURVE_LOG );

	mu.f$par <- c(21.00, 0.9, 0.14);
	mu.f$QQ.par <- c(20.00, 0.8, 0.3);
	mu.f$Qq.par <- c(21.00, 0.9, 0.2);
	mu.f$qq.par <- c(24.00, 0.7, 0.4);

	cov.f <- fg_get_covariance( COV_AR1 );
	cov.f$par <- c(0.8, 2);

	rdata.scan <- "demo3.simu(4).rdata";
	file.simu <- "demo3(4).simu";

	r<-NA;
	try(r<-fg_simulate_test( mu.f, cov.f, 1000, 30, 1:8, file.simu, rdata.scan, n.perm=100 ));
	str(r);
}

#5
if(1)
{
	mu.f <- fg_get_mufunc( CURVE_LOG );

	mu.f$par <- c(21.00, 0.9, 0.14);
	mu.f$QQ.par <- c(20.00, 0.8, 0.3);
	mu.f$Qq.par <- c(21.00, 0.9, 0.2);
	mu.f$qq.par <- c(24.00, 0.7, 0.4);

	cov.f <- fg_get_covariance( COV_AR1 );
	cov.f$par <- c(0.8, 2);

	rdata.scan <- "demo3.simu(5).rdata";
	file.simu <- "demo3(5).simu";

	r<-NA;
	try(r<-fg_simulate_test( mu.f, cov.f, 200, 300, 1:8, file.simu, rdata.scan, n.perm=100 ));
	str(r);
}

#6
if(1)
{
	mu.f <- fg_get_mufunc( CURVE_LOG );

	mu.f$par <- c(21.00, 0.9, 0.14);
	mu.f$QQ.par <- c(20.00, 0.8, 0.3);
	mu.f$Qq.par <- c(21.00, 0.9, 0.2);
	mu.f$qq.par <- c(24.00, 0.7, 0.4);

	cov.f <- fg_get_covariance( COV_AR1 );
	cov.f$par <- c(0.8, 2);

	rdata.scan <- "demo3.simu(6).rdata";
	file.simu <- "demo3(6).simu";

	r<-NA;
	try(r<-fg_simulate_test( mu.f, cov.f, 200, 30, 1:9, file.simu, rdata.scan, n.perm=100 ));
	str(r);
}

#7
if(1)
{
	mu.f <- fg_get_mufunc( CURVE_LOG );

	mu.f$par <- c(21.00, 0.9, 0.14);
	mu.f$QQ.par <- c(20.00, 0.8, 0.3);
	mu.f$Qq.par <- c(21.00, 0.9, 0.2);
	mu.f$qq.par <- c(24.00, 0.7, 0.4);

	cov.f <- fg_get_covariance( COV_AR1 );
	cov.f$par <- c(0.8, 2);

	rdata.scan <- "demo3.simu(7).rdata";
	file.simu <- "demo3(7).mm";

	r<-NA;
	try(r<-fg_simulate_test( mu.f, cov.f, 200, 30, 1:8, file.simu, rdata.scan, n.perm=100 ));
	str(r);
}

#8
if(1)
{
	mu.f <- fg_get_mufunc( CURVE_LOG );

	mu.f$par <- c(21.00, 0.9, 0.14);
	mu.f$QQ.par <- c(20.00, 0.8, 0.3);
	mu.f$Qq.par <- c(21.00, 0.9, 0.2);
	mu.f$qq.par <- c(24.00, 0.7, 0.4);

	cov.f <- fg_get_covariance( COV_AR1 );
	cov.f$par <- c(0.8, 2);

	rdata.scan <- "demo3.simu(8).wmv";
	file.simu <- "demo3(8).simu";

	r<-NA;
	try(r<-fg_simulate_test( mu.f, cov.f, 200, 30, 1:8, file.simu, rdata.scan, n.perm=100 ));
	str(r);
}

#9
if(1)
{
	mu.f <- fg_get_mufunc( CURVE_LOG );

	mu.f$par <- c(21.00, 0.9, 0.14);
	mu.f$QQ.par <- c(20.00, 0.8, 0.3);
	mu.f$Qq.par <- c(21.00, 0.9, 0.2);
	mu.f$qq.par <- c(24.00, 0.7, 0.4);

	cov.f <- fg_get_covariance( COV_AR1 );
	cov.f$par <- c(0.8, 2);

	rdata.scan <- "demo3.simu(9).rdata";
	file.simu <- "demo3(9).simu";

	r<-NA;
	try(r<-fg_simulate_test( mu.f, cov.f, 200, 30, 1:8, file.simu, rdata.scan,gloop=1, n.perm=100 ));
	str(r);
}

#10
if(1)
{
	mu.f <- fg_get_mufunc( CURVE_LOG );

	mu.f$par <- c(21.00, 0.9, 0.14);
	mu.f$QQ.par <- c(20.00, 0.8, 0.3);
	mu.f$Qq.par <- c(21.00, 0.9, 0.2);
	mu.f$qq.par <- c(24.00, 0.7, 0.4);

	cov.f <- fg_get_covariance( COV_AR1 );
	cov.f$par <- c(0.8, 2);

	rdata.scan <- "demo3.simu(10).rdata";
	file.simu <- "demo3(10).simu";

	r<-NA;
	try(r<-fg_simulate_test( mu.f, cov.f, 200, 30, 1:8, file.simu, rdata.scan, n.perm=100 ,sig.snp=201));
	str(r);
}

#11
if(1)
{
	mu.f <- fg_get_mufunc( CURVE_LOG );

	mu.f$par <- c(21.00, 0.9, 0.14);
	mu.f$QQ.par <- c(20.00, 0.8, 0.3);
	mu.f$Qq.par <- c(21.00, 0.9, 0.2);
	mu.f$qq.par <- c(24.00, 0.7, 0.4);

	cov.f <- fg_get_covariance( COV_AR1 );
	cov.f$par <- c(0.8, 2);

	rdata.scan <- "demo3.simu(1).rdata";
	file.simu <- "demo3(1).simu";

	r<-NA;
	try(r<-fg_simulate_test( mu.f, cov.f, 200, 30, 1:8, file.simu, rdata.scan, n.perm=c(2:50)));
	str(r);
}