library(fgwas);

if(1)
{
	mu.f <- fg_get_curve( CURVE_LOG );

	mu.f$par <- c(21.00, 0.9, 0.14);
	mu.f$QQ.par <- c(20.00, 0.8, 0.3);
	mu.f$Qq.par <- c(21.00, 0.9, 0.2);
	mu.f$qq.par <- c(24.00, 0.7, 0.4);

	cov.f <- fg_get_covariance( COV_AR1 );
	cov.f$par <- c(0.8, 2);

	file.simu <- "demo2(1).simu";

	dat<- fg_simulate_data (mu.f, cov.f, 200,30,1:10, file.simu ,NA)

	save(dat,file="F1-6.rdata");
	str(dat)
}


#2
if(1)
{
	mu.f <- fg_get_curve( CURVE_LOG2);

	mu.f$par <- c(21.00, 0.9, 0.14);
	mu.f$QQ.par <- c(20.00, 0.8, 0.3);
	mu.f$Qq.par <- c(21.00, 0.9, 0.2);
	mu.f$qq.par <- c(24.00, 0.7, 0.4);

	cov.f <- fg_get_covariance( COV_AR1 );
	cov.f$par <- c(0.8, 2);

	file.simu <- "demo2(2).simu";

	dat<-NA;
	try(dat<- fg_simulate_data (mu.f, cov.f, 200,30,1:10, file.simu ,NA));
	str(dat);
}

#3
if(1)
{
	mu.f <- fg_get_curve( CURVE_LOG );

	mu.f$par <- c(21.00, 0.9, 0.14);
	mu.f$QQ.par <- c(20.00, 0.8, 0.3);
	mu.f$Qq.par <- c(21.00, 0.9, 0.2);
	mu.f$qq.par <- c(24.00, 0.7, 0.4);

	cov.f <- fg_get_covariance( COV_SAD );
	cov.f$par <- c(0.8, 2);

	file.simu <- "demo2(3).simu";

	dat<-NA;
	try(dat<- fg_simulate_data (mu.f, cov.f, 200,30,1:10, file.simu ,NA));
	str(dat);
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

	file.simu <- "demo2(4).simu";

	dat<-NA;
	try(dat<- fg_simulate_data (mu.f, cov.f, -2,30,1:10, file.simu,NA ));
	str(dat);
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

	file.simu <- "demo2(5).simu";

	dat<-NA;
	try(dat<- fg_simulate_data (mu.f, cov.f, 200,-5,1:10, file.simu ,NA));
	str(dat);
}

#6
if(1)
{
	mu.f <- fg_get_curve( CURVE_LOG );

	mu.f$par <- c(21.00, 0.9, 0.14);
	mu.f$QQ.par <- c(20.00, 0.8, 0.3);
	mu.f$Qq.par <- c(21.00, 0.9, 0.2);
	mu.f$qq.par <- c(24.00, 0.7, 0.4);

	cov.f <- fg_get_covariance( COV_AR1 );
	cov.f$par <- c(0.8, 2);

	file.simu <- "demo2(6).simu";

	dat<-NA;
	try(dat<- fg_simulate_data (mu.f, cov.f, 200,30,11, file.simu ,NA));
	str(dat);
}

#7
if(1)
{
	mu.f <- fg_get_curve( CURVE_LOG );

	mu.f$par <- c(21.00, 0.9, 0.14);
	mu.f$QQ.par <- c(20.00, 0.8, 0.3);
	mu.f$Qq.par <- c(21.00, 0.9, 0.2);
	mu.f$qq.par <- c(24.00, 0.7, 0.4);

	cov.f <- fg_get_covariance( COV_AR1 );
	cov.f$par <- c(0.8, 2);

	file.simu <- diag(1,6);

	dat<-NA;
	try(dat<- fg_simulate_data (mu.f, cov.f, 200,30,1:10, file.simu ,NA));
	str(dat);
}

#8
{
	mu.f <- fg_get_curve( CURVE_LOG );

	mu.f$par <- c(21.00, 0.9, 0.14);
	mu.f$QQ.par <- c(20.00, 0.8, 0.3);
	mu.f$Qq.par <- c(21.00, 0.9, 0.2);
	mu.f$qq.par <- c(24.00, 0.7, 0.4);

	cov.f <- fg_get_covariance( COV_AR1 );
	cov.f$par <- c(0.8, 2);

	file.simu <- "demo2(8).simu";

	dat<-NA;
	try(dat<- fg_simulate_data (mu.f, cov.f, 200,30,1:10, file.simu,sig.snp=-7,NA));
	str(dat);
}
