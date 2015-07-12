source("fg.main.R");

mu.f <- fg_get_mufunc( CURVE_LOG );

mu.f$par <- c(21.00, 0.9, 0.14);
mu.f$QQ.par <- c(20.00, 0.8, 0.3);
mu.f$Qq.par <- c(21.00, 0.9, 0.2);
mu.f$qq.par <- c(24.00, 0.7, 0.4);

cov.f <- fg_get_covariance( COV_AR1 );
cov.f$par <- c(0.8, 2);

rdata.scan <- "demo3.simu.rdata";
file.simu <- "demo3.simu";

r<-fg_simulate_test( mu.f, cov.f, 200, 30, 1:8, file.simu, rdata.scan, n.perm=100 )
