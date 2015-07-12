library(fgwas);

mu.f <- fg_get_curve( CURVE_LOG );

mu.f$par <- c(21.00, 0.9, 0.14);
mu.f$QQ.par <- c(20.00, 0.8, 0.3);
mu.f$Qq.par <- c(21.00, 0.9, 0.2);
mu.f$qq.par <- c(24.00, 0.7, 0.4);

cov.f <- fg_get_covariance( COV_AR1 );
cov.f$par <- c(0.8, 2);

file.simu <- "demo2.simu";
rdata.scan<- "demo2.scan.rdata";
rdata.perm<- "demo2.perm.rdata";
pdf.simu  <- "demo2.simu.pdf";

dat <- fg_simulate_data( mu.f, cov.f, 200, 30, 1:10, file.simu )
r.scan <- fg_snpscan(dat, 1, 30, 2, rdata.scan, cache.reset=T, debug=F);

r.perm <- fg_permutation( dat, 10, rdata.perm, n.loop=2, cache.reset=F, debug=F );

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

fg_report(dat, r.scan, r.perm2, pdf.simu);
fg_report(dat);
