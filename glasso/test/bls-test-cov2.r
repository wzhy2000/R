library(glasso)

phe.out <- "bls.phe.txt"  
snp.out <- "bls.snp.txt"

sigsnp <- c(11, 22, 33, 444, 55);

bls.simulate( phe.out, snp.out, simu_grp=1, simu_n= 400, simu_p=3000, 
		simu_snp_rho = 0.1, 
		simu_rho     = 0.4, 
		simu_sigma2  = 9, 
		simu_mu      = 24, 
		simu_cov_coeff = c( 0, 2 ), 
		simu_a_pos   = c( sigsnp[1], sigsnp[2], sigsnp[3]), 
		simu_a_effect= c( 2.2, -2.5, 2.0 ),  
		simu_d_pos   = c( sigsnp[3], sigsnp[4], sigsnp[5]), 
		simu_d_effect= c( 2.8, 2.0, -2.5 ),
		simu_cov_range=c( 0, 1),
		simu_t_range = c(-1, 1), 
		debug=F );


ret <- bls.simple( phe.out, snp.out, model="y,x_1,x_2,add,dom",  bRefit=TRUE, nMaxIter=2000, fQval.add=0.03, fQval.dom=0.06, debug=T);	

save(ret,sigsnp, file="bls-test-cov2.rdata");

