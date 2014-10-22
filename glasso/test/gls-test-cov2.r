library(glasso)

phe.out <- "gls.phe.txt"  
snp.out <- "gls.snp.txt"

a_effect <- array(0, dim=c(3,5));
d_effect <- array(0, dim=c(3,5));
cov_effect <- array(0, dim=c(2,4));

sigsnp <- ceiling(runif(5,1,300));
a_effect[1,]<-c( sigsnp[1], 1.04, 0.885, -2.055, 0.545);
a_effect[2,]<-c( sigsnp[2], 1.17, -0.20, 0.74, -4.715);
a_effect[3,]<-c( sigsnp[3], 1.40, -2.25, 1.00,  0.00);

d_effect[1,]<-c( sigsnp[3], 1.49, -2.135, 4.82, 1.425);
d_effect[2,]<-c( sigsnp[4], 1.045, 1.320, 1.905,  1.535);
d_effect[3,]<-c( sigsnp[5], 1.265, -1.225, 2.710, -1.96);

cov_effect[1,]<-c( 2.49, -1.135, 0.82, 0.425);
cov_effect[2,]<-c( -1.045, 2.320, 0.905,  0.535);

gls.simulate( phe.out, snp.out, simu_n= 100, simu_p=300, simu_snp_rho=0.4, simu_rho=0.1, simu_sigma2= 4, 
		 simu_mu= c(13.395, -3.08, 1.875, -3.195),  simu_covar_effect=cov_effect, simu_covar_range=c(-1,1),
		 simu_add_effect=a_effect,  simu_dom_effect=d_effect, 
		 simu_z_range = c(30,60), simu_z_count = c(5,12), debug=F);
	
ret<-gls.simple(phe.out, snp.out, model="y,x_1,x_2,z,add,dom",  bRefit=TRUE, nMaxIter=2000, debug=T);	

save(ret,sigsnp, file="gls-test-cov2.rdata");

