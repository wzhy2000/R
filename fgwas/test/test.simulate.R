
library(fgwas);

#d1 <- fg.simulate( 600, 1000, CURVE_LOG, COV_AR1, c(1:8), par.X =c(0.8, 2 ), file.prefix="test-simu1" );

fg.simple( "test-simu1.geno.dat",
		"test-simu1.pheY.csv", 
		"test-simu1.pheX.csv", 
		NULL, 
		no.curve=CURVE_LOG, 
		no.covar=COV_AR1, 
		covar.names=c("X_1", "X_2"), 
		fgwas.filter=T, 
		snp.range=NULL,
		n.perm = 0 ) ;
