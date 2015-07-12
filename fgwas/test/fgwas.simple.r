library(fgwas)

r.simu <- fg.simulate( file.prefix="test.fgwas.simu", "Logistic", "AR1",  n.obs=500, n.snp=100, time.points=c(1:7), par0=NULL, par1=NULL, par2=NULL, par.covar=NULL, par.X=NULL, phe.missing=0.03, snp.missing=0.03, sig.pos=40, plink.format=TRUE );

ret0 <- fgwas.plink(r.simu$file.simple.phe, r.simu$file.plink.bed, r.simu$file.plink.bim, r.simu$file.plink.fam, Y.prefix="Y", Z.prefix="Z", covar.names=c("X_1","X_2"), fgwas.filter = F );

summary(ret0);
plot(ret0);

fg.simulate( file.prefix="test.fgwas.simu", "Logistic", "AR1",  n.obs=500, n.snp=100, time.points=c(1:7), par0=NULL, par1=NULL, par2=NULL, par.covar=NULL, par.X=NULL, phe.missing=0.03, snp.missing=0.03, sig.pos=40, plink.format=FALSE );

ret1 <- fgwas.simple(r.simu$file.simple.phe, r.simu$file.simple.snp, Y.prefix="Y", Z.prefix="Z", covar.names=c("X_1","X_2"), fgwas.filter = F );	

save(ret1, file="test-fgwas-simple.rdata");
summary(ret1);
plot(ret1, fig.prefix="test-fgwas-simple-ret1");
