library(fgwas)

#file.plink.bed <- "/home/zw224/f/bmi/FHS-bmi-v1.bed"  
#file.plink.bim <- "/home/zw224/f/bmi/FHS-bmi-v1.bim"  
#file.plink.fam <- "/home/zw224/f/bmi/FHS-bmi-v1.fam"  

#Yale BulldongN
#file.plink.bed <- "/home/zw224/f/bmi/FHS-bmi-v1-chr2.bed"  
#file.plink.bim <- "/home/zw224/f/bmi/FHS-bmi-v1-chr2.bim"  
#file.plink.fam <- "/home/zw224/f/bmi/FHS-bmi-v1-chr2.fam"  

#TACC Stampede
file.plink.bed <- "/work/03350/tg826494/test/bmidata/FHS-bmi-v1-chr22.bed"  
file.plink.bim <- "/work/03350/tg826494/test/bmidata/FHS-bmi-v1-chr22.bim"  
file.plink.fam <- "/work/03350/tg826494/test/bmidata/FHS-bmi-v1-chr22.fam"  

#file.phe.long  <- "/home/zw224/f/bmi/bmi-phenos-sex-longtime.csv"  
file.phe.long  <- "/work/03350/tg826494/test/bmidata/bmi-phenos-sex-longtime.csv"  

ret1<-ret2<-c();

ret1 <- fgwas.plink( file.phe.long, file.plink.bed, file.plink.bim, file.plink.fam, Y.prefix="Y", Z.prefix="Z", covar.names=c("X"), fgwas.filter = T , options=list(nParallel.cpu=7) );	

save( ret1, ret2, file="fgwas-plink.rdata");

summary(ret1);

plot(ret1, fig.prefix="fgwas-plink");

ret2 <- gls.plink( file.phe.long, file.plink.bed, file.plink.bim, file.plink.fam, Y.prefix="Y", Z.prefix="Z", covar.names=NULL, fgwas.filter = F  , options=list(nParallel.cpu=7) );	


save( ret1, ret2, file="fgwas-plink.rdata");

summary(ret2);

plot(ret2);

