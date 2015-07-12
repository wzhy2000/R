library(fgwas);

#1
if(1)
{
	geno.dat   <- "../rawdata/poplarMF65-201204-geno-v2.dat"
	phe.csv    <- "../rawdata/phe-4.csv" 
        rdata.perm <-  "perm-phe4-sad-demo.rdata";

	d <- fg_dataest( geno.dat, phe.csv, CURVE_ABRK, COV_SAD, snp.range=c(0, 5) );
        r.perm <- fg_permutation( d, 100, rdata.perm );

save(d,r.perm,file="F1-3.rdata");
str(r.perm)

}

#2
if(1)
{
	r.perm<-NA;
	d <-c(2:50);
        rdata.perm <-  "perm-phe4-sad-demo.rdata";
        try(r.perm <- fg_permutation( d, 100, rdata.perm ));
str(r.perm)
}

#3
if(1)
{
        r.perm<-NA;
	geno.dat   <- "../rawdata/poplarMF65-201204-geno-v2.dat"
	phe.csv    <- "../rawdata/phe-4.csv" 
        rdata.perm <-  matrix(1:6,nrow=2);
	d <- fg_dataest( geno.dat, phe.csv, CURVE_ABRK, COV_SAD, snp.range=c(0, 5) );
        try(r.perm <- fg_permutation( d, 100, rdata.perm ));
str(r.perm)
}