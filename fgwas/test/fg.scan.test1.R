library(fgwas);

#1
if(1)
{
	geno.dat   <- "../rawdata/poplarMF65-201204-geno-v2.dat"
	phe.csv    <- "../rawdata/phe-4.csv" 

	d <- fg_dataest( geno.dat, phe.csv, CURVE_ABRK, COV_SAD, snp.range=c(0, 5) );
str(d)
}

#2
if(1)
{
d<-NA

	geno.dat   <- "../rawdata/geno.txt"
	phe.csv    <- "../rawdata/phe-4.csv" 
	try(d <- fg_dataest( geno.dat, phe.csv, CURVE_ABRK, COV_SAD, snp.range=c(0, 5) ));
str(d)
}

#3
if(1)
{
d<-NA
	geno.dat   <- "../rawdata/poplarMF65-201204-geno-v2.dat"
	phe.csv    <- "../rawdata/phe-wrong1.csv" 
	try(d <- fg_dataest( geno.dat, phe.csv, CURVE_ABRK, COV_SAD, snp.range=c(0, 5) ));
str(d)
}

#4
if(1)
{
d<-NA
	geno.dat   <- "../rawdata/poplarMF65-201204-geno-v2.dat"
	phe.csv    <- "../rawdata/phe-4.csv" 
	try(d <- fg_dataest( geno.dat, phe.csv, CURVE_LOG, COV_SAD, snp.range=c(0, 5)));
str(d)
}

#5
if(1)
{
d<-NA
	geno.dat   <- "../rawdata/poplarMF65-201204-geno-v2.dat"
	phe.csv    <- "../rawdata/phe-4.csv" 
	try(d <- fg_dataest( geno.dat, phe.csv, CURVE_LOG2, COV_SAD, snp.range=c(0, 5)));
str(d)
}

#6
if(1)
{
d<-NA
	geno.dat   <- "../rawdata/poplarMF65-201204-geno-v2.dat"
	phe.csv    <- "../rawdata/phe-4.csv" 
	try(d <- fg_dataest( geno.dat, phe.csv, CURVE_ABRK, COV_AR1, snp.range=c(0, 5)));
str(d)
}

#7
if(1)
{
d<-NA
	geno.dat   <- "../rawdata/poplarMF65-201204-geno-v2.dat"
	phe.csv    <- "../rawdata/phe-4.csv" 
	try(d <- fg_dataest( geno.dat, phe.csv, CURVE_ABRK, COV_SAD, snp.range=c(0, 5),times=c(1:20)));
str(d)
}

#8
if(1)
{
d<-NA;
	geno.dat   <- "../rawdata/poplarMF65-201204-geno-v2.dat"
	phe.csv    <- "../rawdata/phe-4.csv" 
	d <- fg_dataest( geno.dat, phe.csv, CURVE_ABRK, COV_SAD, snp.range=c(0, 5),times=24 );
str(d)
}

#9
if(1)
{
d<-NA;
	geno.dat   <- "../rawdata/poplarMF65-201204-geno-v2.dat"
	phe.csv    <- "../rawdata/phe-4.csv" 
	try(d <- fg_dataest( geno.dat, phe.csv, CURVE_ABRK, COV_SAD, snp.range=c(150000, 160000)));
str(d)
}

