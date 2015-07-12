library(fgwas);

#1
if(1)
{
	g.loop2  <<- get_con_param("loop.max")
	if (is.na(g.loop2))  g.loop2 <<- 10 ;
	snpsect <<- get_con_param("snp.sect")
	if (is.na(snpsect))  snpsect<<-1;
	snplen <<- get_con_param("snp.len")
	if (is.na(snplen))   snplen<<-400;
	perm <<- get_con_param("perm.no")

	cat("g.loop2=", g.loop2, " snpsect=", snpsect, " snplen=", snplen, "\n");

	geno.dat   <- "../rawdata/poplarMF65-201204-geno-v2.dat"
	phe.csv    <- "../rawdata/phe-4.csv" 
	rdata.scan <-  "scan-phe4-sad-demo1.rdata";
	rdata.perm <-  "perm-phe4-sad-demo1.rdata";

	dat <- fg_dataest( geno.dat, phe.csv, CURVE_ABRK, COV_SAD, snp.range=c(0, 10) );
	r.snpscan <- fg_snpscan(dat, NA, NA, g.loop2, rdata.scan);
	summary(r.snpscan);
	r.perm <- fg_permutation( dat, 100, rdata.perm );
	summary(r.perm);

	r.detect<-fg_detect_sig (dat, r.snpscan, r.perm);

	save(r.detect,file="F1-5.rdata");
}

#2
if(1)
{
	g.loop2  <<- get_con_param("loop.max")
	if (is.na(g.loop2))  g.loop2 <<- 10 ;
	snpsect <<- get_con_param("snp.sect")
	if (is.na(snpsect))  snpsect<<-1;
	snplen <<- get_con_param("snp.len")
	if (is.na(snplen))   snplen<<-400;
	perm <<- get_con_param("perm.no")

	cat("g.loop2=", g.loop2, " snpsect=", snpsect, " snplen=", snplen, "\n");

	rdata.scan <-  "scan-phe4-sad-demo2.rdata";
	rdata.perm <-  "perm-phe4-sad-demo2.rdata";

	dat <- c(1:200);
	r.snpscan <- fg_snpscan(dat, NA, NA, g.loop2, rdata.scan);
	summary(r.snpscan);
	r.perm <- fg_permutation( dat, 100, rdata.perm );

	r.detect<-NA;
	try(r.detect<-fg_detect_sig (dat, r.snpscan, r.perm));
	str(r.detect);
}

#3
if(1)
{
	g.loop2  <<- get_con_param("loop.max")
	if (is.na(g.loop2))  g.loop2 <<- 10 ;
	snpsect <<- get_con_param("snp.sect")
	if (is.na(snpsect))  snpsect<<-1;
	snplen <<- get_con_param("snp.len")
	if (is.na(snplen))   snplen<<-400;
	perm <<- get_con_param("perm.no")
	cat("g.loop2=", g.loop2, " snpsect=", snpsect, " snplen=", snplen, "\n");

	geno.dat   <- "../rawdata/poplarMF65-201204-geno-v2.dat"
	phe.csv    <- "../rawdata/phe-4.csv" 
	rdata.scan <-  "scan-phe4-sad-demo3.rdata";
	rdata.perm <-  "perm-phe4-sad-demo3.rdata";

	dat <- fg_dataest( geno.dat, phe.csv, CURVE_ABRK, COV_SAD, snp.range=c(0, 10) );
	r.snpscan <- diag(1,6);
	summary(r.snpscan);
	r.perm <- fg_permutation( dat, 100, rdata.perm );
	summary(r.perm);

	r.detect<-NA;
	try(r.detect<-fg_detect_sig (dat, r.snpscan, r.perm));
	str(r.detect);
}

#4
if(1)
{
	g.loop2  <<- get_con_param("loop.max")
	if (is.na(g.loop2))  g.loop2 <<- 10 ;
	snpsect <<- get_con_param("snp.sect")
	if (is.na(snpsect))  snpsect<<-1;
	snplen <<- get_con_param("snp.len")
	if (is.na(snplen))   snplen<<-400;
	perm <<- get_con_param("perm.no")

	cat("g.loop2=", g.loop2, " snpsect=", snpsect, " snplen=", snplen, "\n");

	geno.dat   <- "../rawdata/poplarMF65-201204-geno-v2.dat"
	phe.csv    <- "../rawdata/phe-4.csv" 
	rdata.scan <-  "scan-phe4-sad-demo4.rdata";
	rdata.perm <-  "perm-phe4-sad-demo4.rdata";

	dat <- fg_dataest( geno.dat, phe.csv, CURVE_ABRK, COV_SAD, snp.range=c(0, 10) );
	r.snpscan <- fg_snpscan(dat, NA, NA, g.loop2, rdata.scan);
	summary(r.snpscan);
	r.perm <- matrix(1:6,nrow=2);
	summary(r.perm);

	r.detect<-NA;
	try(r.detect<-fg_detect_sig (dat, r.snpscan, r.perm));
	str(r.detect);
}