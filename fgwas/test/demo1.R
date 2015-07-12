source("fg.main.R");
source("fg.com.R");

geno.dat   <- "../rawdata/poplarMF65-201204-geno-v2.dat"
phe.csv    <- "../rawdata/phe-4.csv" 
pdf.rpt <-     "scan-phe4-sad-demo.pdf";
rdata.scan <-  "scan-phe4-sad-demo.rdata";
rdata.perm <-  "perm-phe4-sad-demo.rdata";

g.loop2  <<- get_con_param("loop.max")
if (is.na(g.loop2))  g.loop2 <<- 10 ;
snpsect <<- get_con_param("snp.sect")
if (is.na(snpsect))  snpsect<<-1;
snplen <<- get_con_param("snp.len")
if (is.na(snplen))   snplen<<-400;
perm <<- get_con_param("perm.no")

cat("g.loop2=", g.loop2, " snpsect=", snpsect, " snplen=", snplen, "\n");

test1 <- function()
{
	d <- fg_dataest( geno.dat, phe.csv, CURVE_ABRK, COV_SAD, snp.range=c(0, 10) );
	if (d$error)
	{
		cat(d$err.info,"\n");
		return(d$error);
	}
	
	summary(d);
	
	r.scan <- fg_snpscan(d, NA, NA, g.loop2, rdata.scan);

	summary(r.scan);

	r.perm <- fg_permutation( d, 100, rdata.perm );

	summary(r.perm);

	fg_report(d, r.scan, r.perm, pdf.rpt);
}

test1();
