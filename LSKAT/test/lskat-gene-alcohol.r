
#snp.sect <- as.integer( get_con_param("snp.sect") );
#if (is.null(snp.sect)) snp.sect<-NA;

#snp.len <- as.integer( get_con_param("snp.len") );
#if (is.null(snp.len)) snp.len<-NA;

snp.sect<-0;
snp.len<-10;

file.gene.set  <- "/home/zw224/g/rawdata/skat-gen-hg19.SetID"
file.phe.long  <- "/home/zw224/g/rawdata/VACS-960-phe.log.nolead.csv"
file.phe.cov   <- "/home/zw224/g/lskat-gene-cov8-201407/data/plink-r10-pp124-cov.csv"

g.rare <- c(0.5,0.5);
g.nCov <- 8;
g.bTime<- 0;

str.pwd="/home/zw224/g/lskat-gene-cov8-201407/pp124-t0123-r2-beta0.5"
#str.pwd <- getwd();

t <- strsplit(str.pwd, "/")[[1]];
str.p0 <- t[length(t)];
tag.cmd <- strsplit(str.p0, "-")[[1]];

if(tag.cmd[4]=="beta25") g.rare<-c(1,25);
if(tag.cmd[4]=="beta5") g.rare<-c(1,5);

file.plink <- "";
if(tag.cmd[1]=="pp124" && tag.cmd[3]=="r10") file.plink<-"/home/zw224/g/lskat-gene-cov8-201407/data/plink-hwe6-r10-pp124"
if(tag.cmd[1]=="pp124" && tag.cmd[3]=="r5") file.plink<-"/home/zw224/g/lskat-gene-cov8-201407/data/plink-hwe6-r5-pp124"
if(tag.cmd[1]=="pp124" && tag.cmd[3]=="r2") file.plink<-"/home/zw224/g/lskat-gene-cov8-201407/data/plink-hwe6-r2-pp124"
if(tag.cmd[1]=="pp2" && tag.cmd[3]=="r10") file.plink<-"/home/zw224/g/lskat-gene-cov8-201407/data/plink-hwe6-r10-pp2"
if(tag.cmd[1]=="pp2" && tag.cmd[3]=="r5") file.plink<-"/home/zw224/g/lskat-gene-cov8-201407/data/plink-hwe6-r5-pp2"
if(tag.cmd[1]=="pp2" && tag.cmd[3]=="r2") file.plink<-"/home/zw224/g/lskat-gene-cov8-201407/data/plink-hwe6-r2-pp2"

if(tag.cmd[1]=="pp2") file.phe.cov   <- "/home/zw224/g/lskat-gene-cov8-201407/data/plink-r10-pp2-cov.csv"

file.plink.bed <- paste(file.plink, ".bed",sep="") 
file.plink.fam <- paste(file.plink, ".fam",sep="")
file.plink.bim <- paste(file.plink, ".bim",sep="")
file.ret.rdata <- str.p0;

#only include dnaid and bl,followup data
if(tag.cmd[2]=="t0123")
{
	file.phe.temp  <- "VACS-960-phe.log.nolead.0123.csv"
	tb<-read.csv(file.phe.long);
	write.csv(tb[,c(1:5)], file=file.phe.temp, quote=F, row.names=F);
	file.phe.long  <- file.phe.temp;
}

library(snpStats);
snp.mat <- read.plink( file.plink.bed,  file.plink.bim, file.plink.fam);

tb<-read.csv(file.phe.long);
dna.ids<-snp.mat$fam[,1];
newtb <- c();
for(i in 1:length(dna.ids))
{
	id0 <- which(tb[,1]==dna.ids[i]);
	if (length(id0)<1)
		stop("ERROR ID", dna.ids[i], "\n");
		
		
	newtb <- rbind(newtb, tb[id0[1],]);	
}

file.phe.long<-"VACS-960-phe.log.nolead.fam.csv"
write.csv(newtb, file=file.phe.long, quote=F, row.names=F);

snp.sect = 0
snp.len = 100

cat(">>>GENE TEST<<<\n  ", str.p0, "Set=", file.gene.set, "Cov=", file.phe.cov, "Number=", g.nCov, "Time=", g.bTime, "snp.sect=", snp.sect, "snp.len=", snp.len, "\n");

library(LSKAT)

ret<-longskat_gene_plink( file.plink.bed, file.plink.bim, file.plink.fam, file.phe.long, file.phe.cov, 
		file.gene.set, 	gene.range=c(1:340),
		options=list(g.cov= g.nCov, g.btime= g.bTime, g.maxiter = 20, w.common=c(0.5,0.5), w.rare=g.rare, debug=T, run.cpp=T, n.cpu=6 ) )

save(ret, file=file.ret.rdata);

summary(ret)

plot(ret)