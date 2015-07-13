#  File src/library/grDevices/R/zzz.R

.noGenerics <- TRUE

msg <- function(...)
{
    date <- date()
    x <- regexpr("[0-9]{4}", date)
    yr <- substr(date, x[1], x[1] + attr(x, "match.length") - 1)
 
    cat("##\n## Funmap Package v.2.2-1\n")
    cat("## Build date: ", date(), "\n")
    cat("## Copyright (C) 2014-", yr, ", http://statgen.psu.edu\n", sep="")
    cat("## Written by Zhong Wang(zw355@cornell.edu)\n\n")
}

.onAttach <- function(...)
{
	#msg();

    cat("========================================\n")
    cat("fGWAS PACKAGE is developed by Zhong Wang.\n")

    LOG_ERR      <<- 0;
    LOG_WARN     <<- 1;
    LOG_INFO     <<- 2;
    LOG_DEBUG    <<- 3;

    fgwas.sys <<-list( curves=list(), covariances=list(), n.seed=1, log=LOG_INFO, try.slient = T, loop.est=10, loop.mle=4, mu.range.loop=100 );
	
	ZZZ.regcurve();
	#ZZZ.regcovar();
    
    cat("========================================\n")
}

.onLoad <- function(libname, pkgname)
{
	#Do Nothing
}


