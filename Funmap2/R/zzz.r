msg <- function(...)
{
    date <- date();
    x <- regexpr("[0-9]{4}", date);
    yr <- substr(date, x[1], x[1] + attr(x, "match.length") - 1);
 
    cat("||\n");
    cat("|| Funmap Package v.2.3-1\n");
    cat("|| Build date: ", date(), "\n");
    cat("|| Copyright (C) 2011-", yr, ", http://statgen.psu.edu\n", sep="");
    cat("|| Written by Zhong Wang(zhong.wang@yale.edu)\n");
    cat("||\n");
}

.onAttach<- function(libname, pkgName)
{
	msg();

FM_sys <<- NULL;
FM2.curves   <<- list();
FM2.covars   <<- list();
FM2.crosss   <<- list();
FM2.curve    <<- NULL;
FM2.covar    <<- NULL;
FM2.cross    <<- NULL;

	FM2.start();
	ZZZ.regcovar();
	ZZZ.regcross();
	ZZZ.regcurve();
	ZZZ.regmodel();
}

