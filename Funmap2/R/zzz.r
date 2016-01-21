.onAttach<- function(libname, pkgName)
{
    date <- date();
    x <- regexpr("[0-9]{4}", date);
    yr <- substr(date, x[1], x[1] + attr(x, "match.length") - 1);
 
    packageStartupMessage("||");
    packageStartupMessage("|| Funmap2 Package v.2.5");
    packageStartupMessage("|| Build date: ", date(), "");
    packageStartupMessage("|| Copyright (C) 2011-", yr, ", http://statgen.psu.edu", sep="");
    packageStartupMessage("|| Written by Zhong Wang(wzhy2000@hotmail.com)");
    packageStartupMessage("||");
}

.onLoad<- function(libname, pkgName)
{
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

