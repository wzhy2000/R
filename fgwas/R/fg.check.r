defmacro <- function(..., expr)
{
	expr <- substitute(expr)
	a <- substitute(list(...))[-1]

	## process the argument list
	nn <- names(a);
	if (is.null(nn)) nn <- rep("", length(a))

	for(i in seq(length=length(a))) 
	{
		if (nn[i] == "") 
		{
			nn[i] <- paste(a[[i]])
			msg <- paste(a[[i]], "not supplied")
			a[[i]] <- substitute(stop(foo),
			list(foo = msg))
		}
	}

	names(a) <- nn
	a <- as.list(a)

	## this is where the work is done
	ff <- eval(substitute(
		function(){
			tmp <- substitute(body)
			eval(tmp, parent.frame())
		},
		list(body = expr)))

	## add the argument list
	formals(ff) <- a
	## create a fake source attribute
	mm <- match.call()
	mm$expr <- NULL
	mm[[1]] <- as.name("macro")

	attr(ff, "source") <- c(deparse(mm),
	deparse(expr))

	## return the fmacrof
	ff
}

check_FG_DAT <- defmacro(fg.dat, expr={
	if (class(fg.dat)!="FG.DAT")
	{
	 	warning("Invalid FG.DAT object(fg.dat).\n")
	 	return(list(error=T, err.info="Invalid FG.DAT object(fg.dat)."))
	}
})
	
check_FG_SCAN <- defmacro(fg.scan, expr={
	if (class(fg.scan)!="FG.SCAN")
	{
		 warning("Invalid FG.SCAN object(fg.scan).\n")
		 return(list(error=T, err.info="Invalid FG.SCAN object(fg.scan)."))
	}
})
	
check_FG_PERM <- defmacro(fg.perm, expr={
	if (class(fg.perm)!="FG.PERM")
	{
		warning("Invalid FG.PERM object(fg.perm).\n")
		return(list(error=T, err.info="Invalid FG.PERM object(fg.perm)."))
	}
})


check_FG_CURVE <- defmacro(fg.curve, expr={
	if (class(fg.curve)!="FG.CURVE")
	{
	 	warning("Invalid data object for the specified curve(fg.curve).\n")
	 	return(list(error=T, err.info="Invalid data object for the specified curve(fg.curve)."))
	}
})

check_FG_COVAR <- defmacro(fg.covar, expr={
	if (class(fg.covar)!="FG.COVAR")
	{
	 	warning("Invalid data object for the specified curve(fg.covar).\n")
	 	return(list(error=T, err.info="Invalid data object for the specified curve(fg.covar)."))
	}
})

check_geno_dat <- defmacro(file.geno.dat, expr={
	if (!is.character(file.geno.dat) || length(file.geno.dat)!=1 || !file.exists(file.geno.dat) )
	{
		warning("Invalid SNP file name(file.geno.dat).\n")
		return(list(error=T, err.info="Invalid SNP file name(file.geno.dat)."))
	}
})

check_phe_csv <- defmacro(file.phe.csv, expr={
	if (!is.character(file.phe.csv) || length(file.phe.csv)!=1 || !file.exists(file.phe.csv) )
	{
		warning("Invalid phenotypic file name(file.phe.csv).\n")
		return(list(error=T,err.info="Invalid phenotypic file name(file.phe.csv)."))
	}
})

check_file_pdf<- defmacro(file.pdf, expr={	
	if (!is.null(file.pdf) && !(is.character(file.pdf) && length(file.pdf)==1) )
	{
	 	warning("Invalid PDF file name(file.pdf).\n")
	 	return(list(error=T, err.info="Invalid PDF file name(file.pdf)."))
	}		
})

check_rdata_perm<- defmacro(rdata.perm, expr={
	if (!is.null(rdata.perm) && !(is.character(rdata.perm) && length(rdata.perm)==1) )
	{
		warning("Invalid Rdata file name for the permutation results.\n");
	 	return(list(error=T, err.info="Invalid Rdata file name for the permutation results."))
	}
})	

check_file_data<- defmacro(file.data, expr={
	if (!is.null(file.data) && !(is.character(file.data) && length(file.data)==1 ) )
	{
	 	warning("Invalid file name for phenotypic data and SNP data(file.data).\n")
	 	return(list(error=T, err.info="Invalid file name for phenotypic data and SNP data(file.data)."))
	}		
})
	
check_file_rdata<- defmacro(file.rdata, expr={
	if (!is.null(file.rdata) && !(is.character(file.rdata) && length(file.rdata)==1 ) )
	{
	 	warning("Invalid file name to save result(file.rdata).\n")
	 	return(list(error=T, err.info="Invalid file name to save result(file.rdata)."))
	}		
})

check_file_rdata_scan<- defmacro(file.rdata.scan, expr={
	if (!is.null(file.rdata.scan) && !(is.character(file.rdata.scan) && length(file.rdata.scan)==1) )
	{
	 	warning("Invalid Rdata file name to save QTL scan result(file.rdata.scan).\n")
	 	return(list(error=T, err.info="Invalid Rdata file name to save QTL scan result(file.rdata.scan)."))
	}
})

check_n_obs <- defmacro(n.obs, expr={
	if (!(is.numeric(n.obs) && length(n.obs)==1) )
	{
	 	warning("Invalid sample size(n.obs).\n")
	 	return(list(error=T, err.info="Invalid sample size(n.obs)."))
	}
})

check_n_snp <- defmacro(n.snp, expr={
	if (!(is.numeric(n.snp) && length(n.snp)==1) )
	{
	 	warning("Invalid SNP count(n.snp).\n")
	 	return(list(error=T, err.info="Invalid SNP count(n.snp)."))
	}
})

check_time_points<- defmacro(time.points, expr={
	if (!(is.null(time.points) || (is.numeric(time.points) && length(time.points)>1)))
	{
		warning("Invalid time points(time.points).\n")
		return(list(error=T,err.info="Invalid time points(time.points)."))
	}
})

check_par_X<- defmacro(par.X, expr={
	if (!(is.numeric(par.X) && length(par.X)>1 || is.null(par.X)) )
	{
		warning("Invalid covariate coefficient(par.X).\n")
		return(list(error=T,err.info="Invalid covariate coefficient(par.X)."))
	}
})

check_snp_range<- defmacro(snp.range, expr={
	if (!(is.numeric(snp.range) && length(snp.range)==2 || is.null(snp.range)) )
	{
		warning("Invalid SNP range(snp.range).\n")
		return(list(error=T,err.info="Invalid SNP range(from, to)."))
	}
})

check_sig_pos<- defmacro(sig.pos, expr={
	if (!(is.numeric(sig.pos) && sig.pos>=1 && sig.pos<=n.snp || is.null(sig.pos)) )
	{
		warning("Invalid significant position(sig.pos).\n")
		return(list(error=T, err.info="Invalid significant position(sig.pos)."))
	}
})
	

check_n_loop<- defmacro(n.loop, expr={
	if (!(is.numeric(n.loop) && length(n.loop)==1) )
	{
		warning("Invalid parameter value(n.loop).\n")
		return(list(error=T,err.info="Invalid parameter value(n.loop)."))
	}
})
	

check_n_perm<- defmacro(n.perm, expr={
	if (!(is.numeric(n.perm)) || length(n.perm)!=1 )
	{
		warning("Invalid permutation count(n.perm).\n")
		return(list(error=T, err.info="Invalid permutation count(n.perm)."))
	}
})

check_n_gloop<- defmacro(n.gloop, expr={
	if (!(is.numeric(n.gloop) && length(n.gloop)==1))
	{
		warning("Invalid option value(n.gloop).\n")
		return(list(error=T, err.info="Invalid option value(n.gloop)."))
	}
})	

check_no_curve<- defmacro(no.curve, expr={
	if (!(is.numeric(no.curve) && length(no.curve)==1))
	{
		warning("Invalid option value(n.curve).\n")
		return(list(error=T, err.info="Invalid option value(n.curve)."))
	}
	
	f.curve <- fg_get_curve( no.curve );
	if ( is.null(f.curve) ) 
	{
		return( list(error=T, err.info=paste("Can't find a curve function[", no.curve, "]" ) ) );
	}
	
})

check_no_covar<- defmacro(no.covar, expr={
	if (!(is.numeric(no.covar) && length(no.covar)==1))
	{
		warning("Invalid option value(n.covar).\n")
		return(list(error=T, err.info="Invalid option value(n.covar)."))
	}
	
	f.covar <- ( no.covar );
	if (is.null(f.covar) ) 
	{
		return( list(error=T, err.info=paste("Can't find a covariance[", no.cov, "]") ) ) ;
	}
	
})

	

