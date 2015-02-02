###############################################################
# 
# New Systems Mapping Application(SysMap1)#
# References:
# 1. Systems mapping: A computational tool for personalized 
#    medicine
#
#
# Routine:
#  1) FM2.set_value
#  2) FM2.get_value
#  3) FM2.get_curve
#  4) FM2.reg_curve
#  5) FM2.param
#  6) summary.FM2.par 
#  7) FM2.simulate
#  8) FM2.load_data
#  9) summary.FM2.dat
# 10) report.FM2.dat
# 11) plot.FM2.dat
# 12) FM2.hp_test
# 13) summary.FM2.ret.hp
# 14) report.FM2.ret.hp
# 15) plot.FM2.ret.hp
# 16) FM2.permutation
# 17) summary.FM2.ret.perm
# 18) plot.FM2.ret.perm
# 19) FM2.evaluate
# 20) summary.FM2.ret.eval
# 21) plot.FM2.ret.eval
# 22) FM2.report
#
# History:
# 03/17/2010 Version 0.1
#
###############################################################

# for rmultnorm function, MSBVAR is necessary.
#library(MSBVAR);

FM2.start <- function()
{
	FM_sys <<- FM2.sys("FM");
	FM_sys$reset();
}

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# usage:
#    par<- FM2.param( par_obj, curve_type, cross_type, covar_type );
#    summary(par);
#
#    dat<- FM2.simulate( par);
#    summary(dat);
#    plot(dat);
#
#    dat<- FM2.load_data( pheno_file, geno_file, marker_file, curve_type, cross_type );
#    summary(dat);
#    plot(dat, plot_type);
#
#    ret<- FM2.hp_test( dat, test=c(1) );
#    summary(ret, dat );
#    plot( ret, dat );
#
#    ret<- FM2.permutation( dat );
#    summary(ret);
#    plot(ret);
#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


#--------------------------------------------------------------
# FM2.set_value
#
#--------------------------------------------------------------
FM2.set_value <- function( key, value)
{
	if (is.null( FM_sys ) )
		FM2.start();
	
	return( FM_sys$set_value(key, value) );
}

FM2.get_value <- function( key, def = NA )
{
	if (is.null( FM_sys ) )
		FM2.start();
	
	r<-NA;
	try( r <- FM_sys$get_value(key), TRUE );
	
	if ( ( is.null(r) || is.na(r) ) && !is.na( def ) )
		return(def);
		
	return(r);
}


FM2.get_curve <- function(curve_type)
{
	for (i in 1:length(FM2.curves))
		if (FM2.curves[[i]]$type==curve_type)
			return(FM2.curves[[i]])
	
	return(NULL);
}

FM2.reg_curve <- function(curve)
{
	if (is.null( FM_sys ) )
		FM2.start();

	cat("Curve Registration:", curve$desc, "(", curve$name, ")\n");

	curve$type <-length(FM2.curves)+1; 
	for( key in names(class_curve) )
	{
		if (!(key %in% names(curve)))
			curve[[key]] <- class_curve[[key]];
			
	}

	FM2.curves[[length(FM2.curves)+1]] <<- curve;
	
	return(length(FM2.curves));
}

FM2.get_cross <- function(cross_type)
{
	for (i in 1:length(FM2.crosss))
		if (FM2.crosss[[i]]$type==cross_type)
			return(FM2.crosss[[i]])
	
	return(NULL);
}

FM2.reg_cross <- function(cross)
{
	if (is.null( FM_sys ) )
		FM2.start();

	cat("Cross Registration:", cross$desc, "(", cross$name, ")\n");

	cross$type <-length(FM2.crosss)+1; 
	FM2.crosss[[length(FM2.crosss)+1]] <<- cross;
	
	return(length(FM2.crosss));
}


FM2.get_covar <- function(covar_type)
{
	for (i in 1:length(FM2.covars))
		if (FM2.covars[[i]]$type==covar_type)
			return(FM2.covars[[i]])
	
	return(NULL);
}

FM2.reg_covar <- function(covar)
{
	if (is.null( FM_sys ) )
		FM2.start();

	cat("Covar Registration:", covar$desc, "(", covar$name, ")\n");

	covar$type <-length(FM2.covars)+1; 
	FM2.covars[[length(FM2.covars)+1]] <<- covar;
	
	return(length(FM2.covars));
}


#--------------------------------------------------------------
# FM2.param
#
# Create a param object for the simulation( LC,BI,PC curve * BC,F2) 
#--------------------------------------------------------------
FM2.param <- function( par_obj, curve_type, cross_type = CROSS_BC, covar_type=COVAR_AR1 )
{
	if (is.null( FM_sys ) ) 
		FM2.start();
	
	curve <- FM2.get_curve(curve_type);
	if (is.null(curve))
	{
		stop("inavlid curve type");
	}
	FM2.curve <<-curve;

	cross <- FM2.get_cross(cross_type);
	if (is.null(cross))
	{
		stop("inavlid cross type");
	}
	FM2.cross <<- cross;

	covar <- FM2.get_covar(covar_type);
	if (is.null(covar))
	{
		stop("inavlid covar type");
	}
	FM2.covar <<- covar;

	par <- dat.get_simuparam( par_obj, cross_type, covar_type );

	class(par)<-"FM2.par";
	return(par);
}

#--------------------------------------------------------------
# summary.FM2.par
#
# used by summary( FM2.par.obj )
#--------------------------------------------------------------
summary.FM2.par <- function( par_obj , file=NA , append = TRUE )
{
	if (is.null(FM2.curve))
	{
		stop("inavlid curve type");
	}

	str<- dat.summary_par( par_obj );

	if (is.na(file))
		return ( cat( str ) )
	else
		return ( cat( str, file= file, append = append) )
}

#--------------------------------------------------------------
# FM2.simulate
#
# Create a simulation data set by parameter object
#--------------------------------------------------------------
FM2.simulate<- function( par_obj )
{
	if ( is.null(FM2.curve) )
	{
		stop("inavlid curve type");
	}

	dat <- dat.simulate( par_obj );
	class(dat) <- "FM2.dat";
	return(dat);
}

#--------------------------------------------------------------
# FM2.load_data
#
# Load a real data, phenotype file, genotype file, marker file 
# is necessary.
#--------------------------------------------------------------
FM2.load_data <- function( pheno_file, geno_file, marker_file, curve_type, cross_type, covar_type=COVAR_AR1  )
{
	if (is.null( FM_sys ) )
		FM2.start();

	curve <- FM2.get_curve(curve_type);
	if (is.null(curve))
	{
		stop("inavlid curve type");
	}
	
	FM2.curve <<- curve;

	cross <- FM2.get_cross(cross_type);
	if (is.null(cross))
	{
		stop("inavlid cross type");
	}
	FM2.cross <<- cross;

	covar <- FM2.get_covar(covar_type);
	if (is.null(covar))
	{
		stop("inavlid covar type");
	}
	FM2.covar <<- covar;

	dat = dat.load( pheno_file, geno_file, marker_file, cross_type, covar_type );
	class(dat)<-"FM2.dat";

	return(dat);
}

#--------------------------------------------------------------
# summary.FM2.dat
#
# used by summary( FM2.dat.obj )
#--------------------------------------------------------------
summary.FM2.dat<- function( dat_obj, file=NA , append = TRUE )
{
	if (is.null(FM2.curve))
	{
		stop("inavlid curve type");
	}

	str<- dat.summary( dat_obj );
	if (is.na(file))
		return ( cat( str ) )
	else
		return ( cat( str, file= file, append = append) )
}

#--------------------------------------------------------------
# report.FM2.dat
#
# used by report( FM2.dat.obj )
#--------------------------------------------------------------
report.FM2.dat<- function( dat_obj )
{
	curve <- dat_obj$curve;
	if (is.null(curve))
	{
		stop("inavlid curve type");
	}

	rprlist<- dat.report( curve, dat_obj );
	return(rprlist);
}


#--------------------------------------------------------------
# plot.FM2.dat
#
# used by plot( FM2.dat.obj )
#--------------------------------------------------------------
plot.FM2.dat<- function( dat_obj, plot_type=NA )
{
	if (is.null(FM2.curve))
	{
		stop("inavlid curve type");
	}

	dat.plot( dat_obj, plot_type );
	invisible();
}

#--------------------------------------------------------------
# FM2.hp_test
#
# Hypothesis test for data object, the test methods should be 
# (10,11,12,13,14,15)
#
# 1) scan_step, default=2, an interval distance used to scan flanking 
#    marker, default is 2cm.
# 2) peak_count, default=5, a number shows how many significant QTL will 
#    be selected.
# 3) plot_doctype, default=pdf, the figure output type for summary command.
#
#--------------------------------------------------------------
FM2.qtlmodel<-function( dat_obj, options=list() )
{
	#scan_step
	set.scan_step <- FALSE;
	if (!(is.null(options$scan_step) || is.na(options$scan_step) ))
	{	
		set.scan_step 	<- TRUE;
		old.scan_step   <- FM2.set_value("scan_step", options$scan_step);
	}

	#peak_count
	set.peak_count 	<- FALSE;
	if (!(is.null(options$peak_count) || is.na(options$peak_count) ) )
	{	
		set.peak_count 	<- TRUE;
		old.peak_count  <- FM2.set_value("peak_count", options$peak_count);
	}

	if (is.null(FM2.curve))
	{
		stop("inavlid curve type");
	}

	ret<- FM2.model$qtlscan( dat_obj );

	if (set.scan_step)
		FM2.set_value("scan_step", old.scan_step);
	if (set.peak_count)
		FM2.set_value("peak_count", old.peak_count);

	return (ret);
}



#--------------------------------------------------------------
# FM2.permutation
#
# Permutation tests.
#--------------------------------------------------------------
FM2.permutation<-function( dat_obj, options=list() )
{
	old.debug <- NULL;
	set.debug <- FALSE;
	if (!is.null(options$debug) && !is.na(options$debug) )
	{
        	old.debug <- FM2.set_value("debug", as.numeric( options$debug));
		set.debug <- TRUE;
	}

	old.cluster_count <- NULL;
	set.cluster_count <- FALSE;
	if (!is.null(options$cluster_count) && !is.na(options$cluster_count) )
	{
        	old.cluster_count <- FM2.set_value("cluster_count", as.numeric( options$cluster_count));
		set.cluster_count <- TRUE;
	}

	old.permu_loop <- NULL;
	set.permu_loop <- FALSE;
	if (!is.null(options$permu_loop) && !is.na(options$permu_loop))
	{
        	old.permu_loop <- FM2.set_value("permu_loop", as.numeric( options$permu_loop) );
		set.permu_loop <- TRUE;
	}

	loops <- as.numeric( FM2.get_value("permu_loop", def=1000) );
	if (loops<100)
	{
		#stop("Permutation loop is too few(<100).\n");
		warning("Permutation loop is too few(<100).\n");
	}

	curve <- dat_obj$curve;
	if (is.null(FM2.curve))
	{
		stop("inavlid curve type");
	}

	ret<- permu.execute( dat_obj );
	class( ret ) <- "FM2.ret.perm";

	if (set.debug)
		FM2.set_value("debug", old.debug);
	if (set.cluster_count)
		FM2.set_value("cluster_count", old.cluster_count);
	if (set.permu_loop)
		FM2.set_value("permu_loop", old.permu_loop);

	return( ret );
}

#--------------------------------------------------------------
# summary.FM2.ret.perm
#
# used by summary( FM.ret.perm.obj )
#--------------------------------------------------------------
summary.FM2.ret.perm<-function( perm_obj, file=NA, append=TRUE )
{
	str <- permu.summary( perm_obj );

	if (is.na(file))
		return ( cat( str ) )
	else
		return ( cat( str, file= file, append = append) )
}

#--------------------------------------------------------------
# plot.FM2.ret.perm
#
# used by plot( FM.ret.perm object )
#--------------------------------------------------------------
plot.FM2.ret.perm<-function( perm_obj )
{
	permu.plot( perm_obj );
}

#--------------------------------------------------------------
# FM2.evaluate
#
# Evaluate the model parameters.
#--------------------------------------------------------------
FM2.evaluate<-function(par_obj, loops=100)
{
	ret <- eval.execute(par_obj, loops);
	return (ret );
}

#--------------------------------------------------------------
# summary.FM2.ret.eval
#
# used by summary( FM2.ret.eval.obj )
#--------------------------------------------------------------
summary.FM2.eval.ret<-function( eval_obj, file=NA, append=TRUE )
{
	str <- eval.summary( eval_obj );

	if (is.na(file))
		return ( cat( str ) )
	else
		return ( cat( str, file= file, append = append) )
}

#--------------------------------------------------------------
# plot.FM2.ret.eval
#
# used by plot( FM.ret.eval object )
#--------------------------------------------------------------
plot.FM2.eval.ret<-function( eval_obj )
{
	eval.plot( eval_obj );
}


#--------------------------------------------------------------
# public: FM2.report()
#               
#--------------------------------------------------------------

FM2.report<-function( pdf.file, dat, res, perm=NULL, options=list() )
{
	Report.new( pdf.file, options );
	Report.title( "Functional Mappping Report", "FunMap Report", "http://statgen.psu.edu" );
	Report.par( "dat.file", dat$pheno_file);
	
	Report.AddHeadline( "Data Summary", level=1 );
	fre.report_dat(dat);

	Report.AddHeadline( "QTL Profile", level=1 );
	
	#output the result of curve(1);
	if ( !is.null(res) )   
		fre.report_res(dat, res, perm) ;

	Report.AddHeadline( "QTL Position", level=1 );
	
	#output the result of curve(2);
	if ( !is.null(res) )	  
		fre.report_res2(dat, res, perm) ;

	if ( !is.null(perm) )	  
	{
		Report.AddHeadline( "Permutation Results", level=1 );
		fre.report_perm(dat, perm) ;
	}		

	Report.Output( pdf.file );
}


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Function: FM2.one_step
#
# The result is saved into [pheno_csv.rdata]
#
# Options:
#1) cluster_count, default=1, the cluster count for parallel permutation.
#2) permu_loop, default=1000, the count of permutation loop.
#3) file.summary, default=NA, a filename where the summary information 
#   which displayed in the console will be saved.
#4) file.rdata, default=NA, a filename for RDATA format file where the 
#   data object, result object and permutation object will be saved.
#5) file.report, default=NA, a filename for PDF report file where the 
#   summaries and figures will be saved.
#6) scan_step, default=2, an interval distance used to scan flanking 
#   marker, default is 2cm.
#7) peak_count, default=5, a number shows how many significant QTL will 
#   be selected.
#8) plot_doctype, default=pdf, the figure output type for summary command.
#9) np.order, default=6, the order of Legendre polynomial for nonparametric 
#   method.
#10) sg.order, default=3, the order of segmental lines.
#11) CM.model1, vector, default=c(CURVE_LC, 3), the 1st model(curve) type 
#    and count of parameters for composite model.
#12) CM.model2, vector, default=c(CURVE_NP, 4) ), the 2nd model(curve) type 
#    and count of parameters for composite model.
#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
FM2.qtlscan<-function(pheno_csv, geno_csv, marker_csv, curve_type, cross_type, covar_type=COVAR_AR1, options=list())
{
	if (is.null( FM_sys ) )
		FM2.start();

	#load data
	dat_os<- FM2.load_data( pheno_csv, geno_csv, marker_csv, curve_type, cross_type, covar_type );
	if (is.null(dat_os))
		stop("STOP:Failed to load the phenotype, genotype and marker file.\n");
	
	file.summary<-"";
	if (!is.null(options$file.summary))
		file.summary <- options$file.summary;
	cat("\nSummary report:\n", file=file.summary);
	
	try( summary( dat_os, file=file.summary ) );

	ret_os<- FM2.qtlmodel( dat_os, options=options );

	##Save 1
	if (is.null( options$file.rdata) )
		options$file.rdata = paste(pheno_csv, ".rdata", sep="")
	save(dat_os, ret_os, file=options$file.rdata );

	#report
	if (is.null( options$file.report))
		options$file.report <- paste(pheno_csv, ".pdf", sep="");
	FM2.report( options$file.report, dat_os, ret_os );
	
	nLoop <- FM2.get_value("permu_loop", def=1000);
	if ( !is.null(options$permu_loop) && !is.na(options$permu_loop) )
		nLoop <- as.numeric( options$permu_loop );

	ret_perm_os <- NULL;
	if (nLoop>0)
	{
		ret_perm_os<- FM2.permutation( dat_os, options=options);
		summary( ret_perm_os, file=file.summary );

		p05 <- 0;
		p05 <- 0;
		idx05 <- which(ret_perm_os$pv_table[,1]==0.05);
		idx01 <- which(ret_perm_os$pv_table[,1]==0.01);
		if (length(idx05)>0 )
			p05 <- ret_perm_os$pv_table[idx05[1],2];
		if (length(idx01)>0 )
			p01 <- ret_perm_os$pv_table[idx01[1],2];

		if (p05>0 && p01>0)
			ret_os <- FM2.model$set_lr2_cutoff( ret_os, p05, p01);

		try( summary(ret_os, dat_os, file=file.summary) );

		##Save 2;
		if (is.null( options$file.rdata) )
			save(dat_os, ret_os, ret_perm_os, file=paste(pheno_csv, ".rdata", sep="") )
		else
			save(dat_os, ret_os, ret_perm_os, file=options$file.rdata );
	}
	else
		try( summary(ret_os, dat_os, file=file.summary) );
	
	#report
	if (is.null( options$file.report))
		options$file.report <- paste(pheno_csv, ".pdf", sep="");
	FM2.report( options$file.report, dat_os, ret_os, ret_perm_os )

	return(ret_os);
}


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# FM2.simu_test
#
# Abstract: any model, any cross by simulation data
# 
# The result is saved into [LC_simu_test_XX_XX.rdata]
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
FM2.simu_test<-function( par_obj, curve_type, cross_type, covar_type= COVAR_AR1 )
{
	if (is.null( FM_sys ) )
		FM2.start();

	par_st<- FM2.param( par_obj, curve_type, cross_type, covar_type );
	summary(par_st);
	
	dat_st<- FM2.simulate( par_st);
	summary(dat_st);

	ret_st<- FM2.qtlmodel( dat_st);
	summary(ret_st, dat_st);

	##Step 1
	rdata.name <- paste( "simu_test_",FM2.curve$name,"_",FM2.cross$name,".rdata", sep="");
	save(dat_st, ret_st, file = rdata.name );

	##Step 2(Report) 
	pdf.name <- paste( "simu_test_",FM2.curve$name,"_",FM2.cross$name,".pdf", sep="");
	FM2.report(pdf.name, dat_st, ret_st);

	##Step 3(Ppermutation) 
	FM2.set_value("debug", TRUE)
	ret_perm_st<- FM2.permutation( dat_st, options=list(permu_loop=100) );
	summary( ret_perm_st );
	try(plot(ret_perm_st));

	FM2.report(pdf.name, dat_st, ret_st, ret_perm_st);
	
	##Step 4(Save All Data)
	save(dat_st, ret_st, ret_perm_st, file = rdata.name );
	cat("\nNotice: The results(dat_st, ret_st and ret_perm_st) are stored into ", rdata.name, "\n");

	return( list(dat_st, ret_st, ret_perm_st) );
}
