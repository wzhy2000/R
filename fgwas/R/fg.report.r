fg.report<-function( fg.dat, fg.scan=NULL, fg.perm=NULL, file.pdf=NULL )
{
	check_FG_DAT( fg.dat );
	
	if( !missing(fg.scan) )	check_FG_SCAN( fg.scan );
	if( !missing(fg.perm) )	check_FG_PERM( fg.perm );
	if( !missing(file.pdf) ) check_file_pdf( file.pdf );

	if (is.null(file.pdf))
		file.pdf <- paste(fg.dat$file.pheY.csv, ".PDF", sep="");
		
	Report.new( file.pdf, options );
	Report.title( "Functional GWAS Report", "fGWAS", "http://ccb.bjfu.edu.cn/" );
	Report.par( "dat.file", fg.dat$file.pheY.csv);
	
	Report.AddHeadline( "Data Summary", level=1 );
	proc_report_dat(fg.dat);

	if( !is.null( fg.scan) ) 
	{
		Report.AddHeadline( "SNP Scan", level=1 );
		proc_report_snpscan(fg.dat, fg.scan, fg.perm) ;
	}

	if( !is.null( fg.perm) ) 
	{
		Report.AddHeadline( "SNP Permutation", level=1 );
		proc_report_perm(fg.dat, fg.perm) ;
	}
	
	if( !is.null( fg.perm) && !is.null( fg.scan) )
	{
		r.sig <- fg_detect_sig(fg.dat, fg.scan, fg.perm)

		Report.AddHeadline( "Significant SNP", level=1 );
		proc_report_sig( fg.dat, r.sig$sig.05, r.sig$sig.01 ) ;
	}

	if( is.null( fg.perm) && !is.null( fg.scan) )
	{
		Report.AddHeadline( "Top SNP", level=1 );
		proc_report_topn(fg.dat, fg.scan) ;
	}

	Report.Output( file.pdf );
}
