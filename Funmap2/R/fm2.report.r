###############################################################
# 
# New Systems Mapping Application(SysMap1)
#
# report utility
#    1). fre.report_dat   
#    2). fre.report_res   
#    3). fre.report_res2   
#
# History:
# 12/15/2011 Version 1.1
#
###############################################################

fre.report_dat<-function( dat ) 
{
	if ( class(dat) != "FM2.dat" )
	{
		sErrMsg<- "Error: Not a data set for Functional Mapping.";
		stop( sErrMsg );
	}

	Report.AddParaLine( tab.pos=c(300, 340, -1 ), tab.align=c("R", "C", "L"), c("Curve:", "", dat$curve_name) );
	Report.AddParaLine( tab.pos=c(300, 340, -1 ), tab.align=c("R", "C", "L"), c("Cross:", "", FM2.cross$name ) );
	Report.AddParaLine( tab.pos=c(300, 340, -1 ), tab.align=c("R", "C", "L"), c("Covariance:", "", FM2.covar$name ) );
	Report.AddParaLine( tab.pos=c(300, 340, -1 ), tab.align=c("R", "C", "L"), c("Pheno. file:", "", dat$pheno_file ) );
	Report.AddParaLine( tab.pos=c(300, 340, -1 ), tab.align=c("R", "C", "L"), c("Geno. file:",  "", dat$geno_file ) );
	Report.AddParaLine( tab.pos=c(300, 340, -1 ), tab.align=c("R", "C", "L"), c("Marker file:", "", dat$marker_file) );
	Report.AddParaLine( tab.pos=c(300, 340, -1 ), tab.align=c("R", "C", "L"), c("Sample size:", "", dat$sample_N) );
	Report.AddParaLine( tab.pos=c(300, 340, -1 ), tab.align=c("R", "C", "L"), c("Sample times:","", length(dat$sample_times) ) );
	Report.AddParaLine( tab.pos=c(300, 340, -1 ), tab.align=c("R", "C", "L"), c("Marker count:","", dim(dat$marker_table)[1] ) );

	Report.AddParaBreak();

	c1 <- call("fpt.plot_tiled_curves", dat );
	Report.AddFigure( c1, " ", c(5, 5)*254, left.margin=0.1*254);

	c2 <- call("fpt.plot_overlapping_curves", dat );
	Report.AddFigure( c2, " ", c(5, 5)*254, left.margin=0.1*254);

	Report.AddParaBreak();

	return();
} 


fre.report_res<-function( dat, res, perm=NULL ) 
{
	if ( is.null(res) )
	{
		sErrMsg<- "Error: Not a result for FunMap.";
		stop( sErrMsg );
	}

	par_QQ = c();
	par_Qq = c();
	par_qq = c();

	par.n <- length(res$LR_peaks[1,])-6;
	par.g <- 2;
	if (res$cross_type ==CROSS_BC || res$cross_type ==CROSS_RIL )
		par.n <- par.n/2
	else
	{
		par.n <- par.n/3;
		par.g <- 3;
	}

	for(i in 1:length(res$LR_peaks[,1]) )
	{
		par.str <- c();
		for(j in 0:(par.g-1) )
		{
			par.str0 <- "";
			for( k in 1:par.n)
			{
				par.str0 <- sprintf("%s%.2f,", par.str0, res$LR_peaks[i, 6 + j*par.n + k ] );
			}
			
			if (nchar(par.str0)>20)
				par.str0 = paste( strtrim(par.str0, 20), "...", sep="");

			par.str <- c(par.str, par.str0);
		}
	
		if (res$cross_type ==CROSS_BC)
		{
			par_QQ <- c( par_QQ, "");
			par_Qq <- c( par_Qq, par.str[1]);
			par_qq <- c( par_qq, par.str[2]);
		}								 

		if (res$cross_type ==CROSS_RIL)
		{
			par_QQ <- c( par_QQ, par.str[1]);
			par_Qq <- c( par_Qq, "");
			par_qq <- c( par_qq, par.str[2]);
		}
		
		if (res$cross_type ==CROSS_F2)
		{
			par_QQ <- c( par_QQ, par.str[1]);
			par_Qq <- c( par_Qq, par.str[2]);
			par_qq <- c( par_qq, par.str[3]);
		}

	}

	sigs_list <- data.frame( 
				Grp = c(res$LR_peaks[,1, drop=F]),
				Pos = sprintf("%.2f", res$LR_peaks[,2, drop=F] ),
				LR  = sprintf("%.2f", res$LR_peaks[,3, drop=F] ),			
				rho = sprintf("%.2f", res$LR_peaks[,5, drop=F] ), 
				s2  = sprintf("%.2f", res$LR_peaks[,6, drop=F] ),
				QQ  = par_QQ,
				Qq  = par_Qq,
				qq  = par_qq);
		


	Report.AddTable( sigs_list, 
			 	title = "LR profile",
				frame = "wang",  
				col.names = c("Grp.", "Pos.", "LR", "rho", "Sig^2", "QQ", "Qq", "qq" ),
				col.width = c(100,    100,    200,   200,  120,  400, 400, 400  ),
				col.align = c("L",    "L",    "R",   "R",  "R",  "L", "L", "L"),
				offset.x = 50, 
				max.show = 20)

	Report.AddParaBreak(50);
	
	
	perm.cutoff.05 <- fin.find_cutoff(perm, 0.05);
	perm.cutoff.01 <- fin.find_cutoff(perm, 0.01);
	
	c1 <- call("fpt.plot_qtl_map", dat, res$full_res, perm.cutoff.05, perm.cutoff.01 );
	Report.AddFigure( c1, "QTL profile map", c(6, 6)*254, left.margin=0.5*254);
	Report.AddParaBreak();

	return();
} 

fin.find_cutoff <- function(perm, p=0.05)
{
	perm.cutoff <- NULL;
	if (!is.null(perm))
	{
		if( class(perm)=="FM2.ret.perm" )
		{
			i.pos <- which(perm$pv_table[,1]==p);
			if( length(i.pos)>0 ) perm.cutoff <- perm$pv_table[ i.pos[1], 2 ]
		}
		else
		{
			if( class(perm)=="list" )
			{
				if (p==0.05) perm.cutoff <- perm[[1]];
				if (p==0.01) perm.cutoff <- perm[[2]];
			}	
		}			
	}
	
	return( perm.cutoff );
}

fre.report_res2<-function( dat, res, perm=NULL ) 
{
	f_curve_mu <- FM2.curve$get_mu;
	if(is.null(f_curve_mu))
	{
		curve <- FM2.get_curve(dat$curve_type);
		if (is.null(curve) || curve$name != dat$curve_name )
		{
			stop("inavlid curve type");
		}
		
		FM2.curve <<- curve;
		f_curve_mu <- FM2.curve$get_mu;
	}

	perm.cutoff.05 <- fin.find_cutoff(perm, 0.05);
	perm.cutoff.01 <- fin.find_cutoff(perm, 0.01);

	nMesa   <- max(dat$sample_times);
	nLong   <- (nMesa*6)%/%5;

	simu_QQ <- NULL;
	simu_Qq <- NULL;
	simu_qq <- NULL;

	par.n <- length(res$LR_peaks[1,])-6;
	par.g <- 2;
	if (res$cross_type ==CROSS_BC || res$cross_type ==CROSS_RIL )
		par.n <- par.n/2
	else
	{
		par.n <- par.n/3;
		par.g <- 3;
	}

	for(i in 1:length(res$LR_peaks[,1]) )
	{
		str <- sprintf("Groupp=%d, Postion=%.2f", res$LR_peaks[i,1], res$LR_peaks[i,2] );
		Report.AddHeadline( str, level=2 );

		c1 <- call("fpt.plot_qtl_pos", res$LR_peaks[i,1], dat, res$full_res, cutoff.05=perm.cutoff.05, cutoff.01=perm.cutoff.01, qtl_ps=res$LR_peaks[i,2] );
		Report.AddFigure( c1, "", c(3.5, 3.5)*254, left.margin=0*254);

		QQ_par  <- NULL;
		Qq_par  <- NULL;
		qq_par  <- NULL;

		#get QQ_par, Qq_par, qq_par
		par.str <- c();
		for(j in 0:(par.g-1) )
		{
			par.str0 <- c();
			for( k in 1:par.n)
			{
				par.str0 <- c(par.str0, res$LR_peaks[i, 6 + j*par.n + k ] );
			}
			
			par.str <- rbind(par.str, par.str0);
		}
	
		if (res$cross_type ==CROSS_BC)
		{
			Qq_par <- c( Qq_par, par.str[1,]);
			qq_par <- c( qq_par, par.str[2,]);
		}

		if (res$cross_type ==CROSS_RIL)
		{
			QQ_par <- c( QQ_par, par.str[1,]);
			qq_par <- c( qq_par, par.str[2,]);
		}
		
		if (res$cross_type ==CROSS_F2)
		{
			QQ_par <- c( QQ_par, par.str[1,]);
			Qq_par <- c( Qq_par, par.str[2,]);
			qq_par <- c( qq_par, par.str[3,]);
		}

		c2 <- call("fpt.plot_com_curve", nMesa, nLong, f_curve_mu, dat, 
				QQ_par=QQ_par, Qq_par=Qq_par, qq_par=qq_par, 
				simu_QQ=simu_QQ, simu_Qq=simu_Qq, simu_qq=simu_qq, 
				xlab="Time", ylab="Model");

		Report.AddFigure( c2, "", c(3.5, 3.5)*254, left.margin=0*254);

		Report.AddParaBreak();
	}

	return();
} 

fre.report_perm<-function(dat, perm)
{
	if ( is.null(perm) )
	{
		sErrMsg<- "Error: Not a permutation result for FunMap.";
		stop( sErrMsg );
	}

	p.pv <- c( 0.1, 0.05, 0.01, 0.005, 0.001, 0.0001, 0.00001);
	p.cut <- perm$pv_table[ match(p.pv, perm$pv_table[,1]), 2];
	
	p.df <- data.frame(N0="Cutoff", N1=p.cut[1], N2=p.cut[2], N3=p.cut[3], N4=p.cut[4], N5=p.cut[5], N6=p.cut[6], N7=p.cut[7])

	Report.AddTable( p.df, 
		 	title = "Permutation Cutoff",
			frame = "wang",  
			col.names = c("pv", "0.1", "0.05", "0.01", "0.005", "0.001", "0.0001", "0.00001" ),
			col.width = c(100,   150,   150,    150,    150,      150,      150,      150  ),
			col.align = c("L",   "R",   "R",    "R",    "R",      "R",       "R",     "R"),
			offset.x = 50, 
			max.show = 20)

	Report.AddParaBreak(50);
	
	c2 <- call("fpt.plot_permutation", perm$pv_table);

	Report.AddFigure( c2, "", c(3.5, 3.5)*254, left.margin=0*254);

	Report.AddParaBreak();
	
	return();
}