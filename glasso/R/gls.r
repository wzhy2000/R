#config param
#
# debug=1
# max_iter=2000
# rho_tuning=0.09
# burn_in_round=0.2
# fQVal=0.05
#simu_a_effect[1] = 1, 1.04, 0.885, -2.055, 0.545
#simu_a_effect[2] = 2, 1.17, -0.20, 0.74, -4.715
#simu_a_effect[3] = 3, 1.40, -2.25, 1.00,  0.00

#simu_d_effect[1] = 3, 1.49, -2.135, 4.82, 1.425
#simu_d_effect[2] = 4, 1.045, 1.320, 1.905,  1.535
#simu_d_effect[3] = 5, 1.265, -1.225, 2.710, -1.96

gls.simulate<-function( phe.out, snp.out, simu_grp=1, simu_n= 200, simu_p, simu_snp_rho, simu_rho, simu_sigma2= 16, 
		 simu_mu= c(13.395, -3.08, 1.875, -3.195), 
		 simu_covar_effect = array(c(0,0,0,0), dim=c(1,4)), simu_covar_range = c(-1, 1),
		 simu_add_effect=NA,  simu_dom_effect=NA, 
		 simu_z_range = c(20,80), simu_z_count = c(5,12), debug=F)
{
	if( is.na(simu_grp) || length(simu_grp) > 1)
		stop("The parameter of simu_grp is not a single valid value.");

	if( is.na(simu_n) || length(simu_n) > 1)
		stop("The parameter of simu_n is not a single valid value.");

	if( is.na(simu_p) || length(simu_p) > 1)
		stop("The parameter of simu_n is not a single valid value.");
	
	if( is.na(simu_snp_rho) || length(simu_snp_rho) > 1)
		stop("The parameter of simu_snp_rho is not a single valid value.");

	if( is.na(simu_rho) || length(simu_rho) > 1)
		stop("The parameter of simu_rho is not a single valid value.");
	
	if( is.na(simu_sigma2) || length(simu_sigma2) > 1)
		stop("The parameter of simu_sigma2 is not a single valid value.");

	if( length(which(is.na(simu_mu)))>0 || length(simu_mu) != 4 )
		stop("The parameter of simu_mu is not a vector with 4 numeric valid values.");

	if( length(which(is.na(simu_covar_range)))>0 || length(simu_covar_range) != 2 )
		stop("The parameter of simu_cov_range is not a range.");

	if( length(which(is.na(simu_z_range)))>0 || length(simu_z_range) != 2 )
		stop("The parameter of simu_z_range is not a range.");

	if( length(which(is.na(simu_covar_range)))>0 || length(simu_covar_range) != 2 )
		stop("The parameter of simu_covar_range is not a range.");

	if(!is.matrix(simu_covar_effect))
		stop("The parameter of simu_covar_effect is not a matrix.");
	if(NCOL(simu_covar_effect)!=4)
		stop("The parameter of simu_covar_effect is not a matrix with 4 columns.");
	if(length(which(is.na(simu_covar_effect)))>0 )
		stop("The parameter of simu_covar_effect has NA values.");
	
	if(!is.matrix(simu_add_effect))
		stop("The parameter of simu_add_effect is not a matrix.");
	if(NCOL(simu_add_effect)!=5)
		stop("The parameter of simu_add_effect is not a matrix with 5 columns.");
	if(length(which(is.na(simu_add_effect)))>0 )
		stop("The parameter of simu_add_effect has NA values.");

	if(!is.matrix(simu_dom_effect))
		stop("The parameter of simu_dom_effect is not a matrix.");
	if(NCOL(simu_dom_effect)!=5)
		stop("The parameter of simu_dom_effect is not a matrix with 5 columns.");
	if(length(which(is.na(simu_dom_effect)))>0 )
		stop("The parameter of simu_dom_effect has NA values.");

	simu_sig_add <- NROW(simu_add_effect);
	simu_sig_dom <- NROW(simu_dom_effect);
	
	sigp<-unique(c(simu_add_effect[,1], simu_dom_effect[,1]))
	simu_sigp <- length(sigp);
	err <- 0;
	
	out <- .C("gls_simulate", 
		   as.character(phe.out),			# char* szPhe_out
  		   as.character(snp.out), 			# char* szSnp_out
  		   as.integer(simu_grp), 			# int nSimu_grp
  		   as.integer(simu_n), 				# int nSimu_n
  		   as.integer(simu_p), 				# int nSimu_p
  		   as.double(simu_snp_rho), 			# double fSimu_snp_rho
  		   as.double(simu_rho), 			# double fSimu_rho
  		   as.double(simu_sigma2), 			# double fSimu_sigma2
		   as.double(as.vector(simu_mu)),		# double* pfSimu_mu
		   as.integer(NROW(simu_covar_effect)), 	# int nSimu_covar_len
		   as.double(as.vector(simu_covar_range)),	# double* pfSimu_covar_range
		   as.double(as.matrix(simu_covar_effect)), 	# double* pfSimu_covar_effect
		   as.integer(simu_sigp),			# int nSimu_sig_p
		   as.integer(simu_sig_add),			# int nSimu_add_len
		   as.double(as.matrix(simu_add_effect)), 	# double* pfSimu_add_effect
		   as.integer(simu_sig_dom),			# int nSimu_dom_len
		   as.double(as.matrix(simu_dom_effect)), 	# double* pfSimu_dom_effect
		   as.double(as.vector(simu_z_range)),		# double* simu_z_range
		   as.integer(as.vector(simu_z_count)), 	# int* pnSimu_z_count
		   as.integer(debug),
		   as.integer(err) );
		   
	return(err);		   
}

gls.simple<-function(phe.file, snp.file, model, bRefit=T,    
	    nMaxIter = 2000,
	    fBurnInRound = 0.3,
	    fRhoTuning = 0.09,
	    fQval.add  = 0.05,
	    fQval.dom  = 0.09,
	    debug=F)
{
	r <- .Call("gls_simple", 
			phe.file,
  		   	snp.file, 
  		   	model, 
  		   	bRefit,
  		   	nMaxIter,
		   	fBurnInRound,
		   	fRhoTuning,
	           	fQval.add,
	           	fQval.dom,
			debug);
	add_c19 <- c("grp", "pos", "add.m1", "add.m2", "add.m3", "add.m4",
		"add_L2.mu", "add.mu1", "add.mu2", "add.mu3", "add.mu4",
		"add_L2.min", "add.min1", "add.min2", "add.min3", "add.min4",
		"add_L2.max", "add.max1", "add.max2", "add.max3", "add.max4" );
	dom_c19 <- c("grp", "pos", "dom.m1", "dom.m2", "dom.m3", "dom.m4",
		"dom_L2.mu", "add.mu1", "add.mu2", "add.mu3", "dom.mu4",
		"dom_L2.min", "dom.min1", "dom.min2", "dom.min3", "dom.min4",
		"dom_L2.max", "dom.max1", "dom.max2", "dom.max3", "dom.max4" );
	cov_c19 <- c("cov.m1", "cov.m2", "cov.m3", "cov.m4",
		"co_L2.mu", "cov.mu1", "cov.mu2", "cov.mu3", "cov.mu4",
		"co_L2.min", "cov.min1", "cov.min2", "cov.min3", "cov.min4",
		"co_L2.max", "cov.max1", "cov.max2", "cov.max3", "cov.max4" );


	if(!is.null(r) && !is.na(r))
	{
		if (!is.null(r$varsel_add)) colnames(r$varsel_add) <- add_c19;
		if (!is.null(r$varsel_dom)) colnames(r$varsel_dom) <- dom_c19;
		if (!is.null(r$varsel_cov)) colnames(r$varsel_cov) <- cov_c19;
		if (!is.null(r$refit_add)) colnames(r$refit_add) <- add_c19;
		if (!is.null(r$refit_dom)) colnames(r$refit_dom) <- dom_c19;
		if (!is.null(r$refit_cov)) colnames(r$refit_cov) <- cov_c19;
		
		row.cov <- c<-("Mu");
		if(NROW(r$varsel_cov)>1)
		{
			for(k in 1:(NROW(r$varsel_cov)-1))
				row.cov <- c(row.cov, paste("Cov_",k,sep=""));
			rownames(r$varsel_cov) <- row.cov;
		}
		
		row.cov <- c<-("Mu");
		if(NROW(r$refit_cov)>1)
		{
			for(k in 1:(NROW(r$refit_cov)-1))
				row.cov <- c(row.cov, paste("Cov_",k,sep=""));
			rownames(r$refit_cov) <- row.cov;
		}
		
		class(r) <- "GLS.ret";
	}

	return(r);		   
}

gls.plink<-function( phe.file, tped.file, tfam.file, model, bRefit=TRUE,    
	    nMaxIter   = 2000,
	    fBurnInRound = 0.3,
	    fRhoTuning = 0.095,
	    fQval.add  = 0.05,
	    fQval.dom  = 0.09,
	    debug=F)
{
	r <- .Call("gls_plink", 
			phe.file,
  		   	tped.file, 
  		   	fam.file, 
  		   	model, 
  		   	bRefit,
  		   	nMaxIter,
		   	fBurnInRound,
		   	fRhoTuning,
	           	fQval.add,
	           	fQval.dom,
			debug);

	add_c19 <- c("grp", "pos", "add.m1", "add.m2", "add.m3", "add.m4",
		"add_L2.mu", "add.mu1", "add.mu2", "add.mu3", "add.mu4",
		"add_L2.min", "add.min1", "add.min2", "add.min3", "add.min4",
		"add_L2.max", "add.max1", "add.max2", "add.max3", "add.max4" );
	dom_c19 <- c("grp", "pos", "dom.m1", "dom.m2", "dom.m3", "dom.m4",
		"dom_L2.mu", "add.mu1", "add.mu2", "add.mu3", "dom.mu4",
		"dom_L2.min", "dom.min1", "dom.min2", "dom.min3", "dom.min4",
		"dom_L2.max", "dom.max1", "dom.max2", "dom.max3", "dom.max4" );
	cov_c19 <- c("cov.m1", "cov.m2", "cov.m3", "cov.m4",
		"co_L2.mu", "cov.mu1", "cov.mu2", "cov.mu3", "cov.mu4",
		"co_L2.min", "cov.min1", "cov.min2", "cov.min3", "cov.min4",
		"co_L2.max", "cov.max1", "cov.max2", "cov.max3", "cov.max4" );

	if(!is.null(r) && !is.na(r))
	{
		if (!is.null(r$varsel_add)) colnames(r$varsel_add) <- add_c19;
		if (!is.null(r$varsel_dom)) colnames(r$varsel_dom) <- dom_c19;
		if (!is.null(r$varsel_cov)) colnames(r$varsel_cov) <- cov_c19;
		if (!is.null(r$refit_add)) colnames(r$refit_add) <- add_c19;
		if (!is.null(r$refit_dom)) colnames(r$refit_dom) <- dom_c19;
		if (!is.null(r$refit_cov)) colnames(r$refit_cov) <- cov_c19;

		row.cov <- c<-("Mu");
		for(k in 1:NROW(r$varsel_cov))
			row.cov <- c(row.cov, paste("Cov_",k,sep=""));
		rownames(r$varsel_cov) <- row.cov;

		row.cov <- c<-("Mu");
		for(k in 1:NROW(r$refit_cov))
			row.cov <- c(row.cov, paste("Cov_",k,sep=""));
		rownames(r$refit_cov) <- row.cov;
		
		class(r) <- "GLS.ret";
	}
			
	return(r);		   
}

summary_output2<-function(re1, re_add, re_dom)
{
	cat("(1) Covariate Estimate:\n");
		
	for(i in 1:NROW(re1))
		cat(sprintf("%s \t %s \t %d%d%d%d \t%.3f\t(%.3f,%.3f,%.3f,%.3f)\n", 
			rownames(re1)[i], ifelse(re1[i,1], "Yes", "---"), 
			re1[i,1], re1[i,2], re1[i,3], re1[i,4], 
			re1[i,5], re1[i,6], re1[i,7], re1[i,8], re1[i,9] ));  

	cat("(2) Significant SNPs Estimate:\n");

	cat("    SNP Name\tGrp/Pos\tAdd\tLR2(Median1,2,3,4)\t\tDom\tLR2(Median1,2,3,4)\n");
	for(i in 1:NROW(re_add))
	{
		cat(sprintf("%12s\t%d/%d\t", rownames(re_add)[i], re_add[i,1], re_add[i,2] ));
  	        cat(sprintf("%s(%d%d%d%d)\t%.3f(%.3f,%.3f,%.3f,%.3f)", ifelse(sum(re_add[i,3:6])>0, "Yes/", "---/"),
		    re_add[i,3], re_add[i,4], re_add[i,5], re_add[i,6], 
		    re_add[i,7], re_add[i,8], re_add[i,9], re_add[i,10], re_add[i,11])) 

  	        cat(sprintf("%s(%d%d%d%d)\t%.3f(%.3f,%.3f,%.3f,%.3f)", ifelse(sum(re_dom[i,3:6])>0, "Yes/", "---/"),
		    re_dom[i,3], re_dom[i,4], re_dom[i,5], re_dom[i,6], 
		    re_dom[i,7], re_dom[i,8], re_dom[i,9], re_dom[i,10], re_dom[i,11])) 
		
		cat("\n");  
	}	
}

summary.GLS.ret<-function(r.gls)
{
	if(!is.null(r.gls$refit_add))
	{
		cat("Refit Result\n");
		summary_output2(r.gls$refit_cov, r.gls$refit_add, r.gls$refit_dom)
	}
	else if(!is.null(r.gls$varsel_add))
	{
		cat("Variable Selection Result\n");
		summary_output2(r.gls$varsel_cov, r.gls$varsel_add, r.gls$varsel_dom); 
	}
}

plot.GLS.ret<-function(r.gls, bRefit=T)
{
	ad <- NA;
	if (bRefit)
	{
		if( !is.null(r.gls$refit_add) )
		{
			ad <- rbind( r.gls$refit_add[,7], r.gls$refit_dom[,7]);
		}
		else
			stop("No refit results\n");		
	}
	else
	{
		if( !is.null(r.gls$varsel_add) )
		{
			ad <- rbind( r.gls$varsel_add[,7], r.gls$varsel_dom[,7]);
		}
		else
			stop("No varible selection results\n");		
	}
	
	par(mar=c(4.5, 4, 0.5, 2) + 0.1);
	plot.new();
  	par(mfrow=c(2,1));
	
	par(mfg=c(1, 1));
	draw_snplist( ad[,1], "Estimated additive effect", sigpos_list=NULL);
	par(mfg=c(2, 1));
	draw_snplist( ad[,2], "Estimated dominant effect", sigpos_list=NULL);
}


gls.outputpdf<-function( r.gls, fig.file, bRefit=T)
{
	ad <- NA;
	if (bRefit)
	{
		if( !is.null(r.gls$refit_add) )
		{
			ad <- rbind( r.gls$refit_add[,7], r.gls$refit_dom[,7]);
		}
		else
			stop("No refit results\n");		
	}
	else
	{
		if( !is.null(r.gls$varsel_add) )
		{
			ad <- rbind( r.gls$varsel_add[,7], r.gls$varsel_dom[,7]);
		}
		else
			stop("No varible selection results\n");		
	}
	
	pdf(fig.file, width=6, height=4.5);

	par(mar=c(4.5, 4, 0.5, 2) + 0.1);
	plot.new();
  	par(mfrow=c(2,1));
	
	par(mfg=c(1, 1));
	draw_snplist( ad[,1], "Estimated additive effect", sigpos_list=NULL);
	par(mfg=c(2, 1));
	draw_snplist( ad[,2], "Estimated dominant effect", sigpos_list=NULL);
	
	dev.off();
}

draw_snplist<-function( values, yLabel, sigpos_list=NULL)
{
	xlim <- c(0, length(values));
	ylim <- range(values);
	nB <- length(values);
	snps <- c(0:nB);

	plot(1:10,1:10, xlim=xlim, ylim=ylim, type="n", xlab="SNP", ylab=yLabel);
	rect(snps[-(nB+1)], 0, snps[-1L], values,  col = "blue", border = "black");
}

#--------------------------------------------------------------
# plot_sig_curve
# 
# Input:pvs[,1]  snp_name
#       pvs[,2]  chromoseom no
#       pvs[,3]  position
#       pvs[,4-7]    QQ or Additive
#       pvs[,8-11]   Qq or Dominant
#       pvs[,12-15]  qq
# fig.file: pdf file name
#
# Used by: GLS,DYN
#--------------------------------------------------------------
plot_sig_curve<-function( pvs.csv, fig.file, n.lgr = 4)
{
	pvs <- read.csv( pvs.csv, header=TRUE);

	n.row <- (length(pvs[,2])-1)%/%3 + 1;
	n.col <- 3;

	pdf(fig.file, width=6, height=n.row*2);
	
  	par( mfrow=c(n.row,n.col) );

	for (i in 1:n.row )
	for (j in 1:n.col )
	{
		if ( (i-1)*3+j > length(pvs[,2]) )
			break;

		n.par <- (i-1)*3+j;
		par(mfg=c(i, j));
		QQ.par<- pvs[n.par ,4:(4+n.lgr-1)];
		Qq.par<- pvs[n.par ,(4+n.lgr):(4+n.lgr*2-1)];
		qq.par <- NULL;
		if ( length(pvs[1,])>=(4+n.lgr*2) )
			qq.par<- pvs[n.par ,(4+n.lgr*2):(4+n.lgr*3-1)];

		if ( length(pvs[1,])>=(4+n.lgr*2) )
			draw_single_curve( pvs[n.par, 1], QQ=QQ.par, Qq=Qq.par, qq=qq.par )
		else
			draw_single_curve( pvs[n.par, 1], add=QQ.par, dom=Qq.par);
	}

	dev.off();
}

draw_single_curve<-function( snp_name, QQ=NULL, Qq=NULL, qq=NULL, add=NULL, dom=NULL )
{
	old.p1 <- par( mar=c(2,2,1,1)+0.1);
	on.exit(par(old.p1),add = T);

	tp <- seq(-1, 1, 0.05);
	#ui <- cbind( rep(1, length(tp)), tp, (3*tp^2-1)/2, (5*tp^3-3*tp)/2 ) ;
	ui <- cbind( rep(1, length(tp)), tp, (3*tp^2-1)/2, (5*tp^3-3*tp)/2, (35*tp^4-30*tp^2+3)/8 ) ;
	
	y <- c();
	if (!is.null(QQ))  y <- cbind(y, ui%*%t(QQ))
	if (!is.null(Qq))  y <- cbind(y, ui%*%t(Qq))
	if (!is.null(qq))  y <- cbind(y, ui%*%t(qq))
	if (!is.null(add)) y <- cbind(y, ui%*%t(add))
	if (!is.null(dom)) y <- cbind(y, ui%*%t(dom))

	ylim.max <- max(y, na.rm=T)*1.1;
	ylim.min <- min(y, na.rm=T)*1.1;
	if (ylim.min > 0)
	   ylim.min - min(y, na.rm=T)*0.9;
	if (ylim.min>0) ylim.min = 0;

	y.num <- length(tp);

	plot( c(0,0), c(0,0), type="n", xaxt="s", yaxt="s", yaxs="i", main=snp_name, 
		  xlab="Time", ylab="Y", xlim=c(-1, 1.2), ylim=c( ylim.min, ylim.max ) );

	cur.lab <- c();
	cur.col <- c();

	if (!is.null(QQ))
	{
		y <- ui%*%t(QQ);
		lines(tp, y, col="red");
		cur.lab <- c( cur.lab, "QQ");
		cur.col <- c( cur.col, "red");
	}

	if (!is.null(Qq))
	{
		y <- ui%*%t(Qq);
		lines(tp, y, col="green");
		cur.lab <- c( cur.lab, "Qq");
		cur.col <- c( cur.col, "green");
	}

	if (!is.null(qq))
	{
		y <- ui%*%t(qq);
		lines(tp, y, col="blue");
		cur.lab <- c( cur.lab, "qq");
		cur.col <- c( cur.col, "blue");
	}

	if (!is.null(add))
	{
		y <- ui%*%t(add);
		lines(tp, y, col="orange");
		cur.lab <- c( cur.lab, "Add");
		cur.col <- c( cur.col, "orange");
	}

	if (!is.null(dom))
	{
		y <- ui%*%t(dom);
		lines(tp, y, col="purple");
		cur.lab <- c( cur.lab, "Dom");
		cur.col <- c( cur.col, "purple");
	}

	legend( "topright", 
				   legend = cur.lab,
	               text.width = strwidth("ABC"),
				   text.col = cur.col,
	               col = cur.col,
				   lty=1,
	               xjust = 1, 
	               yjust = 1,
				   cex=0.8)
}
