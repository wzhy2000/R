bls.simulate<-function( phe.out, snp.out, simu_grp=1, simu_n= 500, simu_p=1000, 
		simu_snp_rho = 0.1, 
		simu_rho     = 0.4, 
		simu_sigma2  = 3, 
		simu_mu      = 26, 
		simu_cov_coeff = c( 0, 2 ), 
		simu_a_pos   = c( 100, 200, 300), 
		simu_a_effect= c( 2.2, -2.5, 2.0 ),  
		simu_d_pos   = c( 300, 500, 700), 
		simu_d_effect= c( 2.8, 2.0, -2.5 ),
		simu_cov_range=c( 0, 1),
		simu_t_range = c(-1, 1), 
		debug=F )
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

	if( is.na(simu_mu) || length(simu_mu) > 1)
		stop("The parameter of simu_mu is not a single valid value.");

	if ( length(simu_a_pos) != length(simu_a_effect ) )
		stop("The length of simu_a_pos is same as simu_a_effect.");

	if ( length(simu_d_pos) != length(simu_d_effect ) )
		stop("The length of simu_d_pos is same as simu_d_effect.");

	if ( length(which(simu_a_pos<=0 | simu_a_pos>simu_p))>0  )
		stop("The parameter of simu_a_pos should be in correct SNP range.");

	if ( length(which(simu_a_pos<=0 | simu_a_pos>simu_p))>0  )
		stop("The parameter of simu_a_pos should be in correct SNP range.");

	if ( length(simu_a_pos)>0 && length(which(simu_d_pos<=0 | simu_d_pos>simu_p))>0  )
		stop("The parameter of simu_a_pos should be in correct SNP range.");

	if ( length(simu_a_effect)>0 && length(which(is.na(simu_a_effect)))>0  )
		stop("The parameter of simu_d_pos has NA values.");

	if ( length(simu_d_effect)>0 && length(which(is.na(simu_d_effect)))>0  )
		stop("The parameter of simu_d_pos has NA values.");

	if ( length(simu_cov_coeff)>0 && length(which(is.na(simu_cov_coeff)))>0  )
		stop("The parameter of simu_cov_coeff has NA values.");

	if ( length(simu_cov_range)>0 && length(which(is.na(simu_cov_range)))>0  )
		stop("The parameter of simu_cov_range has NA values.");

	if ( length(simu_t_range)>0 && length(which(is.na(simu_t_range)))>0  )
		stop("The parameter of simu_t_range has NA values.");
	
	if ( length(simu_t_range)!=2)
		stop("The parameter of simu_t_range should be a valid range.");

	if ( length(simu_cov_range)!=2)
		stop("The parameter of simu_cov_range should be a valid range.");

	if (simu_cov_range[1]>simu_cov_range[2])
		simu_cov_range<-c(simu_cov_range[2], simu_cov_range[1]);
		
	if (simu_t_range[1]>simu_t_range[2])
		simu_t_range<-c(simu_t_range[2], simu_t_range[1]);
	
	sigp<-unique(c(simu_a_pos, simu_d_pos))
	simu_sigp <- length(sigp);
	simu_a_len <- length(simu_a_pos);
	simu_d_len <- length(simu_d_pos);
	simu_covar_count <- length(simu_cov_coeff);

	err <- 0;
	out <- .C("bls_simulate", 
		   as.character(phe.out),		# char* szPhe_out
  		   as.character(snp.out), 		# char* szSnp_out
  		   as.integer(simu_grp), 		# int nSimu_grp
  		   as.integer(simu_n), 			# int nSimu_n
  		   as.integer(simu_p), 			# int nSimu_p
  		   as.double(simu_snp_rho), 		# double fSimu_snp_rho
  		   as.double(simu_rho), 		# double fSimu_rho
  		   as.double(simu_sigma2), 		# double fSimu_sigma2
		   as.double(simu_mu),			# double fSimu_mu
		   as.integer(simu_covar_count),
		   as.double(as.vector(simu_cov_coeff)),# double* pfSimu_cov_coeff
		   as.integer(simu_sigp),		# int nSimu_sig_p
		   as.integer(simu_a_len),
		   as.integer(as.vector(simu_a_pos)),	# int* nSimu_a_pos
		   as.double(as.vector(simu_a_effect)), # double* pfSimu_a_effect
		   as.integer(simu_d_len),
		   as.integer(as.vector(simu_d_pos)),	# int* nSimu_d_pos
		   as.double(as.matrix(simu_d_effect)), # double* pfSimu_d_effect
		   as.double(as.vector(simu_cov_range)),# double* pfSimu_cov_range
		   as.double(as.vector(simu_t_range)),	# double* pfSimu_t_range
		   as.integer(debug),
		   as.integer(err) );
		   
	return(err);		   
}

bls.plink<-function( phe.file, tped.file, tfam.file, model, bRefit=TRUE,    
	    nMaxIter    = 2000,
	    fBurnInRound= 0.3,
	    fRhoTuning  = 0.095,
	    fQval.add   = 0.01,
	    fQval.dom   = 0.02,
	    debug=F)
{
	r <- .Call("bls_plink", 
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
	if(!is.null(r) && !is.na(r))
	{
		if (!is.null(r$varsel))
			colnames(r$varsel) <- c("grp", "pos", "add.sig", "add.mu", "add.min", "add.max", "dom.sig", "dom.mu", "dom.min", "dom.max", "h2");
		if (!is.null(r$refit))
			colnames(r$refit) <- c("grp", "pos", "add.sig", "add.mu", "add.min", "add.max", "dom.sig", "dom.mu", "dom.min", "dom.max", "h2");
		if (!is.null(r$varsel_cov))
			colnames(r$varsel_cov) <- c("cov.sig", "cov.mu", "cov.min", "cov.max");
		if (!is.null(r$refit_cov))
			colnames(r$refit_cov) <- c("cov.sig", "cov.mu", "cov.min", "cov.max");
		class(r) <- "BLS.ret";
	}
	
	return(r);		   
		   
}

bls.simple<-function(phe.file, snp.file, model, bRefit=T,    
	    nMaxIter = 2000,
	    fBurnInRound = 0.3,
	    fRhoTuning = 0.095,
	    fQval.add   = 0.01,
	    fQval.dom   = 0.02,
	    debug=F)
{
	r <- .Call("bls_simple", 
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
	
	if(!is.null(r) && !is.na(r))
	{
		if (!is.null(r$varsel))
			colnames(r$varsel) <- c("grp", "pos", "add.sig", "add.mu", "add.min", "add.max", "dom.sig", "dom.mu", "dom.min", "dom.max", "h2");
		if (!is.null(r$refit))
			colnames(r$refit) <- c("grp", "pos", "add.sig", "add.mu", "add.min", "add.max", "dom.sig", "dom.mu", "dom.min", "dom.max", "h2");
		if (!is.null(r$varsel_cov))
			colnames(r$varsel_cov) <- c("cov.sig", "cov.mu", "cov.min", "cov.max");
		if (!is.null(r$refit_cov))
			colnames(r$refit_cov) <- c("cov.sig", "cov.mu", "cov.min", "cov.max");
		class(r) <- "BLS.ret";
	}
	
	return(r);		   
}

summary_output<-function(re1, re2)
{
	cat("(1) Covariate Estimate:\n");
		
	for(i in 1:NROW(re1))
		cat(sprintf("%s \t %s \t%.3f\t(%.3f,%.3f)\n", rownames(re1)[i], 
		    ifelse(re1[i,1], "Yes", "---"), re1[i,2], re1[i,3], re1[i,4]));  

	cat("(2) Significant SNPs Estimate:\n");

	cat("    SNP Name\tGrp/Pos\tAdd\tMedian    95%CI  \tDom\tMedian    95%CI  \tH2\n");
	for(i in 1:NROW(re2))
	{
		cat(sprintf("%12s\t%d/%d\t", rownames(re2)[i], re2[i,1], re2[i,2] ));
		if (re2[i,3])
		    cat(sprintf("%s\t%.3f (%.3f,%.3f)", "Yes", re2[i,4], re2[i,5], re2[i,6])) 
		else
		    cat(sprintf("---\t----- (-----, -----)" ));  

		if (re2[i,7])
		    cat(sprintf("\t%s\t%.3f (%.3f,%.3f)", "Yes", re2[i,8], re2[i,9], re2[i,10]))  
		else
		    cat(sprintf("\t---\t----- (-----, -----)" ));  
		
		cat(sprintf("\t%.3f\n", re2[i,11]));  
	}	
}

summary.BLS.ret<-function(r.bls)
{
	if(!is.null(r.bls$refit))
	{
		cat("Refit Result\n");
		summary_output(r.bls$refit_cov, r.bls$refit)
	}
	else if(!is.null(r.bls$varsel))
	{
		cat("Variable Selection Result\n");
		summary_output(r.bls$varsel_cov, r.bls$varsel); 
	}
}

bls.outputpdf<-function( r.bls, fig.file, bRefit=T )
{
	adh2 <- NA;
	if (bRefit)
	{
		if( !is.null(r.bls$refit) )
		{
			adh2 <- r.bls$refit[,c(4,8,11)];
		}
		else
			stop("No refit results\n");		
	}
	else
	{
		if( !is.null(r.bls$varsel) )
		{
			adh2 <- r.bls$varsel[,c(4,8,11)];
		}
		else
			stop("No varible selection results\n");		
	}
	
	pdf(fig.file, width=6, height=6);

	par(mar=c(4.5, 4, 0.5, 2) + 0.1);
	plot.new();
  	par(mfrow=c(3,1));
	
	par(mfg=c(1, 1));
	draw_snplist( adh2[,1], "Estimated additive effect", sigpos_list=NULL);
	par(mfg=c(2, 1));
	draw_snplist( adh2[,2], "Estimated dominant effect", sigpos_list=NULL);
	par(mfg=c(3, 1));
	draw_snplist( adh2[,3], "Heritability", sigpos_list=NULL);
	
	dev.off();
}

#--------------------------------------------------------------
# plot_adh2
# 
# Input:pvs[,1]  snp_name
#       pvs[,2]  chromoseom no
#       pvs[,3]  position
#       pvs[,4]  Additive
#       pvs[,5]  Dominant
#       pvs[,6]  H2
# fig.file: pdf file name
#
# Used by: BLS
#--------------------------------------------------------------
plot.BLS.ret<-function( r.bls, bRefit=T )
{
	adh2 <- NA;
	if (bRefit)
	{
		if( !is.null(r.bls$refit) )
		{
			adh2 <- r.bls$refit[,c(4,8,11)];
		}
		else
			stop("No refit results\n");		
	}
	else
	{
		if( !is.null(r.bls$varsel) )
		{
			adh2 <- r.bls$varsel[,c(4,8,11)];
		}
		else
			stop("No varible selection results\n");		
	}
	
	par(mar=c(4.5, 4, 0.5, 2) + 0.1);
	plot.new();
  	par(mfrow=c(3,1));
	
	par(mfg=c(1, 1));
	draw_snplist( adh2[,1], "Estimated additive effect", sigpos_list=NULL);
	par(mfg=c(2, 1));
	draw_snplist( adh2[,2], "Estimated dominant effect", sigpos_list=NULL);
	par(mfg=c(3, 1));
	draw_snplist( adh2[,3], "Heritability", sigpos_list=NULL);
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