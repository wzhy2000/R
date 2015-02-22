
#--------------------------------------------------------------
# eval.execute
#
# Evaluate a MOEL by the specific parameters and loops .
#--------------------------------------------------------------
eval.execute<-function(par_obj, nLoop)
{
	eval <- list(par = par_obj);
	eval$raw_dat <- list();
	eval$qtl_par <- list();
	
	for( i in 1:nLoop)
	{
		dat <- FM2.simulate(par_obj);
		eval$raw_dat[[i]]<- dat$phenos_table;
		
		ret <- FM2.qtlmodel(dat);
		eval$qtl_par[[i]]<- c( ret$qtl_LR, ret$qtl_pos, c( unlist(ret$curve_par ) ) );
	}
	
	me <- array(0, dim=c( nLoop, length(eval$qtl_par[[1]])));
	for (i in 1:nLoop)
	{
		me[i,]<- eval$qtl_par[[i]];
	}

	eval$qtl_par <- me;

	n<-length(eval$qtl_par[1,]);
	m<-length(eval$qtl_par[,1]);
	sums<-array(0, dim=c(n,4));
	
	par_v<-c( 0, eval$par$simu_qtlpos, eval$par$simu_rho, eval$par$simu_s2);
	if (!is.null(eval$par$QQ2)) par_v<-c(par_v, unlist(eval$par$QQ2));
	if (!is.null(eval$par$Qq1)) par_v<-c(par_v, unlist(eval$par$Qq1));
	if (!is.null(eval$par$qq0)) par_v<-c(par_v, unlist(eval$par$qq0));

	for (i in 1:n)
	{
		sums[i,1]<-par_v[i];
		sums[i,2]<-mean(eval$qtl_par[,i]);
		sums[i,3]<-sqrt(var( eval$qtl_par[,i] ) );
		sums[i,4]<-sqrt(sum( (eval$qtl_par[,i] - par_v[i])^2 )/m)/sqrt(m)
	}

	sums[1,4]<-NA;
	sums[1,1]<-NA;
	rownames0<-c("MLE_LR2", "qtl_pos", "rho", "s2");
	if (!is.null(eval$par$QQ2)) rownames0<-c(rownames0, unlist(names(eval$par$QQ2)));
	if (!is.null(eval$par$Qq1)) rownames0<-c(rownames0, unlist(names(eval$par$Qq1)));
	if (!is.null(eval$par$qq0)) rownames0<-c(rownames0, unlist(names(eval$par$qq0)));
	
	rownames(sums)<-rownames0;
	colnames(sums)<-c("Simu.", "mean(MLE)", "SD", "SE" );
	
	eval$sum <- sums;
	
	class(eval) <- "FM2.eval.ret";
	
	return(eval);
}

#--------------------------------------------------------------
# eval.summary
#
# 
#--------------------------------------------------------------
eval.summary<-function( eval_obj, f_curve_mu )
{
	sums <- eval_obj$sum;
	
	m<-length(eval_obj$qtl_par[,1]);

	cat("Results evaluated by ", m , " times test.\n");
	show(sums);
	
	return ();
}

#--------------------------------------------------------------
# eval.plot
#
# 
#--------------------------------------------------------------
eval.plot<-function( eval_obj )
{
	fpt.plot_eval_ret(eval_obj);

}


