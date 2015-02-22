###############################################################
# 
# System management utility
#
# History:
# 03/17/2010 Version 0.1
#
###############################################################


#--------------------------------------------------------------
# FM2.sys
#
#--------------------------------------------------------------
FM2.sys<-function( module )
{
	i_module    <- module
	i_seq_id    <- 0
	i_start_timer<-proc.time()
	i_selpeaks	<- NULL
	i_hash      <- NULL
	
	get_seq_id<-function()
	{
		i_seq_id <<- i_seq_id+1;
		return (i_seq_id);
	}

	reset<-function()
	{
		i_seq_id <<- 0;
		i_start_timer<<-0;

		i_hash <<- new.env(hash=TRUE, parent=emptyenv(), size=100L);

		assign("plot_doctype",  "pdf",   i_hash);
		assign("scan_step",     "2",     i_hash);
		assign("peak_count",    "5",     i_hash);
		assign("permu_loop",    "1000",  i_hash);
		assign("cluster_count", "1",     i_hash);
		assign("np.order",      "6",     i_hash);
	}

	task_start<-function(...)
	{
		i_start_timer<<-proc.time();
		
		msgs <- list(...);
		if (length(msgs)==0)
			cat("The task is started...\r\n")
		else	
		{
			cat(..., sep="");
		}
		
		flush.console();
	}

	task_elapsed<-function(finished = NA, ...)
	{
		nt <- proc.time() - i_start_timer;
		
		if(is.na(finished))
		{
		}
		else
			sys_t <-  round(c(nt[3], nt[3]*( 1 - finished)/finished ));
		
		sTime1<-"";
		sTime2<-"";
		if (sys_t[1]>60)
			sTime1 <- sprintf( "%02g:%02g:%02g", sys_t[1]%/%3600, (sys_t[1]- sys_t[1]%/%3600*3600)%/%60, sys_t[1]%%60 )
		else
			sTime1 <- sprintf( "%02g seconds", sys_t[1] );

		if (sys_t[2]>60)
			sTime2 <- sprintf( "%02g:%02g:%02g", sys_t[2]%/%3600, (sys_t[2]- sys_t[2]%/%3600*3600)%/%60, sys_t[2]%%60 )
		else
			sTime2 <- sprintf( "%02g seconds", sys_t[2] );
		sTime<-paste( sTime1, " has elapsed, left time: ", sTime2, ". ");
		
		msgs <- list(...);
		if ( length(msgs) != 0 )
		{
			for (i in 1:length(msgs))
				if (msgs[[i]]=="$SYS_PROMPT$")
					msgs[[i]] <- sTime;
		}
		else
			msgs[[1]]<-paste( sTime, "\r\n");

		cat( unlist(msgs), sep="" );
		flush.console();
	}

	task_stop<-function(...)
	{
		msgs <- list(...);
		if (length(msgs)==0)
			cat("The task is stopped...\r\n")
		else	
			cat( ... , sep="" );

		flush.console();
	}

	is_LR_peak<-function( grp_idx, qtl ) 
	{
		if (!is.vector(i_selpeaks))
		{
			sel <- which(i_selpeaks[,1]==grp_idx & i_selpeaks[,2]==qtl);
			return(length(sel)>0);
		}
		else
		{
			return(i_selpeaks[1]==grp_idx && i_selpeaks[2]==qtl)
		}
	}
	
	set_LR_peaks<-function( peaks )
	{
		i_selpeaks <<- peaks;
	}
	
	get_LR_peaks<-function()
	{
		return (i_selpeaks);
	}

	set_value<-function( key, value)
	{
		old_value <- NULL;
		try( old_value <-  get( key, i_hash) , silent = TRUE);

		assign( key, value, i_hash);
		return (old_value);
	}

	get_value<-function( key )
	{
		return ( get( key, i_hash) );
	}


	return (
		list(reset         = reset, 
			seq_id        = get_seq_id, 
			task_start    = task_start, 
			task_stop     = task_stop, 
			task_elapsed  = task_elapsed,
			get_LR_peaks  = get_LR_peaks,
			set_LR_peaks  = set_LR_peaks,
			is_LR_peak    = is_LR_peak,
			set_value     = set_value,
			get_value     = get_value ))
}
 
#sys<-FM.sys("abc");
#> sys$seq_id();
#[1] 1
#> sys$seq_id();
#[1] 2
#> sys$seq_id();
#[1] 3
#> sys$reset();
#> sys$seq_id();
#[1] 0
