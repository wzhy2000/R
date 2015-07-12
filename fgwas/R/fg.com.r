
get_regvalue<-function(var.name, def=0)
{
	return(def);
}

reset_seed<-function()
{	
	fg.sys$n.seed <<- fg.sys$n.seed + runif(1, 1, 5);
	set.seed( fg.sys$n.seed );
	if(fg.sys$log>=LOG_DEBUG) cat(".");
}

get_con_param<-function(parm.id)
{
	for (e in commandArgs()) 
	{
		ta = strsplit(e,"=", fixed=TRUE);
		if(! is.na( ta[[1]][2])) 
		{
			temp = ta[[1]][2];
			if( ta[[1]][1] == parm.id) {
				temp = as.integer(temp);
				return (temp);
			}
		}
	}

	return(NA);
}

remove_extname<-function(filename)
{
	rs <- unlist(strsplit(filename, "\\."))
	if (length(rs)==1)
		return(filename);

	rs.last <- rs[length(rs)];
	
	if (grepl(pattern = "[\\/]", x = rs.last))
		return(filename);
		
	rs <- rs[-length(rs)];
	file.noext <- paste(rs, collapse=".");
	
	return(file.noext);
}

get_default_options<-function()
{
	options=list(nParallel.cpu = 0, cache.reset=T, n.loop=5, fgwas.cutoff=0.05, debug=F)
	return(options);	
}

show_options<-function(options)
{
	cat( "* Parallel Computing: ", ifelse( options$nParallel.cpu>1, "Yes,", "No,"), options$nParallel.cpu,  "CPU(s)\n");
	cat( "* Cache Reset: ",  ifelse( options$cache.reset, "Yes", "No"),  "\n");
	cat( "* Optim iteration.: ",  options$n.loop, "\n");
	cat( "* Debug Output: ", ifelse( options$debug, "Yes", "No"),"\n");

	if(options$debug) fg.sys$log <<- LOG_DEBUG else fg.sys$log <<- LOG_INFO;
}
