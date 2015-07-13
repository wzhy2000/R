\name{fg.snpscan}
\alias{fg.snpscan}
\title{
SNP scaning
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
fg.snpscan( fg.dat, fgwas.filter = F, options = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fg.dat}{}
  \item{fgwas.filter}{}
  \item{options}{}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (fg.dat, fgwas.filter = F, options = NULL) 
{
    check_FG_DAT(fg.dat)
    options0 <- get_default_options()
    if (!missing(options)) 
        options0[names(options)] <- options
    options <- options0
    r.filter <- list()
    r <- list()
    subset_simple_op <- function(snpmat, sub.idx) {
        return(snpmat[sub.idx, , drop = F])
    }
    subset_plink_op <- function(snpmat, sub.idx) {
        return(snpmat[sub.idx, , drop = F])
    }
    if (fgwas.filter) {
        r.filter <- snpmat_fgwas_filter(fg.dat$phe.mat, fg.dat$snp.mat, 
            "Y", "T", colnames(fg.dat$pheX), options$nParallel.cpu, 
            options$fgwas.cutoff)
        if (r.filter$error) 
            stop(r.filter$err.info)
        r <- mle_parallel(NROW(r.filter$snp.mat), subset_simple_op, 
            r.filter$snp.mat, fg.dat$pheY, fg.dat$pheX, fg.dat$pheT, 
            fg.dat$phe.est, options$nParallel.cpu, options$debug)
    }
    else {
        r <- mle_parallel(NROW(fg.dat$snp.mat), subset_simple_op, 
            fg.dat$snp.mat, fg.dat$pheY, fg.dat$pheX, fg.dat$pheT, 
            fg.dat$phe.est, options$nParallel.cpu, options$debug)
    }
    r.scan <- list(filter = r.filter, snp.range = fg.dat$snp.range, 
        phe.est = fg.dat$phe.est, results = r, top = proc_top_select(r, 
            10), options = options)
    class(r.scan) <- "FGWAS.SCAN"
    return(r.scan)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
