\name{fg.permutation}
\alias{fg.permutation}
\title{
Permutation
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
fg.permutation(fg.dat, n.perm, file.rdata.perm = NULL, options = list( n.loop = 5, debug = F))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fg.dat}{}
  \item{n.perm}{}
  \item{file.rdata.perm}{}
  \item{method}{}
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
function (fg.dat, n.perm, file.rdata.perm = NULL, options = list(cache.reset = T, 
    n.loop = 5, debug = F)) 
{
    check_FG_DAT(fg.dat)
    check_n_perm(n.perm)
    if (!missing(file.rdata.permp)) 
        check_rdata_perm(file.rdata.perm)
    if (missing(options)) 
        options <- get_default_options()
    else {
        options0 <- get_default_options()
        options0[names(options)] <- options
        options <- options0
    }
    cat(" PERM START n.perm=", n.perm, "\n")
    r.perm <- c()
    for (i in 1:n.perm) {
        p0 <- fg_permu_core(fg.dat, n.loop)
        r.perm <- rbind(r.perm, c(p0$maxLR2, p0$perm.sample))
        cat(" PERM [", i, "] max.LR2=", p0$maxLR2, "\n")
        if (!is.null(file.rdata.perm)) 
            save(r.perm, file = file.rdata.perm)
    }
    p_cut <- fg_permu_cutoff(r.perm)
    if (!is.null(p_cut$pcut.05)) 
        cat(" PERM cutoff(pv-0.05) = ", p_cut$pcut.05, "\n")
    if (!is.null(p_cut$pcut.01)) 
        cat(" PERM cutoff(pv-0.01) = ", p_cut$pcut.01, "\n")
    fg.perm <- list(error = F, pcut.05 = p_cut$pcut.05, pcut.01 = p_cut$pcut.01, 
        pcut.tab = p_cut$pcut.tab, max.lr2 = r.perm[, 1], full = r.perm)
    class(fg.perm) = "FGWAS.PERM"
    return(fg.perm)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
