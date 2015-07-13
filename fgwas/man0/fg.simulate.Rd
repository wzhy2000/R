\name{fg.simulate}
\alias{fg.simulate}
\title{
Data simulation
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
fg.simulate(file.prefix, curve, covariance, n.obs, n.snp, time.points, par0 = NA, par1 = NA, par2 = NA, par.covar = NA, par.X = NA, phe.missing = 0.03, snp.missing = 0.03, sig.pos = NA, plink.format = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file.prefix}{}
  \item{curve}{String indicating curve name.}
  \item{covariance}{String indicating covariance name.}
  \item{n.obs}{}
  \item{n.snp}{}
  \item{time.points}{}
  \item{par0}{}
  \item{par1}{}
  \item{par2}{}
  \item{par.covar}{}
  \item{par.X}{}
  \item{phe.missing}{}
  \item{snp.missing}{}
  \item{sig.pos}{}
  \item{plink.format}{}
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
function (file.prefix, curve, covariance, n.obs, n.snp, time.points, 
    par0 = NA, par1 = NA, par2 = NA, par.covar = NA, par.X = NA, 
    phe.missing = 0.03, snp.missing = 0.03, sig.pos = NA, plink.format = FALSE) 
{
    library(mvtnorm)
    check_n_obs(n.obs)
    check_n_snp(n.snp)
    check_no_curve(no.curve)
    check_no_covar(no.covar)
    check_time_points(time.points)
    if (!missing(par.qq)) 
        check_par_qq(par.qq)
    if (!missing(par.Qq)) 
        check_par_Qq(par.Qq)
    if (!missing(par.QQ)) 
        check_par_QQ(par.QQ)
    if (!missing(par.covar)) 
        check_par_covar(par.covar)
    if (!missing(par.X)) 
        check_par_X(par.X)
    if (!missing(file.prefix)) 
        check_file_data(file.prefix)
    if (!missing(sig.pos)) 
        check_sig_pos(sig.pos)
    if (is.null(sig.pos)) {
        sig.pos <- round(runif(1, n.snp * 0.25, n.snp * 0.75))
        cat(" * A significant SNP is randomly specified to location(", 
            sig.pos, ")\n")
    }
    if (class(curve) == "fgwas.curve") 
        fg_curve <- curve
    else fg_curve <- fg.getCovar(curve)
    if (missing(par0) && missing(par1) && missing(par2)) {
        par0 <- fg_curve@par_simu[1, ]
        par1 <- fg_curve@par_simu[2, ]
        par2 <- fg_curve@par_simu[3, ]
    }
    if (class(covariance) == "fgwas.covar") 
        fg_covar <- covariance
    else fg_covar <- fg.getCovar(covariance)
    if (missing(par.covar)) 
        par.covar <- fg_covar@par_simu
    fg.dat <- proc_dat_simu(n.obs, n.snp, par.X, fg_curve, fg_covar, 
        time.points, sig.pos, snp.missing, phe.missing)
    fg.dat$file.phe.out <- paste(file.prefix, ".phe.csv", sep = "")
    write.csv(data.frame(ID = fg.dat$ids, fg.dat$pheX, fg.dat$pheY), 
        file = fg.dat$file.phe.out, quote = F, row.names = F)
    if (!plink.format) {
        fg.dat$file.geno.dat <- paste(file.prefix, ".geno.dat", 
            sep = "")
        write.table(data.frame(fg.dat$snp.info, fg.dat$gen), 
            file = fg.dat$file.geno.dat, quote = F, row.names = F)
        return(list(err = err, file.simple.phe = fg.dat$file.phe.out, 
            file.simple.snp = fg.dat$file.snp.out))
    }
    else {
        r <- convert_simpe_to_plink(data.frame(fg.dat$snp.info, 
            fg.dat$gen), file.snp.out)
        return(list(err = err, file.simple.phe = fg.dat$file.phe.out, 
            file.plink.bed = r$file.plink.bed, file.plink.bim = r$file.plink.bim, 
            file.plink.fam = r$file.plink.fam))
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
