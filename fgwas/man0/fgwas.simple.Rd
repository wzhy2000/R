\name{fgwas.simple}
\alias{fgwas.simple}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
fgwas.simple(file.phe, file.snp, Y.prefix, Z.prefix, covariate.names, curve = NA, covariance = NA, fgwas.filter = FALSE, options = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file.phe}{}
  \item{file.snp}{}
  \item{Y.prefix}{}
  \item{Z.prefix}{}
  \item{covariate.names}{}
  \item{curve}{}
  \item{covariance}{}
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
function (file.phe, file.snp, Y.prefix, Z.prefix, covariate.names, 
    curve = NA, covariance = NA, fgwas.filter = FALSE, options = NULL) 
{
    cat("[ fGWAS SIMPLE ] Procedure.\n")
    cat("Checking the parameters ......\n")
    if (missing(file.phe) || missing(file.snp) || missing(Y.prefix) || 
        missing(Z.prefix) || missing(covariate.names)) 
        stop("! file.phe, file.snp, Y.prefix, Z.prefix and covariate.names must be assigned with the valid values.")
    if (!(is.character(Y.prefix) && length(Y.prefix) == 1)) 
        stop("! The parameter of Y.prefix should be assigned with a prefix of outcome column in the phenotypic data.")
    if (!(is.character(Z.prefix) && length(Z.prefix) == 1)) 
        stop("! The parameter of Z.prefix should be assigned with a prefix of time column in the phenotypic data.")
    if (!missing("covariate.names") && length(covariate.names) > 
        0 && !is.character(covariate.names)) 
        stop("! The parameter of covariate.names should be assigned with covariate names in the phenotypic data.")
    if (!(is.logical(fgwas.filter) && length(fgwas.filter) == 
        1)) 
        stop("! The parameter of fgwas.filter should be a logical value(TRUE or FALSE).")
    cat("* Phenotypic Data File = ", file.phe, "\n")
    cat("* Simpe SNP File = ", file.snp, "\n")
    show_fgwas_parameters(curve, covariance, Y.prefix, Z.prefix, 
        covariate.names, fgwas.filter)
    if (missing(options)) 
        options <- get_default_options()
    else {
        options0 <- get_default_options()
        options0[names(options)] <- options
        options <- options0
    }
    cat("Checking the optional items......\n")
    show_options(options)
    r.simple <- list()
    r.simple$options <- options
    r.simple$params <- list(file.phe = file.phe, file.snp = file.snp, 
        Y.prefix = Y.prefix, Z.prefix = Z.prefix, covariate.names = covariate.names, 
        fgwas.filter = fgwas.filter)
    r.phe <- read_simple_phenotype(file.phe, Y.prefix, Z.prefix, 
        covariate.names)
    if (r.phe$error) 
        stop(r.phe$err.info)
    r.simpe$phe.mat <- r.phe$phe.mat
    r.est <- fg.estimate(r.simple$phe.mat, curve, covariance)
    if (r.est$error) 
        stop(r.est$err.info)
    else {
        if (is.na(curve)) 
            curve <- r.est$curve
        if (is.na(covariance)) 
            covariance <- r.est$covariance
        r.simple$curve <- curve
        r.simple$covariance <- covariance
        r.simple$est.values <- r.est$est.values
    }
    r.snp <- read_simple_genotype(file.snp, r.simpe$phe.mat)
    if (r.snp$error) 
        stop(r.snp$err.info)
    r.simpe$snp.mat <- r.snp$snp.mat
    if (fgwas.filter) {
        r.filter <- snpmat_lm_filter(r.simple$phe.mat, r.simple$snp.mat, 
            Y.prefix, Z.prefix, covariate.names, options$nParallel.cpu, 
            options$fgwas.cutoff, "LONG")
        if (r.filter$error) 
            stop(r.filter$err.info)
        r.simple$filter <- r.filter$r.fgwas
        r.simple$org.snp.mat <- r.simple$snp.mat
        r.simple$snp.mat <- r.filter$snp.mat
    }
    subset_op <- function(snpmat, sub.idx) {
        return(snpmat[sub.idx, , drop = F])
    }
    r.fgwas <- fgscan_parallel(NROW(r.simple$snp.mat), subset_op, 
        r.simple$phe.mat, Y.prefix, Z.prefix, covariate.names, 
        options$debug, options$nParallel.cpu)
    r.simple$org.snp.mat <- NULL
    r.simple$snp.mat <- NULL
    if (!is.null(r.fgwas) && !is.na(r.fgwas)) {
        r.simple$fgwas <- r.fgwas
        return(r.simple)
    }
    else {
        cat("! No results\n")
        r.simple$fgwas <- "try-error"
        return(r.simple)
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
