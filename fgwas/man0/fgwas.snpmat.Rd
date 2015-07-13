\name{fgwas.snpmat}
\alias{fgwas.snpmat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
fgwas.snpmat(phe.mat, snp.mat, Y.prefix, Z.prefix, covariate.names, curve = NA, covariance = NA, fgwas.filter = FALSE, options = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{phe.mat}{}
  \item{snp.mat}{}
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
function (phe.mat, snp.mat, Y.prefix, Z.prefix, covariate.names, 
    curve = NA, covariance = NA, fgwas.filter = FALSE, options = NULL) 
{
    cat("[ fGWAS SNPMAT ] Procedure.\n")
    cat("Checking the parameters ......\n")
    if (missing(phe.mat) || missing(snp.mat) || missing(Y.prefix) || 
        missing(Z.prefix) || missing(covariate.names)) 
        stop("! phe.mat, snp.mat, file.plink.tfam, Y.prefix, Z.prefix and covariate.names must be assigned with the valid values.")
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
    cat("* Phenotypic Matrix = ", dim(phe.mat), "\n")
    cat("* SNP Matrix = ", dim(snp.mat), "\n")
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
    if (class(phe.mat) == "data.frame") {
        cat("Phenotypic data frame is converted to the matrix class.\n")
        phe.colnames <- colnames(phe.mat)
        phe.rownames <- rownames(phe.mat)
        phe.mat <- matrix(as.numeric(as.matrix(phe.mat, rownames.force = NA)), 
            ncol = NCOL(phe.mat))
        colnames(phe.mat) <- phe.colnames
        rownames(phe.mat) <- phe.rownames
    }
    r.snpmat <- list()
    r.snpmat$options <- options
    r.snpmat$params <- list(file.phe = file.phe, file.snp = file.snp, 
        Y.prefix = Y.prefix, Z.prefix = Z.prefix, covariate.names = covariate.names, 
        fgwas.filter = fgwas.filter)
    r.snpmat$phe.mat <- phe.mat
    r.est <- fg.estimate(r.snpmat$phe.mat, curve, covariance)
    if (r.est$error) 
        stop(r.est$err.info)
    else {
        if (is.na(curve)) 
            curve <- r.est$curve
        if (is.na(covariance)) 
            covariance <- r.est$covariance
        r.snpmat$curve <- curve
        r.snpmat$covariance <- covariance
        r.snpmat$est.values <- r.est$est.values
    }
    if (fgwas.filter) {
        r.filter <- snpmat_fgwas_filter(phe.mat, snp.mat, Y.prefix, 
            Z.prefix, covariate.names, options$nParallel.cpu, 
            options$fgwas.cutoff, "GLS")
        if (r.filter$error) 
            stop(r.filter$err.info)
        r.snpmat$filter <- r.filter$r.fgwas
        snp.mat <- r.filter$snp.mat
    }
    subset_op <- function(snpmat, sub.idx) {
        return(snpmat[sub.idx, , drop = F])
    }
    r.fgwas <- snpmat_parallel(NROW(r.snpmat$snp.mat), subset_op, 
        r.snpmat$snp.mat, phe.mat, Y.prefix, Z.prefix, covariate.names, 
        options$debug, options$nParallel.cpu)
    if (!is.null(r.fgwas) && !is.na(r.fgwas)) {
        r.snpmat$fgwas <- r.fgwas
        return(r.snpmat)
    }
    else {
        cat("! No results\n")
        r.snpmat$fgwas <- "try-error"
        return(r.snpmat)
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
