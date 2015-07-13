\name{fgwas.plink}
\alias{fgwas.plink}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
fgwas.plink(file.phe, file.plink.bed, file.plink.bim, file.plink.fam, Y.prefix, Z.prefix, covariate.names, curve = NA, covariance = NA, fgwas.filter = FALSE, options = NULL, force.split = FALSE, plink.command = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file.phe}{}
  \item{file.plink.bed}{}
  \item{file.plink.bim}{}
  \item{file.plink.fam}{}
  \item{Y.prefix}{}
  \item{Z.prefix}{}
  \item{covariate.names}{}
  \item{curve}{}
  \item{covariance}{}
  \item{fgwas.filter}{}
  \item{options}{}
  \item{force.split}{}
  \item{plink.command}{}
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
function (file.phe, file.plink.bed, file.plink.bim, file.plink.fam, 
    Y.prefix, Z.prefix, covariate.names, curve = NA, covariance = NA, 
    fgwas.filter = FALSE, options = NULL, force.split = FALSE, 
    plink.command = NULL) 
{
    cat("[ fGWAS PLINK ] Procedure.\n")
    cat("Checking the parameters ......\n")
    if (missing(file.phe) || missing(file.plink.bed) || missing(file.plink.bim) || 
        missing(file.plink.fam) || missing(Y.prefix) || missing(Z.prefix) || 
        missing(covariate.names)) 
        stop("! file.phe, file.plink.bed, file.plink.bim, file.plink.fam, Y.prefix, Z.prefix and covariate.names must be assigned with the valid values.")
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
    if (!(is.logical(force.split) && length(force.split) == 1)) 
        stop("! The parameter of force.split should be a logical value(TRUE or FALSE).")
    cat("* Phenotypic Data File = ", file.phe, "\n")
    cat("* PLINK BED File = ", file.plink.bed, "\n")
    cat("* PLINK BIM File = ", file.plink.bim, "\n")
    cat("* PLINK FAM File = ", file.plink.fam, "\n")
    cat("* PLINK Command = ", plink.command, "\n")
    cat("* Force Split by PLINK Command = ", force.split, "\n")
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
    r.plink <- list()
    r.plink$options <- options
    r.plink$params <- list(file.phe = file.phe, file.plink.bed = file.plink.bed, 
        file.plink.bim = file.plink.bim, file.plink.fam = file.plink.fam, 
        Y.prefix = Y.prefix, Z.prefix = Z.prefix, covariate.names = covariate.names, 
        fgwas.filter = fgwas.filter)
    r.phe <- read_simple_phenotype(file.phe, Y.prefix, Z.prefix, 
        covariate.names)
    if (r.phe$error) 
        stop(r.phe$err.info)
    r.plink$phe.mat <- r.phe$phe.mat
    r.est <- fg.estimate(r.plink$phe.mat, Y.prefix, Z.prefix, 
        covariate.names, curve, covariance)
    if (r.est$error) 
        stop(r.est$err.info)
    else {
        if (is.na(curve)) 
            curve <- r.est$curve
        if (is.na(covariance)) 
            covariance <- r.est$covariance
        r.plink$curve <- curve
        r.plink$covariance <- covariance
        r.plink$est.values <- r.est$est.values
    }
    pd <- list()
    if (force.split || !try_load_plink(file.plink.bed, file.plink.bim, 
        file.plink.fam)) {
        r.filter <- plink_fgwas_bigdata(file.plink.bed, file.plink.bim, 
            file.plink.fam, file.phe, plink.command, Y.prefix, 
            Z.prefix, covariate.names, options$nParallel.cpu, 
            options$fgwas.cutoff, "GLS")
        if (r.filter$error) 
            stop(r.filter$err.info)
        fgwas.filter <- TRUE
        r.plink$filter <- r.filter$r.fgwas
        pd <- r.filter$snp.mat
    }
    else {
        pd <- load_plink_binary(file.plink.bed, file.plink.bim, 
            file.plink.fam, file.phe)
        if (is.null(pd)) 
            stop("Failed to load PLINK dataset!")
        if (fgwas.filter) {
            r.filter <- plink_fgwas_filter(pd, Y.prefix, Z.prefix, 
                covariate.names, options$nParallel.cpu, options$fgwas.cutoff, 
                "GLS")
            if (r.filter$error) 
                stop(r.filter$err.info)
            pd <- r.filter$snp.mat
            r.plink$filter <- r.filter$r.fgwas
        }
    }
    if (fgwas.filter) {
        subset_op <- function(snpmat, sub.idx) {
            return(snpmat[sub.idx, , drop = F])
        }
        r.fgwas <- fgscan_parallel(NROW(pd), subset_op, pd, r.plink$phe.mat, 
            Y.prefix, Z.prefix, covariate.names, options$debug, 
            options$nParallel.cpu)
    }
    else {
        subset_op <- function(snpmat, sub.idx) {
            snp.sub <- get_plink_subsnp(snpmat, sub.idx)
            snp.mat <- cbind(snp.sub$info[, c(2, 3)], snp.sub$snp)
            return(snp.mat)
        }
        r.fgwas <- fgscan_parallel(NCOL(pd$snp.mat$genotypes), 
            subset_op, pd$snp.mat, pd$phe.mat, Y.prefix, Z.prefix, 
            covariate.names, options$debug, options$nParallel.cpu)
    }
    if (!is.null(r.fgwas) && !is.na(r.fgwas)) {
        r.plink$fgwas <- r.fgwas
        return(r.plink)
    }
    else {
        cat("! No results\n")
        r.plink$fgwas <- "try-error"
        return(r.plink)
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
