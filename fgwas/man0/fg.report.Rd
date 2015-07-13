\name{fg.report}
\alias{fg.report}
\title{
Report
}
\description{
}
\usage{
fg.report( fg.dat, fg.scan = NULL, fg.perm = NULL, file.pdf = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fg.dat}{}
  \item{fg.scan}{}
  \item{fg.perm}{}
  \item{file.pdf}{}
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
function (fg.dat, fg.scan = NULL, fg.perm = NULL, file.pdf = NULL) 
{
    check_FG_DAT(fg.dat)
    if (!missing(fg.scan)) 
        check_FG_SCAN(fg.scan)
    if (!missing(fg.perm)) 
        check_FG_PERM(fg.perm)
    if (!missing(file.pdf)) 
        check_file_pdf(file.pdf)
    if (is.null(file.pdf)) 
        file.pdf <- paste(fg.dat$file.pheY.csv, ".PDF", sep = "")
    Report.new(file.pdf, options)
    Report.title("Functional GWAS Report", "fGWAS", "http://ccb.bjfu.edu.cn/")
    Report.par("dat.file", fg.dat$file.pheY.csv)
    Report.AddHeadline("Data Summary", level = 1)
    proc_report_dat(fg.dat)
    if (!is.null(fg.scan)) {
        Report.AddHeadline("SNP Scan", level = 1)
        proc_report_snpscan(fg.dat, fg.scan, fg.perm)
    }
    if (!is.null(fg.perm)) {
        Report.AddHeadline("SNP Permutation", level = 1)
        proc_report_perm(fg.dat, fg.perm)
    }
    if (!is.null(fg.perm) && !is.null(fg.scan)) {
        r.sig <- fg_detect_sig(fg.dat, fg.scan, fg.perm)
        Report.AddHeadline("Significant SNP", level = 1)
        proc_report_sig(fg.dat, r.sig$sig.05, r.sig$sig.01)
    }
    if (is.null(fg.perm) && !is.null(fg.scan)) {
        Report.AddHeadline("Top SNP", level = 1)
        proc_report_topn(fg.dat, fg.scan)
    }
    Report.Output(file.pdf)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
