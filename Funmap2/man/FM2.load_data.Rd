\name{FM2.load_data}
\alias{FM2.load_data}

\title{Data loading from real files}

\description{ 
  Load the experimental data from real files.
} 

\usage{ FM2.load_data(pheno_file, geno_file, marker_file, model_type, cross_type, covar_type=COVAR_AR1 ) }

\arguments{
\item{pheno_file}{The phenotype file.The format is described in the details section.}
\item{geno_file}{The genotype file.The format is described in the details section.}
\item{marker_file}{The marker file.The format is described in the details section.}
\item{model_type}{5 models(curves) can be selected. CURVE_LC, CURVE_BI, CURVE_PC, CURVE_NP or CURVE_EXP. see \code{\link{FM2.param}}.}
\item{cross_type}{CROSS_BC, CROSS_F2 or CROSS_RIL. see \code{\link{FM2.param}}.}
\item{cross_type}{COVAR_AR1, COVAR_SAD2. see \code{\link{FM2.param}}.}
}

\value{
A data object can be returned. Its structure is same as the simulation data, 
see \code{\link{FM2.simulate}} for details. You can use summary(\code{\link{summary.FM2.dat}}) or \code{str} command to check it.

}

\details{

The phenotype file, genotype file and marker file must be a CSV file. 

1) \bold{The phenotype file.}
The first column is individual ID and the rest columns are sample data 
for every measurement. It looks like the following file.

\preformatted{
ID.	1st	2nd	3rd	4th	5th	6th	7th
1	2.9033	4.118	6.1495	7.8161	9.8379	12.963	14.918
2	4.3306	5.3783	7.0647	9.3624	11.439	13.379	15.701
3	2.3997	4.052	5.5431	7.6933	9.8471	11.249	12.849
4	3.3044	4.154	5.8924	7.7133	9.2144	10.945	12.403
...
}

2) \bold{The genotype file.}
The first column is individual ID and the rest columns are 
genotype data for every marker. An example is shown in the following table.
Three genotypes (aa=0, Aa=1, AA=2) and missing data(coded as -1)
are valid maker values.

\preformatted{
ID.	marker1	marker2	marker3	marker4	marker5	marker6
1	1	1	0	1	-1	0
2	1	1	1	1	0	0
3	1	1	1	0	1	1
...
}

3) \bold{The marker list file.}
The first column is marker's ID, the rest columns are the marker's name,
distance, group index, and group name for every marker.In the marker file,
the distance field is a distance (in cM) in one chromosome or linkage group. 
The header row should be included.

\preformatted{
id	Marker	Dist.	Grp_idx	Group
1	marker1	0	1	G1
2	marker2	20	1	G1
3	marker3	40	1	G1
...
}

}

\author{Zhong Wang, \email{zzw2@psu.edu} }

\examples{
if(0)
{
	dat <- FM2.load_data("pheno.csv","geno.csv","marker.csv", CURVE_LC, CROSS_BC, COVAR_AR1);
	summary( dat );
	plot(dat);
}
}

\keyword{functions}
