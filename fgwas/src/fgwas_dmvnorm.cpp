#include <R.h>
#include <Rinternals.h>
#include <Rembedded.h>
#include <Rdefines.h>
#include <Rmath.h>
#include <R_ext/Applic.h>
#include <R_ext/Lapack.h>

#include "fm_rlogger.h"
#include "fm_matrix.h"
#include "fm_vector.h"
#include "fm_err.h"

#include "fgwas_dmvnorm.h"

/*
dmvnorm <- function (x, mean, sigma, log = FALSE)
{
    if (is.vector(x)) {
        x <- matrix(x, ncol = length(x))
    }
    if (missing(mean)) {
        mean <- rep(0, length = ncol(x))
    }
    if (missing(sigma)) {
        sigma <- diag(ncol(x))
    }
    if (NCOL(x) != NCOL(sigma)) {
        stop("x and sigma have non-conforming size")
    }
    if (!isSymmetric(sigma, tol = sqrt(.Machine$double.eps))) {
        stop("sigma must be a symmetric matrix")
    }
    if (length(mean) != NROW(sigma)) {
        stop("mean and sigma have non-conforming size")
    }
    distval <- mahalanobis(x, center = mean, cov = sigma)
    logdet <- sum(log(eigen(sigma, symmetric = TRUE, only.values = TRUE)$values))
    logretval <- -(ncol(x) * log(2 * pi) + logdet + distval)/2
    if (log)
        return(logretval)
    exp(logretval)
}

mahalanobis<-function (x, center, cov, inverted = FALSE, ...)
{
    x <- if (is.vector(x))
        matrix(x, ncol = length(x))
    else as.matrix(x)
    x <- sweep(x, 2, center)
    if (!inverted)
        cov <- solve(cov, ...)
    retval <- rowSums((x %*% cov) * x)
    names(retval) <- rownames(x)
    retval
}

eigen <-function (x, symmetric, only.values = FALSE, EISPACK = FALSE)
{
    x <- as.matrix(x)
    if (!is.null(dimnames(x)))
        dimnames(x) <- list(NULL, NULL)
    n <- nrow(x)
    if (!n)
        stop("0 x 0 matrix")
    if (n != ncol(x))
        stop("non-square matrix in 'eigen'")
    complex.x <- is.complex(x)
    if (!complex.x && !is.double(x))
        storage.mode(x) <- "double"
    if (!all(is.finite(x)))
        stop("infinite or missing values in 'x'")
    if (missing(symmetric))
        symmetric <- isSymmetric.matrix(x)
    if (!EISPACK) {
        if (symmetric) {
            z <- if (!complex.x)
                .Call("La_rs", x, only.values, PACKAGE = "base")
            else .Call("La_rs_cmplx", x, only.values, PACKAGE = "base")
            ord <- rev(seq_along(z$values))
        }
        else {
            z <- if (!complex.x)
                .Call("La_rg", x, only.values, PACKAGE = "base")
            else .Call("La_rg_cmplx", x, only.values, PACKAGE = "base")
            ord <- sort.list(Mod(z$values), decreasing = TRUE)
        }
        return(list(values = z$values[ord], vectors = if (!only.values) z$vectors[,
            ord, drop = FALSE]))
    }
}

*/


/* Real, symmetric case of eigen , La_rs */
int get_eigen( CFmMatrix* pFmX, bool only_values, CFmVector* pEigenValue, CFmMatrix* pFmZ)
{
    int n, lwork, info = 0, ov;
    char jobv[1], uplo[1], range[1];
    double *work, *rx, *rvalues, tmp, *rz = NULL;
    int liwork, *iwork, itmp, m;
    double vl = 0.0, vu = 0.0, abstol = 0.0;

    /* valgrind seems to think vu should be set, but it is documented not to be used if range='a' */
    int il, iu, *isuppz;

    uplo[0] = 'L';

    n = pFmX->GetNumRows();
    if ( n != pFmX->GetNumCols() )
    	return( ERR_LAPACK ); // # "'x' must be a square numeric matrix"));

    if ( only_values )
    	jobv[0] = 'N';
    else
    	jobv[0] = 'V';

    /* work on a copy of x, since LAPACK trashes it */
	rx = (double *) R_alloc(n * (size_t) n, sizeof(double));
	Memcpy(rx, pFmX->GetData(), (size_t) n * n);

    pEigenValue->Resize(n);
    rvalues = pEigenValue->GetData();

    range[0] = 'A';
    if (!ov && pFmZ)
    {
		pFmZ->Resize(n, n);
		rz = pFmZ->GetData();
    }

    isuppz = (int *) R_alloc(2*(size_t)n, sizeof(int));

    /* ask for optimal size of work arrays */
    lwork = -1; liwork = -1;
    F77_CALL(dsyevr)(jobv, range, uplo, &n, rx, &n,
		     &vl, &vu, &il, &iu, &abstol, &m, rvalues,
		     rz, &n, isuppz,
		     &tmp, &lwork, &itmp, &liwork, &info);

    if (info != 0)
    	return( ERR_LAPACK ); // #  "error code %d from Lapack routine dsyevr";

    lwork = (int) tmp;
    liwork = itmp;
    work = (double *) R_alloc(lwork, sizeof(double));
    iwork = (int *) R_alloc(liwork, sizeof(int));
    F77_CALL(dsyevr)(jobv, range, uplo, &n, rx, &n,
		     &vl, &vu, &il, &iu, &abstol, &m, rvalues,
		     rz, &n, isuppz,
		     work, &lwork, iwork, &liwork, &info);

    if ( info != 0 )
    	return( ERR_LAPACK ); // #  "error code %d from Lapack routine dsyevr";

    return 0;
}


int dmvnorm_log( CFmMatrix* fmY, CFmVector* pMu, CFmMatrix* pCovar, CFmVector* pPv )
{
	static CFmVector fmEigen( 0, 0.0 );

	int nRet = get_eigen( pCovar, true, &fmEigen, NULL );

	for( int i=0;i< fmEigen.GetLength(); i++ )
		fmEigen.Set( i, log( fmEigen.Get(i) ) );

    double distval = fmEigen.Sum();

	static CFmMatrix fmCovarInv( 0, 0 );
	fmCovarInv = pCovar->GetInverted();

	pPv->Resize(0);
	static CFmVector fmTmp( 0, 0.0 );
	for(int i=0; i< fmY->GetNumRows(); i++)
	{
		fmTmp = fmY->GetRow(i) - *pMu;

	    fmTmp = (fmTmp * fmCovarInv).GetRow(0) * fmTmp;

	    double logretval = - ( pCovar->GetNumCols() * log( 2 * M_PI ) + fmTmp.Sum() + distval )/2;

		pPv->Put( logretval );
	}

	return(0);
}