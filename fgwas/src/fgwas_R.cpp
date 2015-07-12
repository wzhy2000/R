/* fgwas_R.cpp  -	fGWAS MLE estimate
 *	Copyright (C) 2014
 */

#include <R.h>
#include <Rinternals.h>
#include <Rembedded.h>
#include <Rdefines.h>

#include "fm_rlogger.h"
#include "fm_matrix.h"
#include "fm_vector.h"
#include "fm_err.h"

#include "fgwas_R.h"
#include "fgwas_optim.h"



#define BOOLEAN_ELT(x,__i__)	LOGICAL(x)[__i__]
#define INTEGER_ELT(x,__i__)	INTEGER(x)[__i__]
#define NUMERIC_ELT(x,__i__)	REAL(x)[__i__]
#define STRING_ELT(x,__i__)	    CHAR(STRING_ELT(x,__i__))

SEXP X_SnpMle_C( CFmMatrix* pFmYMat, CFmMatrix* pFmXMat, CFmVector* pFmSNP, CFmVector* pFmTime,
			CFmVector* pXInit, CFmVector* pCovarInit, CFmVector* pCurveInit,
			int nCovarType, int nCurveType, int nOptimType, const char* szCurveExp)
{

	CFgOptim fgOptim(pFmYMat, pFmXMat, pFmSNP, pFmTime, nCovarType, nCurveType, nOptimType);

	CFmVector fmPar( 0, 0.0 );
	double fLoglike = 0;

	int nRet = fgOptim.Call(pXInit, pCovarInit, pCurveInit, &fmPar, &fLoglike);

	SEXP sRet, t;

   	PROTECT(sRet = t = allocList(2));

	SEXP expVS = GetSEXP(&fmPar);
	SETCAR( t, expVS );
	SET_TAG(t, install("par") );
	t = CDR(t);

	CFmVector fmLogL(1, fLoglike);
	SEXP expVS1 = GetSEXP(&fmLogL);
	SETCAR( t, expVS1 );
	SET_TAG(t, install("LogLik") );
	t = CDR(t);

	UNPROTECT(1);

    return(sRet);
}

SEXP _SnpMle_C( SEXP spYMat,
			SEXP spXMat,
			SEXP spSnpVec,
			SEXP spTimeVec,
  		   	SEXP spCovarType,
  		   	SEXP spCurveType,
  		   	SEXP spXInit,
  		   	SEXP spCovarInit,
  		   	SEXP spCurveInit,
  		   	SEXP spCurveExp,
  		   	SEXP spOptimMethod )
{
	//int nUsed0, nTotal0;
	//CFmVector::StatCache( &nTotal0, &nUsed0 );
	//int nUsed1, nTotal1;
	//CFmMatrix::StatCache( &nTotal1, &nUsed1 );
	//
	//Rprintf( "Enter C Range, Vec.count=%d, Mat.count=%d\n", nUsed0, nUsed1);

	CFmMatrix fmYMat( 0, 0 );
	CFmMatrix fmXMat( 0, 0 );
	CFmVector fmSNP( 0, 0.0 );
	CFmVector fmTime( 0, 0.0 );
	CFmVector fmXInit( 0, 0.0 );
	CFmVector fmCovarInit( 0, 0.0 );
	CFmVector fmCurveInit( 0, 0.0 );

	GetMatrix( spYMat, &fmYMat );
	GetMatrix( spXMat, &fmXMat );
	GetVector( spSnpVec, &fmSNP );
	GetVector( spTimeVec, &fmTime );
	GetVector( spXInit, &fmXInit );
	GetVector( spCovarInit, &fmCovarInit );
	GetVector( spCurveInit, &fmCurveInit );

	int nCovarType   = INTEGER_ELT( spCovarType, 0 );
	int nCurveType   = INTEGER_ELT( spCurveType, 0 );
	int nOptimType   = INTEGER_ELT( spOptimMethod, 0 );
	const char* szCurveExp = STRING_ELT( spCurveExp, 0 );

	SEXP ret;
	try
	{
		ret = X_SnpMle_C( &fmYMat, &fmXMat, &fmSNP, &fmTime, &fmXInit, &fmCovarInit, &fmCurveInit, nCovarType, nCurveType, nOptimType, szCurveExp);
	}
    catch(const char* str)
    {
        _log_error( _HI_, "Exception=%s", str);
        return( R_NilValue );
    }

	//CFmVector::StatCache( &nTotal0, &nUsed0 );
	//CFmMatrix::StatCache( &nTotal1, &nUsed1 );
	//Rprintf( "Leave C Range, Vec.count=%d, Mat.count=%d\n", nUsed0, nUsed1);

	return(ret);
}


