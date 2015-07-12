/* bls_RI.cpp  -	BLS application
 *	Copyright (C) 2011 THe Center for Statistical Genetics
 *  http://statgen.psu.edu
 */

#include <R.h>
#include <Rinternals.h>
#include <Rembedded.h>
#include <Rdefines.h>

#include "fgwas_R.h"

#define BOOLEAN_ELT(x,__i__)	LOGICAL(x)[__i__]
#define INTEGER_ELT(x,__i__)	INTEGER(x)[__i__]
#define NUMERIC_ELT(x,__i__)	REAL(x)[__i__]

SEXP SnpMle( SEXP spYMat,
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
	SEXP ret = _SnpMle_C( spYMat, spXMat, spSnpVec, spTimeVec, spCovarType, spCurveType, spXInit, spCovarInit, spCurveInit, spCurveExp, spOptimMethod );

	return(ret);
}
