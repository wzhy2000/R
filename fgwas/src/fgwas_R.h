// lskat_R.h: interface for the LSKAT.
//
//////////////////////////////////////////////////////////////////////

#if !defined(LSKAT_R_H__INCLUDED_)
#define LSKAT_R_H__INCLUDED_

#include <stdbool.h>
#include "fm_linux.h"

#ifdef __cplusplus
extern "C" {
#endif

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
  		   	SEXP spOptimMethod );

#ifdef __cplusplus
}
#endif

#endif


