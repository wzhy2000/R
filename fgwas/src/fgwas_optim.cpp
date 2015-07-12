/*
==> Nelder Mead:
void nmmin(int n, double *xin, double *x, double *Fmin, optimfn fn,
	int *fail, double abstol, double intol, void *ex,
	double alpha, double beta, double gamma, int trace,
	int *fncount, int maxit);

fail - receives true if the function failed
intol - user-initialized conversion tolerance
ex - data to pass to the optimization function (fn)
alpha - reflection factor
beta - contraction and reduction factor
gamma - extension factor
fncount - receives the number of times the optimization function was called in the iteration loop

==> BFGS:
void vmmin(int n, double *x, double *Fmin,
	optimfn fn, optimgr gr, int maxit, int trace,
	int *mask, double abstol, double reltol, int nREPORT,
	void *ex, int *fncount, int *grcount, int *fail);

==> Conjugate gradients:
void cgmin(int n, double *xin, double *x, double *Fmin,
	optimfn fn, optimgr gr, int *fail, double abstol,
	double intol, void *ex, int type, int trace,
	int *fncount, int *grcount, int maxit);

==> Limited-memory BFGS with bounds:
void lbfgsb(int n, int lmm, double *x, double *lower,
	double *upper, int *nbd, double *Fmin, optimfn fn,
	optimgr gr, int *fail, void *ex, double factr,
	double pgtol, int *fncount, int *grcount,
	int maxit, char *msg, int trace, int nREPORT);

==>  Simulated annealing:
	void samin(int n, double *x, double *Fmin, optimfn fn, int maxit,
	int tmax, double temp, int trace, void *ex);
*/

/* fgwas_R.cpp  -	fGWAS MLE estimate
 *	Copyright (C) 2014
 */

#include <R.h>
#include <Rinternals.h>
#include <Rembedded.h>
#include <Rdefines.h>
#include <R_ext/Applic.h>

#include "fm_rlogger.h"
#include "fm_matrix.h"
#include "fm_vector.h"
#include "fm_err.h"

#include "fgwas_optim.h"
#include "fgwas_curve.h"
#include "fgwas_covar.h"
#include "fgwas_dmvnorm.h"

CFgOptim::CFgOptim(CFmMatrix* pY, CFmMatrix* pX, CFmVector* pSNP, CFmVector* pTime, int nCovar, int nCurve, int nOptimType)
{
	m_pY = pY;
	m_pX = pX;
	m_pSNP = pSNP;
	m_pTime = pTime;
	m_pCovar = CreateCovariance(nCovar);

	m_pQQCurve = NULL;
	m_pQqCurve = NULL;
	m_pqqCurve = NULL;

	if ( pSNP->Find( 0.0 ) >= 0 ) m_pQQCurve = CreateCurve( nCurve );
	if ( pSNP->Find( 1.0 ) >= 0 ) m_pQqCurve = CreateCurve( nCurve );
	if ( pSNP->Find( 2.0 ) >= 0 ) m_pqqCurve = CreateCurve( nCurve );

	m_opCurve = new CFmMatrix( pSNP->GetLength(), pTime->GetLength() );

}

CFgOptim::~CFgOptim()
{
	if(m_pQQCurve) destroy(m_pQQCurve);
	if(m_pQqCurve) destroy(m_pQqCurve);
	if(m_pqqCurve) destroy(m_pqqCurve);
	if(m_pCovar) destroy(m_pCovar);
	if(m_opCurve) destroy(m_opCurve);
}

double fn_mle(int n, double *pPar, void* pExData)
{
	static CFmMatrix fmCovar( 0, 0 );
	static CFmMatrix fmY( 0, 0 );
	static CFmMatrix fmTmp( 0, 0 );
	static CFmVector fmXinit( 0, 0.0 );
	static CFmVector fQQVal( 0, 0.0 );
	static CFmVector fQqVal( 0, 0.0 );
	static CFmVector fqqVal( 0, 0.0 );
	static CFmVector fmPVec( 0, 0.0 );

	CFgOptim* pO = (CFgOptim*)pExData;

	int k0 = pO->m_pCovar->IntialParam( pPar+0 );

	int k = 0;
	if ( pO->m_pQQCurve ) k = pO->m_pQQCurve->IntialParam( pPar + k0 );
	if ( pO->m_pQqCurve ) k = pO->m_pQqCurve->IntialParam( pPar + k0 );
	if ( pO->m_pqqCurve ) k = pO->m_pqqCurve->IntialParam( pPar + k0 );

	k0 += k;

	fmXinit.Resize(0);
	for(int i=0; i< pO->m_pX->GetNumCols(); i++)
		fmXinit.Put(  pPar[ k0 + i] );

	fmY = *(pO->m_pY);

	fmTmp.Resize(fmY.GetNumRows(), fmY.GetNumCols() );

	for(int i=0; i< pO->m_pX->GetNumCols(); i++)
	{
		for(int j=0; j< fmY.GetNumCols(); j++ )
			fmTmp.SetCol( j, pO->m_pX->GetCol(i)*fmXinit[i] );

		fmY = fmY - fmTmp;
	}

	if( pO->m_pQQCurve ) fQQVal = pO->m_pQQCurve->GetValue( &fQQVal, pO->m_pTime );
	if( pO->m_pQqCurve ) fQqVal = pO->m_pQqCurve->GetValue( &fQqVal, pO->m_pTime );
	if( pO->m_pqqCurve ) fqqVal = pO->m_pqqCurve->GetValue( &fqqVal, pO->m_pTime );

	for(int i=0; i< pO->m_pSNP->GetLength();i++)
	{
		if( pO->m_pSNP->Get(i)==0.0 ) pO->m_opCurve->SetRow( i, fQQVal );
		if( pO->m_pSNP->Get(i)==1.0 ) pO->m_opCurve->SetRow( i, fQqVal );
		if( pO->m_pSNP->Get(i)==2.0 ) pO->m_opCurve->SetRow( i, fqqVal );
	}

	fmY = fmY - *(pO->m_opCurve);

	pO->m_pCovar->GetValue( &fmCovar, pO->m_pTime );

	int nRet = dmvnorm_log( &fmY, NULL, &fmCovar, &fmPVec );

	return( -1*fmPVec.Sum() );
}

int CFgOptim::Call( CFmVector* pXInit, CFmVector* pCovarInit, CFmVector* pCurveInit, CFmVector* pFmPar, double* fLoglike )
{
  	/*
  	 * The following values are based on the help
  	 * page for optim.
  	 */
  	const double abstol = 1e-16;
  	const double reltol = 1e-8;
  	const double alpha = 1.0; /* reflection factor */
  	const double beta = 0.5; /* contraction factor */
  	const double gamm = 2.0; /* expansion factor */
  	const int trace = 0; /* tracing on */
  	int fncount;
  	const int maxit = 10000;
  	double value;
  	int convergenceCode;

	CFmVector fmParInit(0, 0.0);
	fmParInit.Append( *pCovarInit );
	fmParInit.Append( *pCurveInit );
	fmParInit.Append( *pXInit );

	pFmPar->Resize( fmParInit.GetLength() );

  	nmmin(1, fmParInit.GetData(),
  		pFmPar->GetData(),
  		&value,
  		fn_mle,
		&convergenceCode,
		abstol, reltol,
		(void*)this,
		alpha, beta, gamm,trace,
		&fncount,
		maxit);

	*fLoglike = value;

  	Rprintf("fncount: %d\n", fncount);
  	Rprintf("convergence code: %d\n", convergenceCode);

  	return 0;
}


