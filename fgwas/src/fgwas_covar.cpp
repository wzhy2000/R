#include <ctype.h>
#include <string.h>

#include "fm_matrix.h"
#include "fm_vector_str.h"
#include "fm_vector.h"
#include "fm_err.h"
#include "fm_new.h"

#include "fgwas_covar.h"


CFgCovariance::CFgCovariance(int nCovar)
{
	m_nCovar = nCovar;
}

CFgCovariance::~CFgCovariance()
{
}

int CFgCovariance::GetValue(CFmVector* pPar, CFmVector* pTime, CFmVector* pValue)
{
	return (0);
}

CFgCovariance_AR1::CFgCovariance_AR1(int nCovar):CFgCovariance(nCovar)
{
}

CFgCovariance_AR1::~CFgCovariance_AR1()
{
}

int CFgCovariance_AR1::GetValue(CFmVector* pPar, CFmVector* pTime, CFmMatrix* pValue)
{
	double rho = pPar->Get(0);
	double sigma = pPar->Get(1);
	int n= pTime->GetLength();
	pValue->Resize( n, n );

	for(int i=0; i<n; i++)
	for(int j=0; j<n; j++)
	{
		pValue->Set( i, j, pow( rho, abs(i-j) )* sigma * sigma );
	}

	return( 0 );

}

CFgCovariance* CreateCovariance(int nCovar)
{
	CFmNewTemp  fmRef;

	if (nCovar==1) return( new(fmRef) CFgCovariance_AR1(1));
	return NULL;
}

void destroy( CFgCovariance* p)
{
	CFmNewTemp  fmRef;
	p->~CFgCovariance();
	operator delete(p, fmRef);
}

