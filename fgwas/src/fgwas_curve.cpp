#include <ctype.h>
#include <string.h>

#include "fm_matrix.h"
#include "fm_vector_str.h"
#include "fm_vector.h"
#include "fm_err.h"
#include "fm_new.h"

#include "fgwas_curve.h"


CFgCurve::CFgCurve(int nCurve)
{
	m_nCurve = nCurve;
}

CFgCurve::~CFgCurve()
{
}

int CFgCurve::GetValue(CFmVector* pPar, CFmVector* pTime, CFmVector* pValue)
{
	return (0);
}

CFgCurve_Log::CFgCurve_Log(int nCurve):CFgCurve(nCurve)
{
}

CFgCurve_Log::~CFgCurve_Log()
{
}

int CFgCurve_Log::GetValue(CFmVector* pPar, CFmVector* pTime, CFmVector* pValue)
{
	double a = pPar->Get(0);
	double b = pPar->Get(1);
	double r = pPar->Get(2);

	pValue->Resize(0);
	for(int i=0;i<pTime->GetLength();i++)
	{
		double f = a/(1+b*exp(-r * pTime->Get(i) ) );
		pValue->Put( f );
	}

	return( 0 );

}

CFgCurve* CreateCurve(int nCurve)
{
	CFmNewTemp  fmRef;

	if (nCurve==1) return( new(fmRef) CFgCurve_Log(1));
	return NULL;
}

void destroy( CFgCurve* p)
{
	CFmNewTemp  fmRef;
	p->~CFgCurve();
	operator delete(p, fmRef);
}

