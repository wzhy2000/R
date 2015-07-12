#ifndef _FG_OPTIM_H_
#define _FG_OPTIM_H_

class CFmVectorStr;
class CFmVector;
class CFmMatrix;
class CFgCurve; 
class CFgCovariance;

class CFgOptim
{
public:
    CFgOptim(CFmMatrix* pY, CFmMatrix* pX, CFmVector* pSNP, CFmVector* pTime, int nCovar, int nCurve, int nOptimType);
    virtual ~CFgOptim();

    int Call( CFmVector* pXInit, CFmVector* pCovarInit, CFmVector* pCurveInit, CFmVector* pFmPar, double* fLoglike);

    CFmMatrix* m_pY;
    CFmMatrix* m_pX;
    CFmVector* m_pSNP;
    CFmVector* m_pTime;
    CFgCovariance* m_pCovar;

    CFgCurve* m_pQQCurve;
    CFgCurve* m_pQqCurve;
    CFgCurve* m_pqqCurve;

    CFmMatrix* m_opCurve;

private:

};

void destroy( CFgOptim* p);

#endif // _FG_OPTIM_H_

