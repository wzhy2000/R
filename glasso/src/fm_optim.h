#ifndef _FM_OPTIM_H_
#define _FM_OPTIM_H_

class CFmVector;

typedef double (*pOptim_fn)(CFmVector* vctParin, void** pFnPar, int nFnPar);
typedef int (*pOptim_gr)(CFmVector* vctParin, CFmVector* vctDf, void** pFnPar, int nFnPar);

class CFmOptim
{
public:
    CFmOptim();
    virtual ~CFmOptim();

    int do_BFGS( CFmVector& vctParin, double& fVal, void** pFnPar, int nFnPar, pOptim_fn fn, pOptim_gr gr );

    int     m_maxit;
    int     m_trace;
    double  m_abstol;
    double  m_reltol;
    int     m_nREPORT;
    double  m_fnscale;
    int     m_usebounds;
    double* m_lower, *m_upper;
    //A vector of step sizes for the finite-difference approximation to the gradient, on par/parscale scale. Defaults to 1e-3.
    char    m_szError[256];
    int     m_fncount;
    int     m_grcount;
    int     m_ifail;

private:
    static CFmVector *g_pMask;
    static CFmVector *g_parscale;
    static CFmVector* g_pDeps;
    static CFmVector* g_pTmpPar;

    void vmmin(CFmVector& parin, double *Fmin,void** pFnPar, int nFnPar, pOptim_fn pFn, pOptim_gr pGr, int *fncount,int *grcount, int *fail);
    void fmingr(CFmVector& parin, CFmVector& vctDf, pOptim_fn pFn, pOptim_gr pGr, void** pFnPar, int nFnPar);
    double fminfn(CFmVector& parin, pOptim_fn pFn, void** pFnPar, int nFnPar );

};

#endif // _FM_OPTIM_H_
