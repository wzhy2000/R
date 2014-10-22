#include "fm_optim.h"
#include "fm_vector.h"
#include "fm_matrix.h"
#include "fm_rlogger.h"

#define stepredn	0.2
#define acctol		0.0001
#define reltest		10.0

CFmVector* CFmOptim::g_pMask = NULL;
CFmVector *CFmOptim::g_parscale = NULL;
CFmVector* CFmOptim::g_pDeps = NULL;
CFmVector* CFmOptim::g_pTmpPar = NULL;

CFmOptim::CFmOptim()
{
    m_maxit = 100;
    m_trace = 0;
    m_abstol = 1e-8;
    m_reltol = 1e-8;
    m_nREPORT = 10;
    m_fnscale = 1;
    m_usebounds = 0;
    m_lower = NULL;
    m_upper = NULL;
    m_fncount=0,
    m_grcount=0,
    m_ifail=0;

    if (!g_pMask)  g_pMask = new CFmVector(250, 1.0);
    if (!g_parscale)  g_parscale = new CFmVector(250, 1.0);
    if (!g_pDeps)  g_pDeps = new CFmVector(250, 1e-3);
    if (!g_pTmpPar)  g_pTmpPar = new CFmVector(250, 0);
}

CFmOptim::~CFmOptim()
{
}

int CFmOptim::do_BFGS( CFmVector& vctParin,
                       double& fmin,
                       void** pFnPar, int nFnPar,
                       pOptim_fn pFn,
                       pOptim_gr pGr )
{
    int npar = vctParin.GetLength();
    {
        g_pMask->Resize(npar, false);
        for(int n=0;n<npar;n++)g_pMask->Set(n, 1);

        g_parscale->Resize(npar, false);
        for(int n=0;n<npar;n++)g_parscale->Set(n, 1);

        g_pTmpPar->Resize(npar, false);
        for(int n=0;n<npar;n++)g_pTmpPar->Set(n, 0);

    }

    if (!pGr)
    {
        g_pDeps->Resize(npar, false);
        for(int n=0;n<npar;n++)g_pDeps->Set(n, 1e-3);
    }

    m_fncount=0,
    m_grcount=0,
    m_ifail=0;
    try
    {
        vmmin(vctParin, &fmin, pFnPar, nFnPar, pFn, pGr, &m_fncount, &m_grcount, &m_ifail);
    }
    catch(char* szError)
    {
        m_ifail = 3;
    }

    for (int i = 0; i < vctParin.GetLength(); i++)
        vctParin[i] = vctParin[i] * (*g_parscale)[i];

    fmin = fmin * m_fnscale;
    return(m_ifail);
}

/*  BFGS variable-metric method, based on Pascal code
in J.C. Nash, `Compact Numerical Methods for Computers', 2nd edition,
converted by p2c then re-crafted by B.D. Ripley */
void CFmOptim::vmmin(CFmVector& parin,
           double *Fmin,
           void** pFnPar,
           int nFnPar,
           pOptim_fn pFn,
           pOptim_gr pGr,
           int *fncount,
           int *grcount,
           int *fail)
{
    int n0 = parin.GetLength();
    double* b= parin.GetData();

    bool accpoint, enough;
    int   *l, j, ilast, iter = 0, count, funcount, gradcount;;
    double gradproj, s, steplength, D1, D2;

    if (m_maxit <= 0)
    {
        *fail = 0;
        *Fmin = fminfn( parin, pFn, pFnPar, nFnPar);
        *fncount = *grcount = 0;
        return;
    }

    if (m_nREPORT <= 0)
    {
        strcpy( m_szError, "REPORT must be > 0 (method = \"BFGS\")" );
        throw( m_szError );
    }

    l = new int[n0];
    int n = 0;
    for (int i = 0; i < n0; i++)
        if ((*g_pMask)[i])
            l[n++] = i;

    CFmVector vctDf(n0, 0.0);
    CFmVector t(n, 0.0);
    CFmVector X(n, 0.0);
    CFmVector c(n, 0.0);
    CFmMatrix B(n, n);

    double f = fminfn( parin, pFn, pFnPar, nFnPar);
    if (!R_FINITE(f))
    {
        strcpy( m_szError, "initial value in 'vmmin' is not finite" );
        throw(m_szError);
    }

    if (m_trace)
        Rprintf("initial  value %f \n", f);

    *Fmin = f;
    funcount = gradcount = 1;
    fmingr( parin, vctDf, pFn, pGr, pFnPar, nFnPar);
    iter++;
    ilast = gradcount;

    do
    {
        if (ilast == gradcount) {
            for (int i = 0; i < n; i++) {
                for (j = 0; j < i; j++)
                    B.Set(i,j, 0.0);
                B.Set(i,i,1.0);
            }
        }

        for (int i = 0; i < n; i++) {
            X[i] = b[l[i]];
            c[i] = vctDf[l[i]];
        }

        gradproj = 0.0;
        for (int i = 0; i < n; i++) {
            s = 0.0;
            for (int j = 0; j <= i; j++)
                s -= B.Get(i,j) * vctDf[l[j]];
            for (int j = i + 1; j < n; j++)
                s -= B.Get(j,i) * vctDf[l[j]];
            t[i] = s;
            gradproj += s * vctDf[l[i]];
        }

        if (gradproj < 0.0) {	/* search direction is downhill */
            steplength = 1.0;
            accpoint = FALSE;
            do {
                count = 0;
                for (int i = 0; i < n; i++) {
                    b[l[i]] = X[i] + steplength * t[i];
                    if ( reltest + X[i] == reltest + b[l[i]]) /* no change */
                    count++;
                }
                if (count < n)
                {
                    f = fminfn( parin, pFn, pFnPar, nFnPar);
                    funcount++;
                    accpoint = R_FINITE(f) && (f <= *Fmin + gradproj * steplength * acctol);
                    if (!accpoint)
                    {
                        steplength *= stepredn;
                    }
                }
            } while (!(count == n || accpoint));


            enough = (f > m_abstol) && fabs(f - *Fmin) > m_reltol * (fabs(*Fmin) + m_reltol);
            /* stop if value if small or if relative change is low */

            if (!enough)
            {
                count = n;
                *Fmin = f;
            }

            if (count < n) {/* making progress */
                *Fmin = f;
                fmingr( parin, vctDf, pFn, pGr, pFnPar, nFnPar );
                gradcount++;
                iter++;
                D1 = 0.0;
                for (int i = 0; i < n; i++)
                {
                    t[i] = steplength * t[i];
                    c[i] = vctDf[l[i]] - c[i];
                    D1 += t[i] * c[i];
                }

                if (D1 > 0) {
                    D2 = 0.0;

                    for (int i = 0; i < n; i++)
                    {
                        s = 0.0;
                        for (int j = 0; j <= i; j++)
                            s += B.Get( i, j ) * c[j];
                        for (int j = i + 1; j < n; j++)
                            s += B.Get( j, i ) * c[j];
                        X[i] = s;
                        D2 += s * c[i];
                    }

                    D2 = 1.0 + D2 / D1;
                    for (int i = 0; i < n; i++) {
                        for (int j = 0; j <= i; j++)
                            B.Set( i, j , B.Get( i, j ) + (D2 * t[i] * t[j] - X[i] * t[j] - t[i] * X[j]) / D1 );
                    }
                }
                else
                {	/* D1 < 0 */
                    ilast = gradcount;
                }
            }
            else
            {	/* no progress */
                if (ilast < gradcount) {
                    count = 0;
                    ilast = gradcount;
                }
            }
        }
        else
        {		/* uphill search */
            count = 0;
            if (ilast == gradcount)
                count = n;
            else
                ilast = gradcount;
            /* Resets unless has just been reset */
        }

        if (m_trace && (iter % m_nREPORT == 0))
            Rprintf("iter%4d value %f\n", iter, f);

        if (iter >= m_maxit) break;

        if (gradcount - ilast > 2 * n)
            ilast = gradcount;	/* periodic restart */

    } while (count != n || ilast != gradcount);

    if (m_trace)
    {
        Rprintf("final  value %f \n", *Fmin);
        if (iter < m_maxit)
            Rprintf("converged\n");
        else
            Rprintf("stopped after %i iterations\n", iter);
    }

    *fail = (iter < m_maxit) ? 0 : 2;
    *fncount = funcount;
    *grcount = gradcount;

    delete[] l;
}

double CFmOptim::fminfn(CFmVector& parin, pOptim_fn pFn, void** pFnPar, int nFnPar )
{
    for (int i = 0; i < parin.GetLength(); i++)
    {
        if (!R_FINITE( parin[i] ) )
        {
            strcpy(m_szError, "non-finite value supplied by optim");
            throw(m_szError);
        }
        (*g_pTmpPar)[i] = parin[i] * (*g_parscale)[i];
    }

    double s = (*pFn)( g_pTmpPar, pFnPar, nFnPar );
    double val = s/m_fnscale;
    return val;
}

void CFmOptim::fmingr(CFmVector& parin, CFmVector& vctDf, pOptim_fn pFn, pOptim_gr pGr, void** pFnPar, int nFnPar)
{
    int nPar = parin.GetLength();

    if ( pGr )
    {
        /* analytical derivatives */
        for (int i = 0; i < nPar; i++)
        {
            if ( !R_FINITE(parin[i]) )
            {
                strcpy( m_szError, "non-finite value supplied by optim");
                throw(  m_szError );
            }

            (*g_pTmpPar)[i] = parin[i] * (*g_parscale)[i];
        }

        CFmVector s(0, 0.0);
        (*pGr)( g_pTmpPar, &s, pFnPar, nFnPar );
        if( s.GetLength() != nPar)
        {
            sprintf( m_szError, "gradient in optim evaluated to length %d not %d", s.GetLength(), nPar);
            throw(m_szError);
        }

        for (int i = 0; i < nPar; i++)
            vctDf[i] = s[i] * (*g_parscale)[i]/m_fnscale;

    }
    else
    {	 /* numerical derivatives */
        for (int i = 0; i < nPar; i++)
            (*g_pTmpPar)[i] = parin[i] * (*g_parscale)[i];

        if( m_usebounds == 0)
        {
            for (int i = 0; i < nPar; i++)
            {
                int try_count=0;
                while(try_count<=3)
                {
                    double eps = (*g_pDeps)[i];

                    (*g_pTmpPar)[i] = (parin[i] + eps) * ((*g_parscale)[i]);
                    double s = (*pFn)(g_pTmpPar, pFnPar, nFnPar);

                    double val1 = s/(m_fnscale);
                    (*g_pTmpPar)[i] = (parin[i] - eps) * ((*g_parscale)[i]);
                    s = (*pFn)( g_pTmpPar, pFnPar, nFnPar);

                    double val2 = s/(m_fnscale);
                    double df = (val1 - val2)/(2 * eps);

                    if (R_FINITE(df))
                    {
                        vctDf[i] = df;
                        break;
                    }

                    try_count++;
                }

                if (! R_FINITE( vctDf[i] ) )
                {
                    sprintf( m_szError, ("non-finite finite-difference value [%d]"), i+1);
                    throw( m_szError );
                }

                (*g_pTmpPar)[i] = parin[i] * (*g_parscale)[i];
            }
        }
        else
        { /* usebounds */
            for (int i = 0; i < nPar; i++)
            {
                double eps = (*g_pDeps)[i];
                double epsused = eps;
                double tmp = parin[i] + eps;
                if (tmp > m_upper[i])
                {
                    tmp = m_upper[i];
                    epsused = tmp - parin[i] ;
                }

                (*g_pTmpPar)[i] = tmp * (*g_parscale)[i];

                double s = (*pFn)(g_pTmpPar, pFnPar, nFnPar);
                double val1 = s/m_fnscale;

                tmp = parin[i] - eps;
                if (tmp < m_lower[i])
                {
                    tmp = m_lower[i];
                    eps = parin[i] - tmp;
                }

                (*g_pTmpPar)[i] = tmp * (*g_parscale)[i];
                s = (*pFn)(g_pTmpPar, pFnPar, nFnPar);
                double val2 = s/(m_fnscale);
                vctDf[i] = (val1 - val2)/(epsused + eps);

                //if (!R_FINITE( vctDf[i] ))
                //    error("non-finite finite-difference value [%d]", i+1);

                (*g_pTmpPar)[i] = parin[i] * (*g_parscale)[i];
            }
        }
    }
}
