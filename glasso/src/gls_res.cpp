/* GLS_res.cpp  -	BLS Result Object
 *
 *	Copyright (C) 2011 THe Center for Statistical Genetics
 *  http://statgen.psu.edu
 */

#include <stdio.h>
#include <R.h>
#include <Rinternals.h>
#include <Rembedded.h>
#include <Rdefines.h>
#include <Rmath.h>

#include "fm_matrix.h"
#include "fm_vector.h"
#include "fm_vector_str.h"
#include "fm_rlogger.h"
#include "fm_err.h"
#include "fm_rls.h"
#include "fm_filematrix.h"
#include "gls_cfg.h"
#include "gls_res.h"

#define _t(x) ((x).GetTransposed())

GLS_res::GLS_res(CMDOPTIONS *pCmd, GLS_cfg* pCfg)
{
    m_pCfg = pCfg;
    m_pCmd = pCmd;

    m_nSimuRound = 1;
    m_nSnpP = 0;
    m_nSubjN = 0;
    m_nMcmcIter = 0;
    m_nRefitSnp = 0;

    m_pVarsel_SnpName = NULL;
    m_pVarsel_SnpChr  = NULL;
    m_pVarsel_SnpPos  = NULL;
    m_pVarsel_Mu      = NULL;
    m_pVarsel_Alpha   = NULL;
    m_pVarsel_Ra      = NULL;
    m_pVarsel_Rd      = NULL;
    m_pVarsel_BestQ   = NULL;

    m_pRefit_SnpName  = NULL;
    m_pRefit_SnpChr   = NULL;
    m_pRefit_SnpPos   = NULL;
    m_pRefit_Mu       = NULL;
    m_pRefit_Alpha    = NULL;
    m_pRefit_Ra       = NULL;
    m_pRefit_Rd       = NULL;
    m_pRefit_BestQ    = NULL;

    m_pRefitSnps  = NULL;
    m_pSigSnps	  = NULL;
    m_pSigAddSnps = NULL;
    m_pSigDomSnps = NULL;
    m_nTotalSnp = 0;
}

GLS_res::~GLS_res()
{
    if(m_pVarsel_SnpName) delete m_pVarsel_SnpName;
    if(m_pVarsel_SnpChr) delete m_pVarsel_SnpChr;
    if(m_pVarsel_SnpPos) delete m_pVarsel_SnpPos;
    if(m_pVarsel_Mu) delete m_pVarsel_Mu;
    if(m_pVarsel_Alpha) delete m_pVarsel_Alpha;
    if(m_pVarsel_Ra)   delete m_pVarsel_Ra;
    if(m_pVarsel_Rd)   delete m_pVarsel_Rd;
    if(m_pVarsel_BestQ) delete m_pVarsel_BestQ;

    if(m_pRefit_SnpName) delete m_pRefit_SnpName;
    if(m_pRefit_SnpChr) delete m_pRefit_SnpChr;
    if(m_pRefit_SnpPos) delete m_pRefit_SnpPos;
    if(m_pRefit_Mu) delete m_pRefit_Mu;
    if(m_pRefit_Alpha) delete m_pRefit_Alpha;
    if(m_pRefit_Ra) delete m_pRefit_Ra;
    if(m_pRefit_Rd) delete m_pRefit_Rd;
    if(m_pRefit_BestQ) delete m_pRefit_BestQ;

	if(m_pRefitSnps) delete m_pRefitSnps;
    if(m_pSigSnps) delete m_pSigSnps;
    if(m_pSigAddSnps) delete m_pSigAddSnps;
    if(m_pSigDomSnps) delete m_pSigDomSnps;


    _log_debug(_HI_, "GLS_res is released successfully.");
}


// mRound:   Simulation Round,
// nSnpP,    SNP count,
// nSujN,    Individual count,
// mcmcIter  MCMC maximum iteration
int GLS_res::InitVarsel(int mRound, int nSnpP, int SubjN, int mcmcIter)
{
    m_nSimuRound = mRound;
    m_nSnpP      = nSnpP;
    m_nSubjN     = SubjN;
    m_nMcmcIter  = mcmcIter;

    return(0);
}

char* GLS_res::GetVarselSnpName(int idx)
{
    if (idx<0 || idx>=m_pVarsel_SnpName->GetLength())
        return(NULL);

    return (m_pVarsel_SnpName->Get(idx));
}

int GLS_res::GetRefitSnp(CFmVector** pVct)
{
    if (m_pRefitSnps)
    {
        *pVct = m_pRefitSnps;
        return(0);
    }
    else
        return(ERR_NULL_DATA);
}

int GLS_res::SetMcmcResults(bool bRefit, CFmVectorStr* pVctSnpName, CFmVector* pVctChr, CFmVector* pVctPos, CFmFileMatrix* pMatRet, int nCov )
{
    _log_debug(_HI_, "SetMcmcResults: bRefit=%d, pMat.Dim(%d,%d)", bRefit?1:0, pMatRet->GetNumRows(), pMatRet->GetNumCols());

    if (bRefit)
    {
		int nSnp = m_nRefitSnp;

        m_pRefit_SnpName = new CFmVectorStr(pVctSnpName);
        m_pRefit_SnpChr = new CFmVector(pVctChr);
        m_pRefit_SnpPos = new CFmVector(pVctPos);

        m_pRefit_Mu = new CFmMatrix(1, LG*4+3);
        m_pRefit_Alpha = new CFmMatrix(nCov-1, LG*4+3);
        m_pRefit_Ra = new CFmMatrix(nSnp, LG*4+3);
        m_pRefit_Rd = new CFmMatrix(nSnp, LG*4+3);

		SortoutMcmc(pMatRet, m_pRefit_Mu, m_pRefit_Alpha, m_pRefit_Ra, m_pRefit_Rd, nSnp, nCov);

		m_pRefit_Ra->SetRowNames(m_pRefit_SnpName);
		m_pRefit_Rd->SetRowNames(m_pRefit_SnpName);

		m_pRefit_BestQ = new CFmMatrix(nSnp, (LG*4+3)*2 );
		SortoutBestQ( pMatRet, m_pRefit_BestQ , nSnp, nCov );
		m_pRefit_BestQ->SetRowNames( m_pRefit_SnpName );
    }
    else
    {
		int nSnp = m_nSnpP;

        m_pVarsel_SnpName = new CFmVectorStr(pVctSnpName);
        m_pVarsel_SnpChr = new CFmVector(pVctChr);
        m_pVarsel_SnpPos = new CFmVector(pVctPos);

        m_pVarsel_Mu = new CFmMatrix(1, LG*4+3);
        m_pVarsel_Alpha = new CFmMatrix(nCov-1, LG*4+3);
        m_pVarsel_Ra = new CFmMatrix(nSnp, LG*4+3);
        m_pVarsel_Rd = new CFmMatrix(nSnp, LG*4+3);

		SortoutMcmc(pMatRet, m_pVarsel_Mu, m_pVarsel_Alpha, m_pVarsel_Ra, m_pVarsel_Rd, nSnp, nCov);

		m_pVarsel_Ra->SetRowNames(m_pVarsel_SnpName);
		m_pVarsel_Rd->SetRowNames(m_pVarsel_SnpName);

		m_pVarsel_BestQ = new CFmMatrix(nSnp, (LG*4+3)*2);
		SortoutBestQ( pMatRet, m_pVarsel_BestQ , nSnp, nCov );
		m_pVarsel_BestQ->SetRowNames( m_pVarsel_SnpName );
    }

    return(0);
}

int GLS_res::SortoutMcmc(CFmFileMatrix* pMatRet, CFmMatrix* pMu, CFmMatrix* pAlpha, CFmMatrix* pRa, CFmMatrix* pRd, int nSnp, int nCov)
{
	_log_debug(_HI_, "SortoutMcmc: AD matrix file(%s). SNP=%d nCov=%d q_add=%.4f q_dom=%.4f", pMatRet->GetFileName(), nSnp, nCov, m_pCfg->m_fQval_add, m_pCfg->m_fQval_dom );

	CFmVector fmTmp(0, 0.0);
	CFmVector fmModel(0, 0.0);

	GetMcmcInfo( pMatRet, 0, &fmTmp, &fmModel, m_pCfg->m_fQval_add );
	pMu->SetRow(0, fmTmp);

	for (long int i=0; i<nCov-1; i++)
    {
		GetMcmcInfo( pMatRet, i+1, &fmTmp, &fmModel, m_pCfg->m_fQval_add );
		pAlpha->SetRow(i, fmTmp);
	}

	for (long int i=0; i< nSnp; i++)
    {
		GetMcmcInfo( pMatRet, nCov+i*2, &fmTmp, &fmModel, m_pCfg->m_fQval_add );
		pRa->SetRow(i, fmTmp);

		if (fmModel.Sum()>0)
			_log_debug(_HI_, "significant Add [%d], %.4f", i+1, fmTmp[4]);

		GetMcmcInfo( pMatRet, nCov+i*2+1, &fmTmp, &fmModel, m_pCfg->m_fQval_dom );
		pRd->SetRow(i, fmTmp);

		if (fmModel.Sum()>0)
			_log_debug(_HI_, "significant Dom [%d], %.4f", i+1, fmTmp[4]);
    }

    _log_debug(_HI_, "SortoutMcmc: Successful to sort out.");
    return(0);
}

int GLS_res::GetMcmcInfo(CFmFileMatrix* pMatRet, int idx, CFmVector* pFmInfo, CFmVector* pFmModel, double fQval )
{
    int nMcmc= m_pCfg->m_nMaxIter - m_pCfg->GetBurnInRound();
    int q1 = (int)floor( fQval *nMcmc);
    int q2 = (int)floor(( 1 - fQval )*nMcmc -1 );

	pFmInfo->Resize(0);
	pFmModel->Resize(0);

	CFmVector fmVct(0, 0.0);
	CFmVector fmModel(0, 0.0);
    CFmVector fmQ(0, 0.0);
    CFmMatrix fmMat( q2-q1+1, LG);

    for(int k=0; k< LG ; k++)
    {
        int ret = pMatRet->GetCacheCol(idx*LG + k, fmVct);
        if(ret!=0) return(ret);

        fmVct.Sort();

		fmQ.Resize(0);
		for(int qx=q1; qx<=q2; qx++)
			fmQ.Put(fmVct[qx]);

		fmMat.SetCol(k, fmQ );
        if (fmVct[q1]<=0 && fmVct[q2]>=0)
            fmModel.Put( 0.0 );
		else
            fmModel.Put( 1.0 );
	}

	CFmVector fmMedian(0, 0.0);
	CFmVector fmMin(0, 0.0);
	CFmVector fmMax(0, 0.0);
	CFmVector fmTmp(0, 0.0);
    for(int k=0; k< LG ; k++)
    {
		fmTmp = fmMat.GetCol(k);
		fmMedian.Put(fmTmp.GetMedian());
		fmMin.Put(fmTmp.GetMin());
		fmMax.Put(fmTmp.GetMax());
	}

	pFmModel->Append(fmModel);
	pFmInfo->Append(fmModel);

	fmTmp = fmMedian * fmModel;
	fmTmp = fmTmp*fmTmp;
	pFmInfo->Put(fmTmp.Sum());
	pFmInfo->Append(fmMedian);

	fmTmp = fmMin * fmModel;
	fmTmp = fmTmp*fmTmp;
	pFmInfo->Put(fmTmp.Sum());
	pFmInfo->Append(fmMin);

	fmTmp = fmMax * fmModel;
	fmTmp = fmTmp*fmTmp;
	pFmInfo->Put(fmTmp.Sum());
	pFmInfo->Append(fmMax);

	return(0);
}

int GLS_res::InitRefit(int nMcmcIter)
{
    if (m_pRefitSnps) delete m_pRefitSnps;

    m_pRefitSnps = new CFmVector( 0, 0.0 );
    for(int i=0;i<m_nSnpP; i++)
    {
		int nSumA=0;
		for(int j=0;j<LG;j++)
			nSumA = nSumA + (int)m_pVarsel_Ra->Get(i, j);

		int nSumD=0;
		for(int j=0;j<LG;j++)
			nSumD = nSumD + (int)m_pVarsel_Rd->Get(i, j);

		if(nSumA+nSumD>0)
			m_pRefitSnps->Put(i);
	}

    if (m_pRefitSnps->GetLength()<=0)
    {
        _log_prompt( _HI_, "REFIT: No SNP is selected for the refit procedure");
        return(0);
    }

    m_nRefitSnp = m_pRefitSnps->GetLength();
    m_nMcmcIter = nMcmcIter;

    return(0);
}

int GLS_res::GetSigSNP()
{
    if (m_pSigSnps) delete m_pSigSnps;
    if (m_pSigAddSnps) delete m_pSigAddSnps;
    if (m_pSigDomSnps) delete m_pSigDomSnps;

    m_pSigSnps = new CFmVector( 0, 0.0 );
    m_pSigAddSnps = new CFmVector( 0, 0.0 );
    m_pSigDomSnps = new CFmVector( 0, 0.0 );

    for(int i=0;i<m_nRefitSnp; i++)
    {
		int nSumA=0;
		for(int j=0;j<LG;j++)
			nSumA = nSumA + (int)m_pRefit_Ra->Get(i, j);
		if(nSumA>0)
			m_pSigAddSnps->Put(i);

		int nSumD=0;
		for(int j=0;j<LG;j++)
			nSumD = nSumD + (int)m_pRefit_Rd->Get(i, j);
		if(nSumD>0)
			m_pSigDomSnps->Put(i);

		if(nSumA+nSumD>0)
			m_pSigSnps->Put(i);
	}

    if (m_pSigSnps->GetLength()<=0)
    {
        _log_prompt( _HI_, "GetsigSnp: No signicant SNP is selected!");
        return(0);
    }

    return(0);
}

int GLS_res::Summary(char* szOutFile)
{
    char szTemp[MAX_PATH*256] = {0};

    sprintf(szTemp, "%s \nSUMMARY(Result)\n", szTemp );
    sprintf(szTemp, "%s -------------------------------------------\n", szTemp );
    sprintf(szTemp, "%s RESULT Subjects    = %d\n",  szTemp, m_nSubjN);
    sprintf(szTemp, "%s RESULT Snps        = %d\n",  szTemp, m_nSnpP);
    sprintf(szTemp, "%s RESULT Refit level = %.3f\n",szTemp, m_pCfg->m_fDetRefit);
    sprintf(szTemp, "%s RESULT Refit SNP   = %d\n",  szTemp, m_nRefitSnp);
    for(int i=0; i<m_nRefitSnp && i<20 ; i++)
    {
        int idx = (int)m_pRefitSnps->Get(i);
        sprintf(szTemp, "%s RESULT (VARSEL) [%02d]= %s, %.2f, %.2f\n",
                szTemp, i+1, m_pVarsel_SnpName->Get(idx),
                m_pVarsel_Ra->Get(idx,LG + 1), m_pVarsel_Ra->Get(idx, LG + 1) );
    }

    CFmVector* pSigSnpA=m_pSigAddSnps;
    CFmVector* pSigSnpD=m_pSigDomSnps;

    sprintf(szTemp, "%s RESULT (REFIT) Sig. Additive SNP   = %d\n", szTemp, m_pSigAddSnps->GetLength());

    for(int i=0; i<pSigSnpA->GetLength() && i<20 ; i++)
    {
        int idx = (int)pSigSnpA->Get(i);
        sprintf(szTemp, "%s RESULT A[%02d]= %s, %.2f, %.2f\n",
                szTemp, i+1, m_pRefit_SnpName->Get(idx),
                m_pRefit_Ra->Get(idx, LG + 1), m_pRefit_Rd->Get(idx, LG + 1) );
    }

    sprintf(szTemp, "%s RESULT (REFIT) Sig. Dominant SNP   = %d\n", szTemp, m_pSigDomSnps->GetLength());

    for(int i=0; i<pSigSnpD->GetLength() && i<20 ; i++)
    {
        int idx = (int)pSigSnpD->Get(i);
        sprintf(szTemp, "%s RESULT D[%02d]= %s, %.2f, %.2f\n",
                szTemp, i+1, m_pRefit_SnpName->Get(idx),
                m_pRefit_Ra->Get(idx, LG + 1), m_pRefit_Rd->Get(idx, LG + 1) );
    }
    sprintf(szTemp, "%s -------------------------------------------\n", szTemp );

    _log_info(_HI_, szTemp);

    if(szOutFile!=NULL)
    {
        FILE* fp = fopen(szOutFile, "at");
        if (fp==NULL)
        {
            _log_error(_HI_, "Failed to open file(%s) to append the summary information", szOutFile);
        }
        else
        {
            fwrite(szTemp, strlen(szTemp),1, fp);
            fclose(fp);
        }
    }

    return(0);
}

int GLS_res::SaveRData( char* szRdataFile )
{
    _log_info(_HI_, "SaveRData: Start to save the RDATA to file(%s).", szRdataFile);

    CFmRls rls(szRdataFile);

	CFmMatrix matVs1(0,0);
	matVs1.Cbind( *m_pVarsel_SnpChr, (char*)"Chr");
	matVs1.Cbind( *m_pVarsel_SnpPos, (char*)"Pos");
	matVs1.Cbind( *m_pVarsel_Ra);
	matVs1.SetRowNames( m_pVarsel_SnpName);

    rls.SetData( (char*)"varsel.add", &matVs1);

    matVs1.Resize(0,0);
	matVs1.Cbind( *m_pVarsel_SnpChr, (char*)"Chr");
	matVs1.Cbind( *m_pVarsel_SnpPos, (char*)"Pos");
	matVs1.Cbind( *m_pVarsel_Rd);
	matVs1.SetRowNames(m_pVarsel_SnpName);
    rls.SetData( (char*)"varsel.dom", &matVs1);

    if (m_pVarsel_Mu)
	    rls.SetData( (char*)"varsel.mu", m_pVarsel_Mu);

    if (m_pVarsel_Alpha)
	    rls.SetData( (char*)"varsel.cov", m_pVarsel_Mu);

    if (m_pVarsel_BestQ)
		rls.SetData( (char*)"varsel.Qbest", m_pVarsel_BestQ);

    if (m_pRefit_SnpName)
    {
        if (m_pRefit_Ra)
        {
	        CFmMatrix matVs2(0,0);
	        matVs2.Cbind( *m_pRefit_SnpChr, (char*)"Chr");
	        matVs2.Cbind( *m_pRefit_SnpPos, (char*)"Pos");
	        matVs2.Cbind( *m_pRefit_Ra);
	        matVs2.SetRowNames(m_pRefit_SnpName);
			rls.SetData( (char*)"refit.add", &matVs2);
        }
        if (m_pRefit_Rd)
		{
	        CFmMatrix matVs2(0,0);
	        matVs2.Cbind( *m_pRefit_SnpChr, (char*)"Chr");
	        matVs2.Cbind( *m_pRefit_SnpPos, (char*)"Pos");
	        matVs2.Cbind( *m_pRefit_Rd);
	        matVs2.SetRowNames(m_pRefit_SnpName);
			rls.SetData( (char*)"refit.dom", &matVs2);
		}
        if (m_pRefit_Mu)
			rls.SetData( (char*)"refit.mu", m_pRefit_Mu);

        if (m_pRefit_Alpha)
			rls.SetData( (char*)"refit.cov", m_pRefit_Alpha);

        if (m_pRefit_BestQ)
			rls.SetData( (char*)"refit.Qbest", m_pRefit_BestQ);
    }

     int nRet = rls.SaveRData( szRdataFile );

    if(nRet==0)
        _log_debug(_HI_, "SaveRData: Successful to save the RDATA to file(%s).", szRdataFile);

    return(nRet);
}

SEXP GLS_res::GetRObj()
{
    _log_info(_HI_, "GetRObj: Start to save the matrix to R.");

	SEXP sRet, t;
	int nList = 4;
	if (m_pRefit_SnpName) nList = 8;

   	PROTECT(sRet = t = allocList(nList));

    // #NO1: ret.vs
    if (m_pVarsel_Ra)
    {
		CFmMatrix matVs1(0,0);
	    matVs1.Cbind( *m_pVarsel_SnpChr );
	    matVs1.Cbind( *m_pVarsel_SnpPos );
	    matVs1.Cbind( *m_pVarsel_Ra);
		matVs1.SetRowNames(m_pVarsel_SnpName);

		SEXP expVS = GetSEXP(&matVs1);
		SETCAR( t, expVS );
		SET_TAG(t, install("varsel_add") );
		t = CDR(t);
     }

    if (m_pVarsel_Rd)
    {
		CFmMatrix matVs1(0,0);
	    matVs1.Cbind( *m_pVarsel_SnpChr );
	    matVs1.Cbind( *m_pVarsel_SnpPos );
	    matVs1.Cbind( *m_pVarsel_Rd);
		matVs1.SetRowNames(m_pVarsel_SnpName);

		SEXP expVS = GetSEXP(&matVs1);
		SETCAR( t, expVS );
		SET_TAG(t, install("varsel_dom") );
		t = CDR(t);
	}

    if (m_pVarsel_Mu)
    {
		CFmMatrix matVs1(0,0);
		CFmVector fmTmp(0, 0.0);
		matVs1.Cbind(m_pVarsel_Mu->GetRow(0) );
        if (m_pVarsel_Alpha)
        {
			for(int k=0; k<m_pVarsel_Alpha->GetNumRows();k++)
			{
				fmTmp = m_pVarsel_Alpha->GetRow(k);
				matVs1.Cbind( fmTmp );
			}
        }

		matVs1.Transpose();
  		SEXP expCoeff = GetSEXP(&matVs1);
		SETCAR( t, expCoeff );
		SET_TAG(t, install("varsel_cov") );
		t = CDR(t);
	}

    if (m_pVarsel_BestQ)
    {
  		SEXP expQBest = GetSEXP(m_pVarsel_BestQ);
		SETCAR( t, expQBest );
		SET_TAG(t, install("varsel_Qbest") );
		t = CDR(t);
	}


    // #NO2: ret.refit
    if (m_pRefit_SnpName)
    {
        if (m_pRefit_Ra)
        {
			CFmMatrix matVs1(0,0);
			matVs1.Cbind( *m_pRefit_SnpChr );
			matVs1.Cbind( *m_pRefit_SnpPos );
			matVs1.Cbind( *m_pRefit_Ra);
			matVs1.SetRowNames(m_pRefit_SnpName);

			SEXP expVS = GetSEXP(&matVs1);
			SETCAR( t, expVS );
			SET_TAG(t, install("refit_add") );
			t = CDR(t);
        }

        if (m_pRefit_Rd)
        {
			CFmMatrix matVs1(0,0);
			matVs1.Cbind( *m_pRefit_SnpChr );
			matVs1.Cbind( *m_pRefit_SnpPos );
			matVs1.Cbind( *m_pRefit_Rd);
			matVs1.SetRowNames(m_pRefit_SnpName);

			SEXP expVS = GetSEXP(&matVs1);
			SETCAR( t, expVS );
			SET_TAG(t, install("refit_dom") );
			t = CDR(t);
        }

        if (m_pRefit_Mu)
        {
			CFmMatrix matVs1(0,0);
			CFmVector fmTmp(0, 0.0);
			matVs1.Cbind(m_pRefit_Mu->GetRow(0) );
	        if (m_pRefit_Alpha)
	        {
				for(int k=0; k<m_pRefit_Alpha->GetNumRows();k++)
				{
					fmTmp = m_pRefit_Alpha->GetRow(k);
					matVs1.Cbind( fmTmp );
				}
	        }
			matVs1.Transpose();

	  		SEXP expCoeff = GetSEXP(&matVs1);
			SETCAR( t, expCoeff );
			SET_TAG(t, install("refit_cov") );
			t = CDR(t);
		}

		if (m_pRefit_BestQ)
		{
			SEXP expQBest = GetSEXP(m_pRefit_BestQ);
			SETCAR( t, expQBest );
			SET_TAG(t, install("refit_Qbest") );
		}
    }

	UNPROTECT(1);

    _log_debug(_HI_, "GetRObj: Successful to save the matrix to R.");

    return( sRet );
}

int GLS_res::SaveCsv4Fig( char* szCsvFile, bool bRefit )
{
    _log_debug( _HI_, "SaveCsv4Fig varsel(%d), refit(%d)",
                (m_pVarsel_SnpName?m_pVarsel_SnpName->GetLength():0) ,
                (m_pRefit_SnpName?m_pRefit_SnpName->GetLength():0) );

    CFmMatrix fm(0,0);
    if(!bRefit)
    {
        fm.Cbind( *m_pVarsel_SnpChr, (char*)"Chr" );
        fm.Cbind( *m_pVarsel_SnpPos, (char*)"Pos" );
        fm.Cbind( *m_pVarsel_Ra );
        fm.Cbind( *m_pVarsel_Rd );
        fm.SetRowNames(m_pVarsel_SnpName);
    }
    else
    {
        fm.Cbind( *m_pRefit_SnpChr, (char*)"Chr" );
        fm.Cbind( *m_pRefit_SnpPos, (char*)"Pos" );
        fm.Cbind( *m_pRefit_Ra );
        fm.Cbind( *m_pRefit_Rd );
        fm.SetRowNames(m_pRefit_SnpName);
    }

    fm.WriteAsCSVFile(szCsvFile, false);
    return(0);
}

int GLS_res::SaveSigFile( char* szCsvFile, bool bRefit )
{
    _log_info(_HI_, "SaveSigFile: CSV File=%s. bRefit=%d", szCsvFile, bRefit?1:0 );

    CFmMatrix fmSig(0,0);
    CFmVector VctSigA(0, 0.0);
    CFmVector VctSigD(0, 0.0);

    if (bRefit)
    {
        fmSig.Cbind( *m_pRefit_SnpChr );
        fmSig.Cbind( *m_pRefit_SnpPos );
        fmSig.Cbind( *m_pRefit_Ra  );
        fmSig.Cbind( *m_pRefit_Rd  );
	    fmSig.SetRowNames( m_pRefit_SnpName );
    }
    else
    {
        fmSig.Cbind( *m_pVarsel_SnpChr );
        fmSig.Cbind( *m_pVarsel_SnpPos );
        fmSig.Cbind( *m_pVarsel_Ra  );
        fmSig.Cbind( *m_pVarsel_Rd  );
	    fmSig.SetRowNames( m_pVarsel_SnpName );
    }

    fmSig.WriteAsCSVFile(szCsvFile, false);
    return(0);

}

int GLS_res::SortoutBestQ(CFmFileMatrix* pMatRet, CFmMatrix* pBestQ, int nSnp, int nCov)
{
	_log_debug(_HI_, "SortoutBestQ: AD matrix file(%s). SNP=%d nCov=%d", pMatRet->GetFileName(), nSnp, nCov);

	CFmVector fmQBest(LG, 0.0);
	CFmVector fmQPosMean(LG, 0.0);
	CFmVector fmQPosMin(LG, 0.0);
	CFmVector fmQPosMax(LG, 0.0);
	CFmVector fmQ(LG*4, 0.0);
	CFmVector fmTmp(0, 0.0);

	GetBestQInfo( pMatRet, 0, fmQBest, fmQPosMean, fmQPosMin, fmQPosMax );
	{
		fmQ.Append(fmQBest);
		fmTmp = fmQPosMean*fmQPosMean;
		fmQ.Put(fmTmp.Sum());
		fmQ.Append(fmQPosMean);
		fmTmp = fmQPosMin*fmQPosMin;
		fmQ.Put(fmTmp.Sum());
		fmQ.Append(fmQPosMin);
		fmTmp = fmQPosMax*fmQPosMax;
		fmQ.Put(fmTmp.Sum());
		fmQ.Append(fmQPosMax);
	}

	for (long int i=0; i<nCov-1; i++)
    {
		GetBestQInfo( pMatRet, i+1, fmQBest, fmQPosMean, fmQPosMin, fmQPosMax );

		fmQ.Append(fmQBest);
		fmTmp = fmQPosMean*fmQPosMean;
		fmQ.Put(fmTmp.Sum());
		fmQ.Append(fmQPosMean);
		fmTmp = fmQPosMin*fmQPosMin;
		fmQ.Put(fmTmp.Sum());
		fmQ.Append(fmQPosMin);
		fmTmp = fmQPosMax*fmQPosMax;
		fmQ.Put(fmTmp.Sum());
		fmQ.Append(fmQPosMax);
	}

	for (long int i=0; i< nSnp; i++)
    {
		fmQ.Resize(0);
		GetBestQInfo( pMatRet, nCov+i*2, fmQBest, fmQPosMean, fmQPosMin, fmQPosMax );

		fmQ.Append(fmQBest);
		fmTmp = fmQPosMean*fmQPosMean;
		fmQ.Put(fmTmp.Sum());
		fmQ.Append(fmQPosMean);
		fmTmp = fmQPosMin*fmQPosMin;
		fmQ.Put(fmTmp.Sum());
		fmQ.Append(fmQPosMin);
		fmTmp = fmQPosMax*fmQPosMax;
		fmQ.Put(fmTmp.Sum());
		fmQ.Append(fmQPosMax);

		GetBestQInfo( pMatRet, nCov+i*2+1, fmQBest, fmQPosMean, fmQPosMin, fmQPosMax );

		fmQ.Append(fmQBest);
		fmTmp = fmQPosMean*fmQPosMean;
		fmQ.Put(fmTmp.Sum());
		fmQ.Append(fmQPosMean);
		fmTmp = fmQPosMin*fmQPosMin;
		fmQ.Put(fmTmp.Sum());
		fmQ.Append(fmQPosMin);
		fmTmp = fmQPosMax*fmQPosMax;
		fmQ.Put(fmTmp.Sum());
		fmQ.Append(fmQPosMax);

		pBestQ->SetRow(i, fmQ);
    }

    _log_debug(_HI_, "SortoutBestQ: Successful to sort out.");
    return(0);
}

int GLS_res::GetBestQInfo(CFmFileMatrix* pMatRet, int idx, CFmVector& fmQBest, CFmVector& fmQPosMean, CFmVector& fmQPosMin, CFmVector& fmQPosMax)
{
    int nMcmc= m_pCfg->m_nMaxIter - m_pCfg->GetBurnInRound();

	fmQBest.Resize( LG );
	fmQPosMean.Resize( LG );
	fmQPosMin.Resize( LG );
	fmQPosMax.Resize( LG );

	CFmVector fmVct(nMcmc, 0.0);
	CFmVector fmTmp(nMcmc, 0.0);

    for(int k=0; k< LG ; k++)
    {
		fmQBest[k]    = R_NaN;
		fmQPosMean[k] = R_NaN;
		fmQPosMin[k]  = R_NaN;
		fmQPosMax[k]  = R_NaN;
        int ret = pMatRet->GetCacheCol(idx*LG + k, fmVct);
        if(ret!=0) return(ret);

        fmVct.Sort();

		for(int qx=0 ; qx < nMcmc/2; qx++)
			if(( fmVct[qx]>0 && fmVct[nMcmc-1-qx]>0 ) || ( fmVct[qx]<0 && fmVct[nMcmc-1-qx]<0 ) )
			{
				fmQBest[k]    = qx/(nMcmc*1.0 );

				fmTmp.Resize(0);
				for(int qxi = qx; qxi <= nMcmc-1-qx; qxi++)
					fmTmp.Put(fmVct[qxi]);

				fmQPosMean[k] = fmTmp.GetMedian();
				fmQPosMin[k]  = fmTmp.GetMin();
				fmQPosMax[k]  = fmTmp.GetMax();
				break;
			}

		if(isnan(fmQBest[k]))
		{
			fmQBest[k]    = 0.5;
			fmQPosMean[k] = fmVct[nMcmc/2];
			fmQPosMin[k]  = fmVct[nMcmc/2];
			fmQPosMax[k]  = fmVct[nMcmc/2];
		}
	}

	return(0);
}