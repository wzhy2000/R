/* bls_app.cpp  -	BLS application
 *	Copyright (C) 2011 THe Center for Statistical Genetics
 *  http://statgen.psu.edu
 */

#include <R.h>
#include <Rinternals.h>
#include <Rembedded.h>
#include <Rdefines.h>

#include "fm_rlogger.h"
#include "fm_pcf.h"
#include "fm_err.h"
#include "bls_R.h"
#include "bls_model.h"

int _blasso_simulate( CMDOPTIONS *pCmd, BLS_par *pPar)
{
    _log_prompt(_HI_, "Simulation will be performed.\n Total rounds:%d parameter:%s, pcf_file=%s",
              pCmd->nSimuRound, pCmd->szParFile,  pCmd->szPcfFile);

    CFmPcf pcf(pCmd->szPcfFile);
    pcf.UpdatePcfFile( PCF_START );

    int status = 0;
    try
    {
        BLS sm;
        status = sm.LoadSimulate( pCmd, pPar );
        if (status!=0)
            goto _Abort1;

    }
    catch( char const* str )
    {
        _log_error( _HI_, "Simulation has an exception(%s).", str);
        status = ERR_EXCEPTION;
        goto _Abort1;
    }

    pcf.UpdatePcfFile( PCF_END );
    _log_prompt(_HI_, "Simulation is done successfully.");
    return(0);

_Abort1:
    pcf.UpdatePcfFile( PCF_EXCEPT, 0, 0, 0, 0, status );
    _log_prompt( _HI_, "Simulation is exit abnormally with code(%d),", status);
    return(status);
}

int blasso_simulate( const char* szPhe_out,
					 const char* szSnp_out,
					 int nSimu_grp,
					 int nSimu_n,
					 int nSimu_p,
					 double fsimu_snp_rho,
					 double fSimu_rho,
					 double fSimu_sigma2,
					 double fSimu_mu,
					 int nSimu_cov_len,
					 double* fSimu_cov_coeff,
					 int nSimu_sig_p,
					 int nSimu_a_len,
					 int* pnSimu_a_pos,
					 double* pfSimu_a_effect,
					 int nSimu_d_len,
					 int* pnSimu_d_pos,
					 double* pfSimu_d_effect,
					 double* pfSimu_covar_range,
					 double* pfSimu_t_range,
					 int bDebug )
{
	CMDOPTIONS cmd;
	memset(&cmd, 0, sizeof(CMDOPTIONS));

	strcpy( cmd.szPheoutFile, szPhe_out);
	strcpy( cmd.szSnpoutFile, szSnp_out);
	cmd.bDebug = bDebug;

    start_log( cmd.bDebug );

	BLS_par *pPar = new BLS_par(&cmd);

    pPar->simu_grps    = nSimu_grp;
    pPar->simu_n       = nSimu_n;
	pPar->simu_p       = nSimu_p;
    pPar->simu_snp_rho = fsimu_snp_rho;
	pPar->simu_rho     = fSimu_rho;
	pPar->simu_sigma2  = fSimu_sigma2;
	pPar->simu_mu      = fSimu_mu;

	pPar->simu_covar_len   = nSimu_cov_len;
	for(int i=0;i<nSimu_cov_len;i++)
		pPar->simu_covar_coefs[i] = fSimu_cov_coeff[i];

	pPar->sig_p		   = nSimu_sig_p;

	pPar->simu_a_len   = nSimu_d_len;
	for(int i=0;i<nSimu_a_len;i++)
	{
		pPar->simu_a_pos[i] = pnSimu_a_pos[i];
		pPar->simu_a_effect[i] = pfSimu_a_effect[i];
	};

	pPar->simu_d_len   = nSimu_d_len;
	for(int i=0;i<nSimu_d_len;i++)
	{
		pPar->simu_d_pos[i] = pnSimu_d_pos[i];
		pPar->simu_d_effect[i] = pfSimu_d_effect[i];
	};

	pPar->simu_t_range[0] = round(pfSimu_t_range[0]);
	pPar->simu_t_range[1] = round(pfSimu_t_range[1]);

	pPar->simu_covar_range[0] = round(pfSimu_covar_range[0]);
	pPar->simu_covar_range[1] = round(pfSimu_covar_range[1]);


    int nRet = _blasso_simulate( &cmd, pPar );
	delete pPar;

    return(nRet);

}

SEXP _blasso_plink(  CMDOPTIONS* pCmd, BLS_cfg *pCfg  )
{
    _log_prompt(_HI_, "PLINK will be performed.\n");

	SEXP sRet = R_NilValue;
    CFmPcf pcf(pCmd->szPcfFile );
    pcf.UpdatePcfFile( PCF_START );

    int status = 0;
    try
    {
        BLS sm;
        status = sm.LoadPlink( pCmd );
        if (status!=0)
            goto _Abort2;

		status = sm.Varsel(pCfg);
		if (status!=0)
			goto _Abort2;

        if(pCmd->nRunmode != _RUN_VARSEL)
        {
            status = sm.Refit(pCfg);
            if ( status!=0 && status!=ERR_NO_ITEMS )
                goto _Abort2;
        }

        sRet = sm.GetRObj();
    }
    catch( char const* str )
    {
        _log_error( _HI_, "PLINK has an exception(%s).", str);
        status = ERR_EXCEPTION;
        goto _Abort2;
    }

    pcf.UpdatePcfFile( PCF_END );
    _log_prompt(_HI_, "PLINK is done successfully.");
    return( sRet );

_Abort2:
    pcf.UpdatePcfFile( PCF_EXCEPT );
    _log_prompt( _HI_, "PLINK is exit abnormally with code(%d),", status);
    return(R_NilValue);
}

SEXP blasso_plink(  const char* pszPheFile,
  		   	const char*  pszTpedFile,
  		   	const char*  pszTfamFile,
  		   	const char*  pszModel,
  		   	bool bRefit,
  		   	int nMaxIter,
		   	double fBurnInRound,
		   	double fRhoTuning,
	        double fQval_add,
	        double fQval_dom,
	        bool   bDebug)
{
	CMDOPTIONS cmd;
	memset(&cmd, 0, sizeof(CMDOPTIONS));

	strcpy( cmd.szPheFile,  pszPheFile);
	strcpy( cmd.szTpedFile, pszTpedFile);
	strcpy( cmd.szTfamFile, pszTfamFile);
	strcpy( cmd.szModel,    pszModel);
	cmd.bDebug = bDebug;

    start_log( cmd.bDebug );

	BLS_cfg *pCfg = new BLS_cfg();
	pCfg->m_nMaxIter  	  = nMaxIter,
    pCfg->m_fRhoTuning    = fRhoTuning;
    pCfg->m_fBurnInRound  = fBurnInRound;
	pCfg->m_fQval_add	  = fQval_add;
	pCfg->m_fQval_dom	  = fQval_dom;

    if ( strlen(cmd.szTpedFile)==0 ||
         strlen(cmd.szTfamFile)==0 ||
         strlen(cmd.szPheFile)==0 )
    {
        fprintf(stderr, "\nThe TPED, TFAM and phenotype file are required for a PLINK command.( options: -tped -tfam -phe ).\n\n");
        fprintf(stderr, "\nTPED=%s\n", cmd.szTpedFile);
        fprintf(stderr, "\nTFAM=%s\n", cmd.szTfamFile);
        fprintf(stderr, "\nPHE=%s\n",  cmd.szPheFile);
        return( R_NilValue );
    }

    if ( cmd.nRunmode==_RUN_REFIT && strlen(cmd.szVsretFile)==0)
    {
        fprintf(stderr, "\nThe result file of variable selection procedure is required.( options: -vsret ).\n\n");
        return( R_NilValue );
    }

    SEXP sRet;
    sRet = _blasso_plink(&cmd, pCfg);
    delete pCfg;

    return(sRet);
}

SEXP _blasso_simple(  CMDOPTIONS* pCmd, BLS_cfg *pCfg )
{
    _log_prompt(_HI_, "SIMPLE will be performed. \nTotal rounds:%d parameter:%s,%s,%s, pcf_file=%s",
              pCmd->nSimuRound, pCmd->szTpedFile, pCmd->szTfamFile, pCmd->szPheFile, pCmd->szPcfFile);

    SEXP sRet = R_NilValue;
	CFmPcf pcf( pCmd->szPcfFile );
    pcf.UpdatePcfFile( PCF_START );

    int status = 0;
    try
    {
        BLS sm;
        status = sm.LoadSimple( pCmd );
        if (status!=0)
            goto _Abort3;

		status = sm.Varsel(pCfg);
		if (status!=0)
			goto _Abort3;

        if(pCmd->nRunmode != _RUN_VARSEL)
        {
            status = sm.Refit( pCfg );
            if ( status!=0 && status!=ERR_NO_ITEMS )
                goto _Abort3;
        }

        sRet = sm.GetRObj();

    }
    catch( char const* str )
    {
        _log_error( _HI_, "SIMPLE has an exception(%s).", str);
        status = ERR_EXCEPTION;
        goto _Abort3;
    }

    pcf.UpdatePcfFile( PCF_END );
    _log_prompt(_HI_, "SIMPLE is done successfully.");
    return(sRet);

_Abort3:
    pcf.UpdatePcfFile( PCF_EXCEPT, 0, 0, 0, 0, status );
    _log_prompt( _HI_, "SIMPLE is exit abnormally with code(%d),", status);
    return(R_NilValue);
}

SEXP blasso_simple( const char* pszPheFile,
  		   	const char*  pszSnpFile,
  		   	const char*  pszModel,
  		   	bool bRefit,
  		   	int nMaxIter,
		   	double fBurnInRound,
		   	double fRhoTuning,
	        double fQval_add,
	        double fQval_dom,
	        bool   bDebug)
{
	CMDOPTIONS cmd;
	memset(&cmd, 0, sizeof(CMDOPTIONS));

	strcpy( cmd.szSnpFile,  pszSnpFile);
	strcpy( cmd.szPheFile,  pszPheFile);
	strcpy( cmd.szModel,    pszModel);
	cmd.bDebug = bDebug;

    start_log( cmd.bDebug );

	BLS_cfg *pCfg = new BLS_cfg();
	pCfg->m_nMaxIter  	  = nMaxIter,
    pCfg->m_fRhoTuning    = fRhoTuning;
    pCfg->m_fBurnInRound  = fBurnInRound;
	pCfg->m_fQval_add	  = fQval_add;
	pCfg->m_fQval_dom	  = fQval_dom;

    if ( strlen(cmd.szSnpFile)==0 || strlen(cmd.szPheFile)==0 )
    {
        fprintf(stderr, "\nThe SNP and phenotype file are required for a SIMPLE command.( options: -snp -phe ).\n\n");
        return( R_NilValue );
    }

    if ( cmd.nRunmode==_RUN_REFIT && strlen(cmd.szVsretFile)==0)
    {
        fprintf(stderr, "\nThe result file of variable selection procedure is required.( options: -vsret ).\n\n");
        return( R_NilValue );
    }

    SEXP sRet;
    sRet = _blasso_simple(&cmd, pCfg);
    delete pCfg;

    return(sRet);
}
