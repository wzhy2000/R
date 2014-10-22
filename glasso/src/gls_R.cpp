/* gls_app.cpp  -	GLS application
 *	Copyright (C) 2011 THe Center for Statistical Genetics
 *  http://statgen.psu.edu
 */

#include "fm_rlogger.h"
#include "fm_pcf.h"
#include "fm_err.h"
#include "gls_R.h"
#include "gls_model.h"
#include "gls_cfg.h"
#include "gls_par.h"

int _glasso_simulate( CMDOPTIONS *pCmd, GLS_par *pPar)
{
    _log_prompt(_HI_, "Simulation will be performed.\n Total rounds:%d parameter:%s, pcf_file=%s",
              pCmd->nSimuRound, pCmd->szParFile, pCmd->szPcfFile);

    CFmPcf pcf(pCmd->szPcfFile);
    pcf.UpdatePcfFile( PCF_START );

    int status = 0;
    try
    {
        GLS sm;
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


SEXP _glasso_plink( CMDOPTIONS* pCmd, GLS_cfg* pCfg )
{
    _log_prompt(_HI_, "PLINK will be performed.\n Total rounds:%d parameter:%s,%s,%s, pcf_file=%s",
              pCmd->nSimuRound, pCmd->szTpedFile, pCmd->szTfamFile, pCmd->szPheFile, pCmd->szPcfFile);

	SEXP sRet = R_NilValue;
    CFmPcf pcf(pCmd->szPcfFile );
    pcf.UpdatePcfFile( PCF_START );

    int status = 0;
    try
    {
        GLS sm;
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
    return(sRet);

_Abort2:
    pcf.UpdatePcfFile( PCF_EXCEPT );
    _log_prompt( _HI_, "PLINK is exit abnormally with code(%d),", status);
    return(sRet);
}

SEXP _glasso_simple(  CMDOPTIONS* pCmd, GLS_cfg* pCfg )
{
    _log_prompt(_HI_, "SIMPLE will be performed. \nTotal rounds:%d parameter:%s,%s,%s, pcf_file=%s",
              pCmd->nSimuRound, pCmd->szTpedFile, pCmd->szTfamFile, pCmd->szPheFile, pCmd->szPcfFile);

	SEXP sRet = R_NilValue;
    CFmPcf pcf( pCmd->szPcfFile );
    pcf.UpdatePcfFile( PCF_START );

    int status = 0;
    try
    {
        GLS sm;
        status = sm.LoadSimple( pCmd );
        if (status!=0)
            goto _Abort3;

		status = sm.Varsel(pCfg);
		if (status!=0)
			goto _Abort3;

        if(pCmd->nRunmode != _RUN_VARSEL)
        {
            status = sm.Refit(pCfg);
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
    return(sRet);
}

int glasso_simulate( const char* szPhe_out,
					 const char* szSnp_out,
					 int nSimu_grp,
					 int nSimu_n,
					 int nSimu_p,
					 double fSimu_snp_rho,
					 double fSimu_rho,
					 double fSimu_sigma2,
					 double* pfSimu_mu,
					 int nSimu_covar_len,
					 double* pfSimu_covar_effect,
					 double* pfSimu_covar_range,
					 int nSimu_sig_p,
					 int nSimu_a_len,
					 double* pfSimu_a_effect,
			 		 int nSimu_d_len,
			 		 double* pfSimu_d_effect,
			 		 double* pfSimu_z_range,
			 		 int* pnSimu_z_count,
					 int bDebug)
{
	CMDOPTIONS cmd;
	memset(&cmd, 0, sizeof(CMDOPTIONS));

	strcpy( cmd.szPheoutFile, szPhe_out);
	strcpy( cmd.szSnpoutFile, szSnp_out);
	cmd.bDebug = bDebug;

    start_log( cmd.bDebug );

	GLS_par *pPar = new GLS_par(&cmd);

    pPar->simu_grps    = nSimu_grp;
    pPar->simu_n       = nSimu_n;
	pPar->simu_p       = nSimu_p;
    pPar->simu_snp_rho = fSimu_snp_rho;
	pPar->simu_rho     = fSimu_rho;
	pPar->simu_sigma2  = fSimu_sigma2;

	memcpy(pPar->simu_mu, pfSimu_mu, sizeof(double)*4);

	pPar->simu_covar_len = nSimu_covar_len;
	for(int i=0;i<nSimu_covar_len;i++)
	{
		pPar->simu_covar_effect[i][0] = pfSimu_covar_effect[nSimu_covar_len*0+i];
		pPar->simu_covar_effect[i][1] = pfSimu_covar_effect[nSimu_covar_len*1+i];
		pPar->simu_covar_effect[i][2] = pfSimu_covar_effect[nSimu_covar_len*2+i];
		pPar->simu_covar_effect[i][3] = pfSimu_covar_effect[nSimu_covar_len*3+i];
	};

	pPar->sig_p			 = nSimu_sig_p;

	pPar->simu_a_len   = nSimu_a_len;
	for(int i=0;i<nSimu_a_len;i++)
	{
		pPar->simu_a_pos[i] = (int)round(pfSimu_a_effect[i]);
		pPar->simu_a_effect[i][0] = pfSimu_a_effect[nSimu_a_len+i];
		pPar->simu_a_effect[i][1] = pfSimu_a_effect[nSimu_a_len*2+i];
		pPar->simu_a_effect[i][2] = pfSimu_a_effect[nSimu_a_len*3+i];
		pPar->simu_a_effect[i][3] = pfSimu_a_effect[nSimu_a_len*4+i];
	};

	pPar->simu_d_len   = nSimu_d_len;
	for(int i=0;i<nSimu_d_len;i++)
	{
		pPar->simu_d_pos[i] = (int)round(pfSimu_d_effect[i]);
		pPar->simu_d_effect[i][0] = pfSimu_d_effect[nSimu_d_len+i];
		pPar->simu_d_effect[i][1] = pfSimu_d_effect[nSimu_d_len*2+i];
		pPar->simu_d_effect[i][2] = pfSimu_d_effect[nSimu_d_len*3+i];
		pPar->simu_d_effect[i][3] = pfSimu_d_effect[nSimu_d_len*4+i];
	};

	memcpy(pPar->simu_z_count, pnSimu_z_count, sizeof(int)*2);

	pPar->simu_z_range[0] = (int)round(pfSimu_z_range[0]);
	pPar->simu_z_range[1] = (int)round(pfSimu_z_range[1]);

	pPar->simu_covar_range[0] = pfSimu_covar_range[0];
	pPar->simu_covar_range[1] = pfSimu_covar_range[1];

    int nRet = _glasso_simulate( &cmd, pPar );
	delete pPar;

    return(nRet);
}

SEXP glasso_simple( const char* pszPhefile,
  		   	const char*  pzSnpfile,
  		   	const char*  pzModel,
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

	strcpy( cmd.szSnpFile,  pzSnpfile);
	strcpy( cmd.szPheFile,  pszPhefile);
	strcpy( cmd.szModel,    pzModel);
	cmd.bDebug = bDebug;

    start_log( cmd.bDebug );

	GLS_cfg *pCfg = new GLS_cfg();
	pCfg->m_nMaxIter  	  = nMaxIter,
    pCfg->m_fBurnInRound  = fBurnInRound;
	pCfg->m_fQval_add	  = fQval_add;
	pCfg->m_fQval_dom	  = fQval_dom;
    pCfg->m_fRhoTuning    = fRhoTuning;

    if ( strlen(cmd.szSnpFile)==0 ||
         strlen(cmd.szPheFile)==0 )
    {
        fprintf(stderr, "\nThe SNP and phenotype file are required for a SIMPLE command.( options: -snp -phe ).\n\n");
        //return( ERR_OPTION_VALUE );
        return( R_NilValue );
    }

    if (strlen(cmd.szMatadFile)==0)
        CFmSys::GetTempFile(cmd.szMatadFile, "gls.fmat", MAX_PATH);

	SEXP sRet;
    sRet = _glasso_simple( &cmd, pCfg );
    delete pCfg;

    return(sRet);
}

SEXP glasso_plink( const char* pszPhefile,
  		   	const char*  pzTpedfile,
  		   	const char*  pzTfamfile,
  		   	const char*  pzModel,
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

	strcpy( cmd.szPheFile,  pszPhefile);
	strcpy( cmd.szTpedFile, pzTpedfile);
	strcpy( cmd.szTfamFile, pzTfamfile);
	strcpy( cmd.szModel,    pzModel);
	cmd.bDebug = bDebug;

    start_log( cmd.bDebug );

	GLS_cfg *pCfg = new GLS_cfg();
	pCfg->m_nMaxIter  	  = nMaxIter,
    pCfg->m_fBurnInRound  = fBurnInRound;
	pCfg->m_fQval_add     = fQval_add;
	pCfg->m_fQval_dom     = fQval_dom;
    pCfg->m_fRhoTuning    = fRhoTuning;

    if ( strlen(cmd.szSnpFile)==0 ||
         strlen(cmd.szPheFile)==0 )
    {
        fprintf(stderr, "\nThe SNP and phenotype file are required for a SIMPLE command.( options: -snp -phe ).\n\n");
        //return( ERR_OPTION_VALUE );
        return(R_NilValue);
    }

    if (strlen(cmd.szMatadFile)==0)
        CFmSys::GetTempFile(cmd.szMatadFile, "gls.fmat", MAX_PATH);

	SEXP sRet;
    sRet = _glasso_plink( &cmd, pCfg );
    delete pCfg;

    return(sRet);
}