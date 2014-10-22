// bls_app.h: interface for the blasso command.
//
//////////////////////////////////////////////////////////////////////

#if !defined(BLS_APP_H__INCLUDED_)
#define BLS_APP_H__INCLUDED_

#include <stdbool.h>
#include "fm_linux.h"

#ifdef __cplusplus
extern "C" {
#endif

#define _RUN_BOTH       0   //Step 1+2: both
#define _RUN_VARSEL     1   //Step 1: variable selection
#define _RUN_REFIT      2   //Step 2: Refit

#define MAX_RESIN_CNT      128

typedef struct cmd_param{
    bool bHelp;
    bool bVersion;
    bool bDebug;
    bool bZNormalized;
    int  nRunmode;
    int  nSimuRound;
    char szAppName[MAX_PATH];
    char szParFile[MAX_PATH];
    char szTpedFile[MAX_PATH];
    char szTfamFile[MAX_PATH];
    char szVsretFile[MAX_PATH];
    char szLogFile[MAX_PATH];
    char szPheFile[MAX_PATH];
    char szSnpFile[MAX_PATH];
    char szModel[MAX_PATH];
    char szMifFile[MAX_PATH];
    char szPcfFile[MAX_PATH];
    char szCfgFile[MAX_PATH];
    char szPreSigFile[MAX_PATH];
    char szResFile[MAX_PATH];
    char szRdataFile[MAX_PATH];
    char szPheoutFile[MAX_PATH];
    char szSnpoutFile[MAX_PATH];
    char szSigOutFile[MAX_PATH];
    char szFigOutFile[MAX_PATH];
    char szRtdump[MAX_PATH];
    char szTmpFmat[MAX_PATH];
} CMDOPTIONS;

int blasso_simulate( const char* szPhe_out,
			 const char* szSnp_out,
			 int nSimu_grp,
			 int nSimu_n,
			 int nSimu_p,
			 double fsimu_snp_rho,
			 double fSimu_rho,
			 double fSimu_sigma2,
			 double fSimu_mu,
			 int nSimu_cov_count,
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
			 int bDebug );

SEXP blasso_plink(  const char* pszPhefile,
  		   	const char*  pzTpedFile,
  		   	const char*  szTfamFile,
  		   	const char*  pzModel,
  		   	bool bRefit,
  		   	int nMaxIter,
		   	double fBurnInRound,
		   	double fRhoTuning,
	        	double fQval_add,
	        	double fQval_dom,
	        	bool   bDebug);

SEXP blasso_simple( const char* pszPhefile,
  		   	const char*  pzSnpfile,
  		   	const char*  pzModel,
  		   	bool bRefit,
  		   	int nMaxIter,
		   	double fBurnInRound,
		   	double fRhoTuning,
	        	double fQval_add,
	        	double fQval_dom,
		        bool   bDebug);



#ifdef __cplusplus
}
#endif

#endif


