#ifndef _FG_COVAR_H_
#define _FG_COVAR_H_

class CFmVectorStr;
class CFmVector;
class CFmMatrix;

class CFgCovariance
{
public:
    CFgCovariance(int nCovar);
    virtual ~CFgCovariance();

    virtual int GetValue(CFmMatrix* pValue, CFmVector* pTime);	
    virtual int IntialParam(double* pPar);
private:
    int m_nCovar;		
};

class CFgCovariance_AR1 : public CFgCovariance
{
public:
    CFgCovariance_AR1 (int nCovar );
    virtual ~CFgCovariance_AR1();

    virtual int GetValue(CFmMatrix* pValue, CFmVector* pTime);	
    virtual int IntialParam(double* pPar);
private:
	
};

CFgCovariance* CreateCovariance( int nCovar );

void destroy( CFgCovariance* p);

#endif // _FG_COVAR_H_
