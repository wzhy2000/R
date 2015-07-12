#ifndef _FG_CURVE_H_
#define _FG_CURVE_H_

class CFmVectorStr;
class CFmVector;
class CFmMatrix;

class CFgCurve
{
public:
    CFgCurve(int nCurve, CFmVector* pTime);
    virtual ~CFgCurve();

    virtual int GetValue(CFmVector* pValue, CFmVector* pTime);	
    virtual int IntialParam(double* pPar);
private:
    int m_nCurve;	
};

class CFgCurve_Log : public CFgCurve
{
public:
    CFgCurve_Log(int nCurve, CFmVector* pTime);
    virtual ~CFgCurve_Log();

    virtual int GetValue(CFmVector* pValue, CFmVector* pTime);	
    virtual int IntialParam(double* pPar);
private:
	
};

CFgCurve* CreateCurve(int nCurve);

void destroy( CFgCurve* p);

#endif // _FG_CURVE_H_
