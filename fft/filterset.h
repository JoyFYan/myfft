// filterset.h: interface for the filterset class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_FILTERSET_H__3D7A399F_1E10_47F8_BC4D_A01A9A91F6DF__INCLUDED_)
#define AFX_FILTERSET_H__3D7A399F_1E10_47F8_BC4D_A01A9A91F6DF__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "filter.h"

class filterset  
{
public:
	int symmetric;//对称性
	
	//底高通分解和合成滤波器
	filter *analysisLow, *analysisHigh, *synthesisLow,*synthesisHigh;
	
	filterset();//默认构造函数
	filterset (const filterset &filterset);
	filterset (int symmetric,double *anLow, int anLowSize, 
		int anLowFirst,double *synLow = NULL, 
		int synLowSize =0,int synLowFirst = 0);
	virtual ~filterset();
	
	void copy (const filterset& filterset);

	
protected:

};

#endif // !defined(AFX_FILTERSET_H__3D7A399F_1E10_47F8_BC4D_A01A9A91F6DF__INCLUDED_)
