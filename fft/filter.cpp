// filter.cpp: implementation of the filter class.
//
//////////////////////////////////////////////////////////////////////

#include "../main/stdafx.h"
#include "../main/IFTrack.h"
#include "filter.h"


#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////
filter::filter()
{
	coeff = NULL; size = 0; firstIndex = 0;
}

filter::filter(int size, int firstIndex, double *coeff)
{
	init(size,firstIndex,coeff);
}

filter::filter(const filter &filter)
{
	coeff=NULL; 
	copy(filter);
}

filter::~filter()
{
	if (coeff != NULL)
		delete [] coeff;
}

//////////////////////////////////////////////////////////////////////
//功能：滤波器的初始化
//输入参数：filterSize是小波的长度，data是小波的系数，filterFirst滤波开始数
//输出参数：size 大小，firstIndex 初始数，coeff 是系数数列
//没有返回值
//////////////////////////////////////////////////////////////////////
void filter::init(int filterSize, int filterFirst, double *data)
{
	size = filterSize;
	firstIndex = filterFirst;
	center = -firstIndex;
	
	coeff = new double [size];
	if (data != NULL) 
	{
		for (int i = 0; i < size; i++)
			coeff[i] = data[i];
	}
	else//判断是否传入小波系数
	{
		for (int i = 0; i < size; i++)
			coeff[i] = 0;
	}
}

//////////////////////////////////////////////////////////////////////
//功能：滤波器的初始化
//输入参数：filter一个小波实例
//输出参数：另一个小波类实例
//没有返回值
//////////////////////////////////////////////////////////////////////
void filter::copy(const filter &filter)
{
	if (coeff != NULL)
		delete [] coeff;//对全局小波系数指针的判定
	//初始化
	init (filter.size, filter.firstIndex, filter.coeff);
}



