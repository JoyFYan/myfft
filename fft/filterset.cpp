// filterset.cpp: implementation of the filterset class.
//
//////////////////////////////////////////////////////////////////////

#include "../main/stdafx.h"
#include "../main/IFTrack.h"
#include "filterset.h"
#include <assert.h>


#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

double Daub4Coeffs [] = { 0.4829629131445341,  0.8365163037378077,
		          0.2241438680420134, -0.1294095225512603 };

filterset::filterset()
{
	symmetric = FALSE;
	analysisLow = analysisHigh =
		synthesisLow = synthesisHigh = NULL;
}

filterset::~filterset()
{
  delete analysisLow;
  delete analysisHigh;
  delete synthesisLow;
  delete synthesisHigh;
}

//////////////////////////////////////////////////////////////////////
//功能：通过一个低通分解滤波器和一个低通合成滤波器构造一个具有高低分解和高低通合成的
//输入参数：symmetric系数的对称性，anLow低通分解, anLowSize低通大小，anLowFirst低通初始点，
//synLow合成低通,synLowSize合成大小，synLowFirst初始点
//输出参数：无
//没有返回值
//////////////////////////////////////////////////////////////////////
filterset::filterset(int symmetric, double *anLow, int anLowSize, 
					 int anLowFirst, double *synLow, int synLowSize, int synLowFirst)
{
 int i, sign;

  analysisLow = new filter (anLowSize, anLowFirst, anLow);
  
  // If no synthesis coeffs are given, assume wavelet is orthogonal
  if (synLow == NULL)  {
    synthesisLow = new filter (*analysisLow);

    // For orthogonal wavelets, compute the high pass filter using
    // the relation g_n = (-1)^n h_{1-n}^*
    // (or equivalently g_{1-n} = (-1)^{1-n} h_n^*)

    analysisHigh = new filter (analysisLow->size, 2 - analysisLow->size -
				analysisLow->firstIndex,anLow);
      
    // Compute (-1)^(1-n) for first n
    if (analysisLow->firstIndex % 2)
      sign = 1;
    else sign = -1;
    
    for (i = 0; i < analysisLow->size; i++)  {
      analysisHigh->coeff[1 - i - analysisLow->firstIndex - 
		           analysisHigh->firstIndex] = 
	sign * analysisLow->coeff[i];
      assert (1 - i - analysisLow->firstIndex - 
		           analysisHigh->firstIndex >= 0);
      assert (1 - i - analysisLow->firstIndex - 
		           analysisHigh->firstIndex < analysisHigh->size);
      sign *= -1;
    }

    // Copy the high pass analysis filter to the synthesis filter
    synthesisHigh = new filter (*analysisHigh);    
    
  } else {
    // If separate synthesis coeffs given, assume biorthogonal
    
    synthesisLow = new filter (synLowSize, synLowFirst, synLow);

    // For orthogonal wavelets, compute the high frequency filter using
    // the relation g_n = (-1)^n complement (h~_{1-n}) and
    //              g~_n = (-1)^n complement (h_{1-n})
    // (or equivalently g_{1-n} = (-1)^{1-n} complement (h~_n))
    
    analysisHigh = new filter (synthesisLow->size, 2 - synthesisLow->size -
			 synthesisLow->firstIndex,synLow);
    
    // Compute (-1)^(1-n) for first n
    if (synthesisLow->firstIndex % 2)
      sign = 1;
    else sign = -1;
    
    for (i = 0; i < synthesisLow->size; i++)  {
      analysisHigh->coeff[1 - i - synthesisLow->firstIndex -
			  analysisHigh->firstIndex] = 
	sign * synthesisLow->coeff[i];
      assert (1 - i - synthesisLow->firstIndex - 
		           analysisHigh->firstIndex >= 0);
      assert (1 - i - synthesisLow->firstIndex - 
		           analysisHigh->firstIndex < analysisHigh->size);
      sign *= -1;
    }

    synthesisHigh = new filter 
                           (analysisLow->size, 2 - analysisLow->size -
			    analysisLow->firstIndex,anLow); 
    
    // Compute (-1)^(1-n) for first n
    if (analysisLow->firstIndex % 2)
      sign = 1;
    else sign = -1;

    for (i = 0; i < analysisLow->size; i++)  {
      synthesisHigh->coeff[1 - i - analysisLow->firstIndex -
			   synthesisHigh->firstIndex] = 
	sign * analysisLow->coeff[i];
      assert (1 - i - analysisLow->firstIndex - 
		           synthesisHigh->firstIndex >= 0);
      assert (1 - i - analysisLow->firstIndex - 
		           synthesisHigh->firstIndex < synthesisHigh->size);
      sign *= -1;
    }
  }
}

//////////////////////////////////////////////////////////////////////
//功能：一套滤波设定实例对另外一个滤波设定实例的复制
//输入参数：filterset 一个滤波设定实例
//输出参数：无
//没有返回值
//////////////////////////////////////////////////////////////////////
void filterset::copy(const filterset &filterset)
{
  symmetric = filterset.symmetric;
  analysisLow = new filter (*(filterset.analysisLow));
  analysisHigh = new filter (*(filterset.analysisHigh));
  synthesisLow = new filter (*(filterset.synthesisLow));
  synthesisHigh = new filter (*(filterset.synthesisHigh));
}

//////////////////////////////////////////////////////////////////////
//功能：一套滤波设定实例对另外一个滤波设定实例的复制
//输入参数：filterset 一个滤波设定实例
//输出参数：无
//没有返回值
//////////////////////////////////////////////////////////////////////
filterset::filterset(const filterset &filterset)
{
	copy (filterset);//直接调用copy函数来完成实例的初始化
}
