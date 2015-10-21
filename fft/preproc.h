/************************************************************************
* 版权所有 (C)2011, 电子科技大学.图形图像与信号处理应用研究室
* 
* 文件名称：  preproc.h
* 文件标识：  
* 内容摘要：  图像预处理模块头文件
* 其它说明：  
* 当前版本：  
* 作    者：  
* 完成日期：  
* 
************************************************************************/
#ifndef C_PREPROC_H
#define C_PREPROC_H

#include "../main/stdafx.h"
#include "../main/IFTrack.h"
#include "../os/PublicFunction.h"
#include "math.h"
#include "../feature/wavelett.h"
#include "filter.h"
#include "../preproc/filterset.h"
#include <afxtempl.h>
#include <vector>
using namespace std;


#include "../segmentation/segmentation.h"
#include "../feature/feature.h"
#include "../detect/detect.h"
#include "../target/target.h"
#include "../track/track.h"
#include "../track/trackstrategy.h"


#define LINEBYTES(bits)  (((bits)+31)/32*4)
#define PI 3.1415926535

/**************************************************************************
*                        常量                                            *
**************************************************************************/


/**************************************************************************
*                          宏定义                                         *
**************************************************************************/
class FloatImage
{
public:
	float *  lp_data;
	float ** lp_AddRow;
	unsigned int Width;   // x indices range from 0 .. width-1
	unsigned int Height;  // y indices range from 0 .. height-1	
	unsigned int ImageSize;
public:
	FloatImage();
	~FloatImage();
	FloatImage(unsigned int w, unsigned int h);
	void Construct(unsigned int w, unsigned int h);
	void DeleteData();
};
 

 
class CImage
{
public:
	long       	  m_ImageSize;			 // 图像的大小 
    long          m_ImageWidth;			 // 图像的宽度
	long          m_ImageHeight;			 // 图像的高度
	unsigned char  **RowAddress;
	unsigned char  *m_lpDibArray;     // 数据的指针
protected:
private:

public:
	CImage();
	~CImage();
	CImage(unsigned int w, unsigned int h);
	void Construct1(unsigned int w, unsigned int h);
	void DeleteData1();
};
 


/**************************************************************************
*                            数据类型                                     *
**************************************************************************/




/**************************************************************************
*                             类声明                                      *
**************************************************************************/




/**************************************************************************
*                           模板                                         *
**************************************************************************/



/**************************************************************************
*                         全局变量声明                                    *
**************************************************************************/



/**************************************************************************
*                        全局函数原型                                     *
**************************************************************************/
/************************************************************************/
/*                                                                      */
/*                    一、红外图像非均匀校正模块                        */
/*                                                                      */
/************************************************************************/
void TwopointNUC(unsigned char * image,int width,int height,double * G , double * O );  //两点非均匀校正
void HighpassNUC(unsigned char* image , int width , int height, double * F,double M);  //高通非均匀校正
void ANNNUC(unsigned char * image, int width ,int height,float * G,float *O);          //神经网络非均匀校正
void Static_KalmanNUC(unsigned char * image , int width ,int height , double ** X);    //稳态kalman滤波的非均匀校正


/************************************************************************/
/*                                                                      */
/*                    二、红外图像噪声抑制与增强模块                    */
/*                                                                      */
/************************************************************************/

void Medianfilter_3x3(unsigned char* image,int width,int height);                //中值滤波
void Near_Mean(unsigned char* image,int width,int height, int Ksize=3);          //均值滤波
void Gaussian(unsigned char* image,int width,int height, int  Ksize = 3, 
			  float sigma = 1.0);                                                //高斯滤波
void DifGaussian(unsigned char* image,int width,int height, int  Ksize = 3 );    //各向异性高斯滤波
void HMT_denoise(unsigned char *image,int width,int height,int wavelet_style,
				 int layer);                                                         //基于马尔可夫小波域去噪
void Image_SUSAN_Filter(BYTE* InputImage , int ImgW , int ImgH, int bt, double dt) ;  //Susan滤波

void Double_Dimension_histEqual(unsigned char* image,int width,int height);      //直方图双向均衡增强
void LateralInhibition(unsigned char* image,int width,int height, int  Ksize=3); //侧抑制网络增强


void Image_Reserve(unsigned char *image,int width,int height);    //图像反色;

/************************************************************************/
/*                       边缘检测 										*/
/************************************************************************/	
void  Roberts(unsigned char *image,int width,int height);
void  Sobel(unsigned char *image,int width,int height);
BOOL  Canny(LPSTR lpBits,long lWidth, long lHeight, long lThresh, long hThresh);
void  Prewitt(unsigned char *image,int width,int height);
void  Laplacian(unsigned char *image,int width,int height);
// 旋转不变矩算子提取边缘
void  Edge_RotateInvariantOperator(BYTE *InputImage, int ImgW, int ImgH, int Masksize);

BOOL VGaussianKernel		(int ncoeffs, float *coeffs, float s);
BOOL VGaussianD1Kernel		(int ncoeffs, float *coeffs, float s);
BOOL VGaussianSplineKernel  (int ncoeffs, float *coeffs, float s);
void Make_Rfilter_Mask(int masksize, float lambda,  
					   float *lpDmask, float *lpmask= NULL, float maxresponse= 1.0f);
void Make_Second_Rfilter_Mask(int masksize, float lambda, float tau, 
							  float *lpDmask, float *lpmask= NULL, float maxresponse= 1.0f);
void Make_GED_filter_mask(int masksize, float lambda, float tau, 
						  float *lpDmask, float *lpmask= NULL, float maxresponse= 1.0f);
void  Image_SUSAN_Principle(BYTE* InputImage , int ImgW , int ImgH ,int bt);


float Image_Edge_Canny_EstimateNoise(FloatImage &mag, float amax, float q);

void    Image_Edge_Canny_NonMaximalSuppress 
(CImage     &dest, 
 FloatImage &xgrad, 
 FloatImage &ygrad,
 FloatImage &mag,
 float       amax);
void    Image_Edge_Canny(CImage		&InImage,
						 float		sigma,
						 float		&noise,
						 float		lambda= 8,
						 float		tau= 0,
						 int		MaskType= 0,
						 FloatImage *ori= NULL);

/************************************************************************/
/*                       图像锐化										*/
/************************************************************************/
//梯度锐化
void  Gradient_Sharper(unsigned char *image,int width,int height,int thresh);
//Roberts边缘增强
void  Roberts_Sharper(BYTE *InputImage,int ImgW,int ImgH);
//laplassian边缘增强
void  Laplacian_sharper(BYTE *InputImage,int ImgW,int ImgH);


//图像金字塔显示
void Image_Pyramid_Show (BYTE *InputImage, int ImgW, int ImgH);


/***********************************************************************/
/*            图像细化算法  四种不同的细化算法                          */
/************************************************************************/
void beforethin(unsigned char *ip, unsigned char *jp,  
				unsigned long lx, unsigned long ly); 
void ThinnerHilditch(void *image, unsigned long lx, unsigned long ly); 
void ThinnerPavlidis(void *image, unsigned long lx, unsigned long ly); 
void ThinnerRosenfeld(void *image, unsigned long lx, unsigned long ly); 
//注意该函数lWidth应该是Height； 

/************************************************************************/
/*                    图像手动二值化                                 */
/************************************************************************/
void Image_binary(unsigned char*image,int width,int height,int thresh);


#endif  /* C_PREPROC_H */

