/************************************************************************
* ��Ȩ���� (C)2011, ���ӿƼ���ѧ.ͼ��ͼ�����źŴ���Ӧ���о���
* 
* �ļ����ƣ�  preproc.h
* �ļ���ʶ��  
* ����ժҪ��  ͼ��Ԥ����ģ��ͷ�ļ�
* ����˵����  
* ��ǰ�汾��  
* ��    �ߣ�  
* ������ڣ�  
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
*                        ����                                            *
**************************************************************************/


/**************************************************************************
*                          �궨��                                         *
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
	long       	  m_ImageSize;			 // ͼ��Ĵ�С 
    long          m_ImageWidth;			 // ͼ��Ŀ��
	long          m_ImageHeight;			 // ͼ��ĸ߶�
	unsigned char  **RowAddress;
	unsigned char  *m_lpDibArray;     // ���ݵ�ָ��
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
*                            ��������                                     *
**************************************************************************/




/**************************************************************************
*                             ������                                      *
**************************************************************************/




/**************************************************************************
*                           ģ��                                         *
**************************************************************************/



/**************************************************************************
*                         ȫ�ֱ�������                                    *
**************************************************************************/



/**************************************************************************
*                        ȫ�ֺ���ԭ��                                     *
**************************************************************************/
/************************************************************************/
/*                                                                      */
/*                    һ������ͼ��Ǿ���У��ģ��                        */
/*                                                                      */
/************************************************************************/
void TwopointNUC(unsigned char * image,int width,int height,double * G , double * O );  //����Ǿ���У��
void HighpassNUC(unsigned char* image , int width , int height, double * F,double M);  //��ͨ�Ǿ���У��
void ANNNUC(unsigned char * image, int width ,int height,float * G,float *O);          //������Ǿ���У��
void Static_KalmanNUC(unsigned char * image , int width ,int height , double ** X);    //��̬kalman�˲��ķǾ���У��


/************************************************************************/
/*                                                                      */
/*                    ��������ͼ��������������ǿģ��                    */
/*                                                                      */
/************************************************************************/

void Medianfilter_3x3(unsigned char* image,int width,int height);                //��ֵ�˲�
void Near_Mean(unsigned char* image,int width,int height, int Ksize=3);          //��ֵ�˲�
void Gaussian(unsigned char* image,int width,int height, int  Ksize = 3, 
			  float sigma = 1.0);                                                //��˹�˲�
void DifGaussian(unsigned char* image,int width,int height, int  Ksize = 3 );    //�������Ը�˹�˲�
void HMT_denoise(unsigned char *image,int width,int height,int wavelet_style,
				 int layer);                                                         //��������ɷ�С����ȥ��
void Image_SUSAN_Filter(BYTE* InputImage , int ImgW , int ImgH, int bt, double dt) ;  //Susan�˲�

void Double_Dimension_histEqual(unsigned char* image,int width,int height);      //ֱ��ͼ˫�������ǿ
void LateralInhibition(unsigned char* image,int width,int height, int  Ksize=3); //������������ǿ


void Image_Reserve(unsigned char *image,int width,int height);    //ͼ��ɫ;

/************************************************************************/
/*                       ��Ե��� 										*/
/************************************************************************/	
void  Roberts(unsigned char *image,int width,int height);
void  Sobel(unsigned char *image,int width,int height);
BOOL  Canny(LPSTR lpBits,long lWidth, long lHeight, long lThresh, long hThresh);
void  Prewitt(unsigned char *image,int width,int height);
void  Laplacian(unsigned char *image,int width,int height);
// ��ת�����������ȡ��Ե
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
/*                       ͼ����										*/
/************************************************************************/
//�ݶ���
void  Gradient_Sharper(unsigned char *image,int width,int height,int thresh);
//Roberts��Ե��ǿ
void  Roberts_Sharper(BYTE *InputImage,int ImgW,int ImgH);
//laplassian��Ե��ǿ
void  Laplacian_sharper(BYTE *InputImage,int ImgW,int ImgH);


//ͼ���������ʾ
void Image_Pyramid_Show (BYTE *InputImage, int ImgW, int ImgH);


/***********************************************************************/
/*            ͼ��ϸ���㷨  ���ֲ�ͬ��ϸ���㷨                          */
/************************************************************************/
void beforethin(unsigned char *ip, unsigned char *jp,  
				unsigned long lx, unsigned long ly); 
void ThinnerHilditch(void *image, unsigned long lx, unsigned long ly); 
void ThinnerPavlidis(void *image, unsigned long lx, unsigned long ly); 
void ThinnerRosenfeld(void *image, unsigned long lx, unsigned long ly); 
//ע��ú���lWidthӦ����Height�� 

/************************************************************************/
/*                    ͼ���ֶ���ֵ��                                 */
/************************************************************************/
void Image_binary(unsigned char*image,int width,int height,int thresh);


#endif  /* C_PREPROC_H */

