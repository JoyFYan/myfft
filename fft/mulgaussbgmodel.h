// MulGaussBGModel.h: interface for the CMulGaussBGModel class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_MULGAUSSBGMODEL_H__997E0B4D_96C6_4CD8_8505_C0CC2C40995D__INCLUDED_)
#define AFX_MULGAUSSBGMODEL_H__997E0B4D_96C6_4CD8_8505_C0CC2C40995D__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000



#include<vector>
using std::vector;
/************************************************************************/
/*                                                                      */
/*                    �ġ�����Ŀ����ģ��                              */
/*                                                                      */
/************************************************************************/

/************************************************************************/
/*                    4.2 ��չĿ������غ���(��)                      */
/************************************************************************/ 

/************************************************************************/
/*                    A.���ڶ��˹������ģ���˶�Ŀ����                */
/************************************************************************/ 

/************************************************************************/
/*                          ��Ҫ��ͽṹ�������                        */
/************************************************************************/
#define BGFG_MOG_MAX_NGAUSSIANS 500

/* default parameters of gaussian background detection algorithm */
#define BGFG_MOG_BACKGROUND_THRESHOLD     0.7     /* threshold sum of weights for background test */
#define BGFG_MOG_STD_THRESHOLD            2.5     /* lambda=2.5 is 99% */
#define BGFG_MOG_WINDOW_SIZE              200     /* Learning rate; alpha = 1/CV_GBG_WINDOW_SIZE */
#define BGFG_MOG_NGAUSSIANS               5       /* = K = number of Gaussians in mixture */
#define BGFG_MOG_WEIGHT_INIT              0.05
#define BGFG_MOG_SIGMA_INIT               30

#define BGFG_MOG_NCOLORS                  3

//��˹����ģ�Ͳ����ṹ��
typedef struct GaussBGStatModelParams
{    
    int     win_size;                        /* = 1/alpha */
    int     n_gauss;                         //��˹�ֲ��ĸ���
    double  bg_threshold, std_threshold;
    double  weight_init, variance_init;
}GaussBGStatModelParams;

//��˹����ģ���еĸ�˹�ֲ������ṹ��
typedef struct GaussBGValues
{
    int         match_sum;                   //ͳ�Ƹ�˹�ֲ���ƥ��Ĵ���
    double      weight;                      //��˹�ֲ���Ȩ��
    double      variance[BGFG_MOG_NCOLORS];  //�����˹�ֲ��ķ���
    double      mean[BGFG_MOG_NCOLORS];      //�����˹�ֲ��ľ�ֵ
}
GaussBGValues;

//��˹����ģ�����ص�
typedef struct GaussBGPoint
{
    GaussBGValues* g_values;                  //��˹�ֲ�������ָ�루����Ϊ�����˹�ֲ���
}
GaussBGPoint;

//��˹�����ṹ��Ķ���
typedef struct GaussBGModel
{
	unsigned char*  background;                //���Ƴ��ı���ֵ 
	unsigned char*  foreground;                //��ֵǰ��ͼ(1ͨ��)                        
    GaussBGStatModelParams   params;           //��˹����ģ�Ͳ���
    GaussBGPoint*            g_point;          //��˹��������
    int                      countFrames;      //�Ѿ������֡������
	int width;                                 //ͼ��Ŀ��
	int height;                                //ͼ���֡��
}
GaussBGModel;
/************************************************************************/
/*                          ��Ҫ����������                              */
/************************************************************************/ 
GaussBGModel* CreateGaussianBGModel( unsigned char*Image,int Width,int Height,         //��ʼ֡������˹����ģ��     
						int nChannels, GaussBGStatModelParams* parameters );

void  UpdateGaussianBGModel(  unsigned char*Image,int Width,int Height,                //ͨ��������ͼ��ʱ���¸�˹����ģ��
							int nChannels, GaussBGModel*  bg_model );

void  ReleaseGaussianBGModel( GaussBGModel** _bg_model );                              //�������ͷŸ�˹ģ�������Դ
            
int   MatchTest( double* src_pixel, int nChannels, int* match,                         //ƥ����
	   const GaussBGPoint* g_point,const GaussBGStatModelParams *bg_model_params );

void  UpdateFullWindow( double* src_pixel, int nChannels, int* match,                  //ȫ������¶����и�˹�ֲ���ƥ�䣬��˹�ֲ����������ĸ���
		GaussBGPoint* g_point,const GaussBGStatModelParams *bg_model_params );

void  UpdatePartialWindow( double* src_pixel, int nChannels, int* match,               //���ִ�����¶����и�˹�ֲ���ƥ�䣬��˹�ֲ����������ĸ���
		GaussBGPoint* g_point, const GaussBGStatModelParams *bg_model_params );

void  UpdateFullNoMatch( unsigned char* gm_image, int nChannels,int p, int* match,     //ȫ�������û�и�˹�ֲ���ƥ�䣬��˹�ֲ����������ĸ���
		GaussBGPoint* g_point,const GaussBGStatModelParams *bg_model_params);

void  UpdatePartialNoMatch(double *pixel,int nChannels,int* /*match*/,                 //���ִ������û�и�˹�ֲ���ƥ�䣬��˹�ֲ����������ĸ���
			GaussBGPoint* g_point,const GaussBGStatModelParams *bg_model_params);

void  GetSortKey( const int nChannels, double* sort_key, const GaussBGPoint* g_point,  //�����˹�ֲ��ľ�ֵ�뷽��ı�ֵ
			const GaussBGStatModelParams *bg_model_params );

void  InsertionSortGaussians( GaussBGPoint* g_point, double* sort_key,                 //�Ը�˹�ֲ��ľ�ֵ�뷽��ı�ֵ���н�������
							 GaussBGStatModelParams *bg_model_params );
void  BackgroundTest( const int nChannels, int n, int p, int *match,                   //��⵱ǰ�����ͼ������ֵΪǰ��ֵ���Ǳ���ֵ�������Ʊ���
					 GaussBGModel* bg_model );

/************************************************************************/
/*                    B.���ڱ�Ҷ˹ģ�͵ı�����ģ���˶�Ŀ����          */
/************************************************************************/

/************************************************************************/
/*                          ��Ҫ��ͽṹ�������                        */
/************************************************************************/
#define  LC              128
#define  N1C             15
#define  N2C             25

#define  LCC             64
#define  N1CC            25
#define  N2CC            40

/* BG reference image update parameter */
#define  ALPHA_1         0.02f

/* stat model update parameter
0.002f ~ 1K frame(~45sec), 0.005 ~ 18sec (if 25fps and absolutely static BG) */
#define  ALPHA_2         0.005f

/* start value for alpha parameter (to fast initiate statistic model) */
#define  ALPHA_3         0.1f

#define  DELTA           1

#define  FGD_T           0.9f

typedef struct FGDStatModelParams
{
    int           Lc, N1c, N2c, Lcc, N1cc, N2cc;
    float         alpha1, alpha2, alpha3, delta, T;
}
FGDStatModelParams;

typedef struct BGPixelCStatTable
{
    float          Pv, Pvb;
    unsigned char  v[3];
}
BGPixelCStatTable;

typedef struct BGPixelCCStatTable
{
    float          Pv, Pvb;
    unsigned char  v[6];
}
BGPixelCCStatTable;

typedef struct BGPixelStat
{
    float                 Pbc;
    float                 Pbcc;
    BGPixelCStatTable*    ctable;
    BGPixelCCStatTable*   cctable;
    unsigned char         is_trained_st_model;
    unsigned char         is_trained_dyn_model;
}
BGPixelStat;

typedef struct FGDStatModel
{                                                   
	unsigned char*       background;   //3ͨ��
	unsigned char*       foreground;   //1ͨ��
	
    BGPixelStat*         pixel_stat;
    unsigned char*       Ftd;//1ͨ��
    unsigned char*       Fbd;//1ͨ��
    unsigned char*       prev_frame; //3ͨ��
    FGDStatModelParams   params;
	int width;
	int height;
}
FGDStatModel;
/************************************************************************/
/*                          ��Ҫ����������                              */
/************************************************************************/             
FGDStatModel* CreateFGDStatModel( unsigned char*Image,int Width,int Height,             //��ʼ�洴����Ҷ˹����ͳ��ģ����ز���
								 int nChannels, FGDStatModelParams* parameters );    
void UpdateFGDStatModel( unsigned char*Image,int Width,int Height,int nChannels,        //ͨ��������ͼ��ʱ�ĸ��µ�ǰ��Ҷ˹ģ��
						FGDStatModel*  model );
void ReleaseFGDStatModel( FGDStatModel* _model );                                       //�������ͷű�Ҷ˹ģ�������Դ

/************************************************************************/
/*                    C.���ڼ򵥱�����ģ���˶�Ŀ����                  */
/************************************************************************/

/************************************************************************/
/*                          ��Ҫ��ͽṹ�������                        */
/************************************************************************/
enum { SIMPLEBG_AVG,SIMPLEBG_MEDIAN ,SIMPLEBG_HISTPEAK};	                            //�򵥱�����ģ���͵�ö��
typedef struct SimpleBGModel
{                                                   
	unsigned char*   background;   //
	unsigned char*   foreground;   //1ͨ��
	short int*       pBackHist;
	
	int              SimBGType;    //�򵥱�����ģ������ 1: ������ֵ�˲� 2��������ֵ�˲� 3������ֱ��ͼ
    int              countFrames;  //�Ѿ������֡����
	int              thresh;       //�ָ���ֵ
	int              width;        //ͼ����
	int              height;       //ͼ��߶�
}
SimpleBGModel;

/************************************************************************/
/*                          ��Ҫ����������                              */
/************************************************************************/
SimpleBGModel* CreateSimpleBGModel( unsigned char*Image,int Width,int Height,              //��ʼ֡��ʼ������ģ��
								   int Thresh,int SimBGType);    
void UpdateSimpleModel( unsigned char*Image,int Width,int Height,SimpleBGModel*  model );  //ͨ��������ͼ��ʱ�ĸ��µ�ǰ����ģ��
void ReleaseSimpleModel(SimpleBGModel* _model );                                           //�������ͷű�Ҷ˹ģ�������Դ
void SetBgHistgram(short int *pHistgram, BYTE *sGray,int nWidth, int nHeight);             //��������ʱ��ͳ��ֱ��ͼ
void GetPeakBgByHistgram(short int *pHistgram,BYTE *sBg,int nWidth,int nHeight);           //ͨ��ʱ��ֱ��ͼ������ֵ����
void GetMedianBgByHistgram(short int *pHistgram,BYTE *sBg,int nWidth,                      //ͨ��ʱ��ֱ��ͼ������ֵ����
						   int nHeight,int countFrames);
void GetRuningAvgBg(unsigned char* Bg,unsigned char*Image,int Width,                       //ͨ��ʱ�򱳾���ֵ��������
					int Height,int countFrames);
/************************************************************************/
/*                                 ����                                 */
/************************************************************************/   



/************************************************************************/
/*                    D.���ں��ܶȱ�����ģ���˶�Ŀ����                */
/************************************************************************/

/************************************************************************/
/*                          ��Ҫ��ͽṹ�������                        */
/************************************************************************/

// typedef struct KDEdetectModel
// {
// 	int SampleNum;//��������֡��
//   
// 	double ThreValue;//�ֶ���ֵ����
// 
// }KDEdetectModel;
// 

/************************************************************************/
/*                          ��Ҫ����������                              */
/************************************************************************/


/*************************************************************************
**	�������ƣ�
**		excle()
**	����
**       i   ֵΪi������ 
**	
**	����ֵ��
**		ֵΪi�������ڲ��ұ��ж�Ӧ��ֵ 
**	˵��:
**		�������ұ�ӿ�����ٶ�
************************************************************************/

void excle(double* exclee);
/*************************************************************************
**	�������ƣ�
**		createkdemodel()
**	����
**       SampleImage   ǰ10֡Ϊ����ͼ�����һ֡Ϊ��ǰ�����ͼ�񣬹�11֡
height        ÿһ֡ͼ�������
lineByte      ÿ֡֡ͼ�������
I_sample      ����ͼ��
N_L           ����������֡��+��ǰ֡��N_L֡
exclee        ������
fore          ����Ŀ����ȡ���
**	
**	����ֵ��
**		��
**	˵��:
**		�ú��ܶ�ģ�ͶԵ�ǰ֡��������ȡ��Ŀ��
************************************************************************/
void createkdemodel(vector<unsigned char*> &SampleImage,
					int height,
					int lineByte,
					int N_L,
					double* exclee,
					unsigned char*fore,
					double thre);				   






/************************************************************************/
/*                                 ����                                 */
/************************************************************************/   
#endif 
 
 
 