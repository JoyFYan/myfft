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
/*                    四、红外目标检测模块                              */
/*                                                                      */
/************************************************************************/

/************************************************************************/
/*                    4.2 扩展目标检测相关函数(二)                      */
/************************************************************************/ 

/************************************************************************/
/*                    A.基于多高斯背景建模的运动目标检测                */
/************************************************************************/ 

/************************************************************************/
/*                          主要宏和结构体的声明                        */
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

//高斯背景模型参数结构体
typedef struct GaussBGStatModelParams
{    
    int     win_size;                        /* = 1/alpha */
    int     n_gauss;                         //高斯分布的个数
    double  bg_threshold, std_threshold;
    double  weight_init, variance_init;
}GaussBGStatModelParams;

//高斯背景模型中的高斯分布参数结构体
typedef struct GaussBGValues
{
    int         match_sum;                   //统计高斯分布被匹配的次数
    double      weight;                      //高斯分布的权重
    double      variance[BGFG_MOG_NCOLORS];  //保存高斯分布的方差
    double      mean[BGFG_MOG_NCOLORS];      //保存高斯分布的均值
}
GaussBGValues;

//高斯背景模型像素点
typedef struct GaussBGPoint
{
    GaussBGValues* g_values;                  //高斯分布函数的指针（可以为多个高斯分布）
}
GaussBGPoint;

//高斯背景结构体的定义
typedef struct GaussBGModel
{
	unsigned char*  background;                //估计出的背景值 
	unsigned char*  foreground;                //二值前景图(1通道)                        
    GaussBGStatModelParams   params;           //高斯背景模型参数
    GaussBGPoint*            g_point;          //高斯背景像素
    int                      countFrames;      //已经处理的帧数计数
	int width;                                 //图像的宽度
	int height;                                //图像的帧数
}
GaussBGModel;
/************************************************************************/
/*                          主要函数的声明                              */
/************************************************************************/ 
GaussBGModel* CreateGaussianBGModel( unsigned char*Image,int Width,int Height,         //初始帧创建高斯背景模型     
						int nChannels, GaussBGStatModelParams* parameters );

void  UpdateGaussianBGModel(  unsigned char*Image,int Width,int Height,                //通过新输入图像即时更新高斯背景模型
							int nChannels, GaussBGModel*  bg_model );

void  ReleaseGaussianBGModel( GaussBGModel** _bg_model );                              //检测结束释放高斯模型相关资源
            
int   MatchTest( double* src_pixel, int nChannels, int* match,                         //匹配检测
	   const GaussBGPoint* g_point,const GaussBGStatModelParams *bg_model_params );

void  UpdateFullWindow( double* src_pixel, int nChannels, int* match,                  //全窗情况下而且有高斯分布被匹配，高斯分布各个参数的更新
		GaussBGPoint* g_point,const GaussBGStatModelParams *bg_model_params );

void  UpdatePartialWindow( double* src_pixel, int nChannels, int* match,               //部分窗情况下而且有高斯分布被匹配，高斯分布各个参数的更新
		GaussBGPoint* g_point, const GaussBGStatModelParams *bg_model_params );

void  UpdateFullNoMatch( unsigned char* gm_image, int nChannels,int p, int* match,     //全窗情况下没有高斯分布被匹配，高斯分布各个参数的更新
		GaussBGPoint* g_point,const GaussBGStatModelParams *bg_model_params);

void  UpdatePartialNoMatch(double *pixel,int nChannels,int* /*match*/,                 //部分窗情况下没有高斯分布被匹配，高斯分布各个参数的更新
			GaussBGPoint* g_point,const GaussBGStatModelParams *bg_model_params);

void  GetSortKey( const int nChannels, double* sort_key, const GaussBGPoint* g_point,  //计算高斯分布的均值与方差的比值
			const GaussBGStatModelParams *bg_model_params );

void  InsertionSortGaussians( GaussBGPoint* g_point, double* sort_key,                 //对高斯分布的均值与方差的比值进行降序排列
							 GaussBGStatModelParams *bg_model_params );
void  BackgroundTest( const int nChannels, int n, int p, int *match,                   //检测当前处理的图像像素值为前景值还是背景值，并估计背景
					 GaussBGModel* bg_model );

/************************************************************************/
/*                    B.基于贝叶斯模型的背景建模的运动目标检测          */
/************************************************************************/

/************************************************************************/
/*                          主要宏和结构体的声明                        */
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
	unsigned char*       background;   //3通道
	unsigned char*       foreground;   //1通道
	
    BGPixelStat*         pixel_stat;
    unsigned char*       Ftd;//1通道
    unsigned char*       Fbd;//1通道
    unsigned char*       prev_frame; //3通道
    FGDStatModelParams   params;
	int width;
	int height;
}
FGDStatModel;
/************************************************************************/
/*                          主要函数的声明                              */
/************************************************************************/             
FGDStatModel* CreateFGDStatModel( unsigned char*Image,int Width,int Height,             //初始真创建贝叶斯背景统计模型相关参数
								 int nChannels, FGDStatModelParams* parameters );    
void UpdateFGDStatModel( unsigned char*Image,int Width,int Height,int nChannels,        //通过新输入图像即时的更新当前贝叶斯模型
						FGDStatModel*  model );
void ReleaseFGDStatModel( FGDStatModel* _model );                                       //检测结束释放贝叶斯模型相关资源

/************************************************************************/
/*                    C.基于简单背景建模的运动目标检测                  */
/************************************************************************/

/************************************************************************/
/*                          主要宏和结构体的声明                        */
/************************************************************************/
enum { SIMPLEBG_AVG,SIMPLEBG_MEDIAN ,SIMPLEBG_HISTPEAK};	                            //简单背景建模类型的枚举
typedef struct SimpleBGModel
{                                                   
	unsigned char*   background;   //
	unsigned char*   foreground;   //1通道
	short int*       pBackHist;
	
	int              SimBGType;    //简单背景建模的类型 1: 背景均值滤波 2：背景中值滤波 3：背景直方图
    int              countFrames;  //已经处理的帧计数
	int              thresh;       //分割阈值
	int              width;        //图像宽度
	int              height;       //图像高度
}
SimpleBGModel;

/************************************************************************/
/*                          主要函数的声明                              */
/************************************************************************/
SimpleBGModel* CreateSimpleBGModel( unsigned char*Image,int Width,int Height,              //初始帧初始化背景模型
								   int Thresh,int SimBGType);    
void UpdateSimpleModel( unsigned char*Image,int Width,int Height,SimpleBGModel*  model );  //通过新输入图像即时的更新当前背景模型
void ReleaseSimpleModel(SimpleBGModel* _model );                                           //检测结束释放贝叶斯模型相关资源
void SetBgHistgram(short int *pHistgram, BYTE *sGray,int nWidth, int nHeight);             //建立背景时域统计直方图
void GetPeakBgByHistgram(short int *pHistgram,BYTE *sBg,int nWidth,int nHeight);           //通过时域直方图建立峰值背景
void GetMedianBgByHistgram(short int *pHistgram,BYTE *sBg,int nWidth,                      //通过时域直方图建立中值背景
						   int nHeight,int countFrames);
void GetRuningAvgBg(unsigned char* Bg,unsigned char*Image,int Width,                       //通过时域背景均值建立背景
					int Height,int countFrames);
/************************************************************************/
/*                                 结束                                 */
/************************************************************************/   



/************************************************************************/
/*                    D.基于核密度背景建模的运动目标检测                */
/************************************************************************/

/************************************************************************/
/*                          主要宏和结构体的声明                        */
/************************************************************************/

// typedef struct KDEdetectModel
// {
// 	int SampleNum;//样本所含帧数
//   
// 	double ThreValue;//手动阈值设置
// 
// }KDEdetectModel;
// 

/************************************************************************/
/*                          主要函数的声明                              */
/************************************************************************/


/*************************************************************************
**	函数名称：
**		excle()
**	参数
**       i   值为i的像素 
**	
**	返回值：
**		值为i的像素在查找表中对应的值 
**	说明:
**		建立查找表加快计算速度
************************************************************************/

void excle(double* exclee);
/*************************************************************************
**	函数名称：
**		createkdemodel()
**	参数
**       SampleImage   前10帧为样本图像最后一帧为当前处理的图像，共11帧
height        每一帧图像的行数
lineByte      每帧帧图像的列数
I_sample      样本图像
N_L           样本所含的帧数+当前帧共N_L帧
exclee        索引表
fore          最后的目标提取结果
**	
**	返回值：
**		无
**	说明:
**		用核密度模型对当前帧做处理，提取出目标
************************************************************************/
void createkdemodel(vector<unsigned char*> &SampleImage,
					int height,
					int lineByte,
					int N_L,
					double* exclee,
					unsigned char*fore,
					double thre);				   






/************************************************************************/
/*                                 结束                                 */
/************************************************************************/   
#endif 
 
 
 