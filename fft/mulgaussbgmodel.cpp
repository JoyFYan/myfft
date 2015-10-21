// MulGaussBGModel.cpp: implementation of the CMulGaussBGModel class.
//
//////////////////////////////////////////////////////////////////////

#include "../main/stdafx.h"
#include "../main/IFTrack.h"
#include "MulGaussBGModel.h"
#include "stdlib.h"
#include "math.h"
#include <memory.h>
#include <assert.h>
#include "../os/PublicFunction.h"

#include "../preproc/preproc.h"
#include "../segmentation/segmentation.h"
#include "../feature/feature.h"
#include "../detect/detect.h"
#include "../target/target.h"
#include "../track/track.h"
#include "../track/trackstrategy.h"

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif

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

/*************************************************************************
**	函数名称：
**		CreateGaussianBGModel()
**	参数
**      unsigned char*Image   初始帧图像
**		int Width             图像的宽度
**		int Height            图像的高度
**      int nChannels         图像的通道数
**		GaussBGStatModelParams* parameters  高斯背景模型参数
**	返回值：
**		高斯背景统计模型 
**	说明:
**		初始时创建高斯背景统计模型相关参数
************************************************************************/
GaussBGModel* CreateGaussianBGModel( unsigned char*Image,int Width,int Height,int nChannels, GaussBGStatModelParams* parameters )
{
    GaussBGModel* bg_model = 0;
    double var_init;
    GaussBGStatModelParams params;
    int i, j, k, n, m, p;
    
    //init parameters
    if(parameters == NULL)
    {
        params.win_size = BGFG_MOG_WINDOW_SIZE;                        // 初始化阶段的帧数；用户自定义模型学习率a=1/win_size;
        params.bg_threshold = BGFG_MOG_BACKGROUND_THRESHOLD;
        params.std_threshold = BGFG_MOG_STD_THRESHOLD;
        params.weight_init = BGFG_MOG_WEIGHT_INIT;
        params.variance_init = BGFG_MOG_SIGMA_INIT*BGFG_MOG_SIGMA_INIT;//方差
        params.n_gauss = BGFG_MOG_NGAUSSIANS;                          //高斯分布函数的个数
    }
    else
    {
        params = *parameters;                                           //用户自定义参数
    }
    
	bg_model = (GaussBGModel*)malloc(sizeof(*bg_model));
    memset( bg_model, 0, sizeof(*bg_model) );
    bg_model->params = params;
	bg_model->width = Width;
	bg_model->height = Height;
    
    //prepare storages
	bg_model->g_point = (GaussBGPoint*)malloc(sizeof(GaussBGPoint)*(Width*Height));
	for (i=0;i<Width*Height;i++)
		bg_model->g_point[i].g_values = (GaussBGValues*)malloc(sizeof(GaussBGValues)*params.n_gauss);
	bg_model->background = (unsigned char*)malloc(Width*Height*nChannels);
	bg_model->foreground = (unsigned char*)malloc(Width*Height*1);
    //initializing
    var_init = 2 * params.std_threshold * params.std_threshold;            //初始化方差
   
    for( i = 0, p = 0, n = 0; i < Height; i++ )        //行
    {
        for( j = 0; j < Width; j++, n++ )               //列
        {
			//以下几步是对第一个高斯函数做初始化
            bg_model->g_point[n].g_values[0].weight = 1;     //权值赋为1
            bg_model->g_point[n].g_values[0].match_sum = 1;  //高斯函数被匹配的次数
            for( m = 0; m < nChannels; m++)
            {
                bg_model->g_point[n].g_values[0].variance[m] = var_init;
				//均值赋为当前像素的值
                bg_model->g_point[n].g_values[0].mean[m] = (unsigned char)Image[p + m];
            }

			//除第一以外的高斯分布函数的初始化(均值、权值和匹配次数都置零)
            for( k = 1; k < params.n_gauss; k++)
            {
                bg_model->g_point[n].g_values[k].weight = 0;
                bg_model->g_point[n].g_values[k].match_sum = 0;
                for( m = 0; m < nChannels; m++){
                    bg_model->g_point[n].g_values[k].variance[m] = var_init;
                    bg_model->g_point[n].g_values[k].mean[m] = 0;
                }
            }
            p += nChannels;
        }
    }//g_point[]:像素，g_values[]:高斯分布函数，mean[]:通道

    bg_model->countFrames = 0; //设置初始帧为0
	return bg_model;
}

/*************************************************************************
**	函数名称：
**		UpdateGaussianBGModel()
**	参数
**      unsigned char*Image      当前要处理的图像
**		int Width                图像的宽度
**		int Height               图像的高度
**      int nChannels            图像的通道数
**		GaussBGModel*  bg_model  高斯背景模型
**	返回值：
**		无
**	说明:
**		通过新输入图像即时的更新当前高斯背景模型，背景模型的更新，
**      不仅要更新高斯分布函数的参数，还要更新各高斯函数的权重
************************************************************************/
void UpdateGaussianBGModel(  unsigned char*Image,int Width,int Height,int nChannels, GaussBGModel*  bg_model )
{
    int i, j, k;
    bg_model->countFrames++;
    
    for( i = 0; i < Height; i++ )
    {
        for( j = 0; j < Width; j++ )
        {
            int match[BGFG_MOG_MAX_NGAUSSIANS];                //对高斯函数做标记，match[m]=1表示函数m为匹配的高斯分布函数
            double sort_key[BGFG_MOG_MAX_NGAUSSIANS];          //此数组存贮每个高斯函数的均值与方差比值
            const int n = i*Width+j;
            const int p = n*nChannels;
            
            // A few short cuts
            GaussBGPoint* g_point = &bg_model->g_point[n];
            const GaussBGStatModelParams bg_model_params = bg_model->params;
            double pixel[4];                                 // pixel[]存贮当前像素的各通道RGB值
            int no_match;
            
            for( k = 0; k < nChannels; k++ )
                pixel[k] = (unsigned char)Image[p+k];

            //检查是否有与当前像素匹配的高斯函数
            no_match = MatchTest( pixel, nChannels, match, g_point, &bg_model_params );
            if( bg_model->countFrames == bg_model->params.win_size )   //参数学习因子自定义                 
            {
                UpdateFullWindow( pixel, nChannels, match, g_point, &bg_model->params );//存在匹配函数时，相关高斯参数的更新
                if( no_match == -1)        //没有高斯函数匹配时，增加新的高斯函数
                    UpdateFullNoMatch( Image, nChannels, p, match, g_point, &bg_model_params );
            }
            else                                                       //参数学习因子取决于高斯分布的匹配次数
            {
                UpdatePartialWindow( pixel, nChannels, match, g_point, &bg_model_params );//存在匹配函数时，相关高斯参数的更新
                if( no_match == -1)       //没有高斯函数匹配时，增加新的高斯函数
                    UpdatePartialNoMatch( pixel, nChannels, match, g_point, &bg_model_params );
            }
			//计算高斯分布的均值和方差比值
            GetSortKey( nChannels, sort_key, g_point, &bg_model_params );
			//对高斯分布的均值和方差比值进行降序排列
            InsertionSortGaussians( g_point, sort_key, (GaussBGStatModelParams *)&bg_model_params );
			//背景检测，前景的更新
            BackgroundTest( nChannels, n, p, match, bg_model );
        }
    }
}
/*************************************************************************
**	函数名称：
**		ReleaseGaussianBGModel()
**	参数
**		GaussBGModel** _bg_model  高斯背景模型
**	返回值：
**		无 
**	说明:
**		检测结束释放高斯背景模型相关资源
************************************************************************/
void ReleaseGaussianBGModel( GaussBGModel** _bg_model )
{
    if( *_bg_model )
    {
        GaussBGModel* bg_model = *_bg_model;
        if( bg_model->g_point )
        {
			int nValues = bg_model->width*bg_model->height;
			for (int i=0;i<nValues;i++)
				free(bg_model->g_point[i].g_values);
			free(bg_model->g_point);
        }
		free(bg_model->background);
        free(bg_model->foreground);
        memset( bg_model, 0, sizeof(*bg_model) );
		free(*_bg_model);
    }
}
/*************************************************************************
**	函数名称：
**		MatchTest()
**	参数
**      double* src_pixel                              即当前帧中该像素的RGB值
**		int nChannels                                  图像的通道数
**      int* match                                     当前处理的图像像素值的匹配值指针
**		GaussBGPoint* g_point                          高斯背景点指针
**		const GaussBGStatModelParams *bg_model_params  高斯背景模型参数
**	返回值：
**		为-1则没有匹配的位置，否则返回高斯匹配函数的位置
**	说明:
**		拿当前像素的值与已存在的高斯分布函数比较，查找是否存在匹配的
**      的高斯分布函数，如果有则返回 k值(高斯分布函数的序号)
************************************************************************/
int MatchTest( double* src_pixel, int nChannels, int* match,
						const GaussBGPoint* g_point,
						const GaussBGStatModelParams *bg_model_params)
{
    int k;
    int matchPosition=-1;
    for ( k = 0; k < bg_model_params->n_gauss; k++) match[k]=0;
    
    for ( k = 0; k < bg_model_params->n_gauss; k++) 
	{
        double sum_d2 = 0.0;
        double var_threshold = 0.0;
        for(int m = 0; m < nChannels; m++)
		{
            double d = g_point->g_values[k].mean[m]- src_pixel[m];//通道m的原始模型值与当前像素的值之差
            sum_d2 += (d*d);
            var_threshold += g_point->g_values[k].variance[m];
        }  //difference < STD_LIMIT*STD_LIMIT or difference**2 < STD_LIMIT*STD_LIMIT*VAR

		//当前sum_d2为d0,d1,d2的平方和，var_threshold的值为像素各通道方差之和

        var_threshold = bg_model_params->std_threshold*bg_model_params->std_threshold*var_threshold;

        if(sum_d2 < var_threshold)    //查看是否可以与某高斯分布匹配
		{
            match[k] = 1;
            matchPosition = k;
            break;                   //如果和第k个高斯函数匹配，则终止与后续函数的匹配
        }
    }
    
    return matchPosition;
}
/*************************************************************************
**	函数名称：
**		UpdateFullWindow()
**	参数
**      double* src_pixel                              即当前帧中该像素的RGB值
**		int nChannels                                  图像的通道数
**      int* match                                     当前处理的图像像素值的匹配值指针
**		GaussBGPoint* g_point                          高斯背景点指针
**		const GaussBGStatModelParams *bg_model_params  高斯背景模型参数
**	返回值：
**		无
**	说明:
**		更新每个高斯分布的权值(对匹配的高斯函数k加大权值，其余的则减小权值)，如果前面的结果中
**      存在匹配的高斯分布函数k,则需要再对第k个高斯分布函数的均值mean和方差variance做修正
************************************************************************/
void UpdateFullWindow( double* src_pixel, int nChannels, int* match,
					  GaussBGPoint* g_point,
					  const GaussBGStatModelParams *bg_model_params )
{
	//用户自定义模型学习率a
    const double learning_rate_weight = (1.0/(double)bg_model_params->win_size);

    for(int k = 0; k < bg_model_params->n_gauss; k++)
	{
		//对每个高斯分布的权值做修正：w=(1-a)w+a*m (a:模型学习率，m是匹配，匹配就是1，不匹配就是0)
        g_point->g_values[k].weight = g_point->g_values[k].weight +
            (learning_rate_weight*((double)match[k] -g_point->g_values[k].weight));

		//如果存在匹配的高斯分布函数k（当前像素为背景像素）,则需要再对第k个高斯分布函数的均值mean和方差variance更新
        if(match[k])
		{
            double learning_rate_gaussian = (double)match[k]/(g_point->g_values[k].weight*
                (double)bg_model_params->win_size);                 //参数学习率p(p=a/w)
            for(int m = 0; m < nChannels; m++)
			{
				//参数更新公式：u=(1-p)*u0+p*x; o*o=(1-p)*o*o+p*tmpDiff*tmpDiff
				//当前像素的通道m的值与原始模型值之差
                const double tmpDiff = src_pixel[m] - g_point->g_values[k].mean[m];
                g_point->g_values[k].mean[m] = g_point->g_values[k].mean[m] +
                    (learning_rate_gaussian * tmpDiff);
                g_point->g_values[k].variance[m] = g_point->g_values[k].variance[m]+
                    (learning_rate_gaussian*((tmpDiff*tmpDiff) - g_point->g_values[k].variance[m]));
            }
        }
    }
}
/*************************************************************************
**	函数名称：
**		UpdatePartialWindow()
**	参数
**      double* src_pixel                              即当前帧中该像素的RGB值
**		int nChannels                                  图像的通道数
**      int* match                                     当前处理的图像像素值的匹配值指针
**		GaussBGPoint* g_point                          高斯背景点指针
**		const GaussBGStatModelParams *bg_model_params  高斯背景模型参数
**	返回值：
**		无
**	说明:
**		对所有的高斯分布函数做更新.至少每个高斯分布的权值必须修正，如果前面的结果中存在匹配的高斯分布函数k,
**      则需要再对第k个高斯分布函数的match_sum修改，最终对那些匹配的高斯分布函数k的参数match_sum>0的做均值mean和方差variance修正
************************************************************************/
void UpdatePartialWindow( double* src_pixel, int nChannels, int* match, GaussBGPoint* g_point, const GaussBGStatModelParams *bg_model_params )
{
    int k, m;
    int window_current = 0;
    
    for( k = 0; k < bg_model_params->n_gauss; k++ )
        window_current += g_point->g_values[k].match_sum;//window_current为k个高斯分布函数的match_sum值之和
    
    for( k = 0; k < bg_model_params->n_gauss; k++ )
    {
        g_point->g_values[k].match_sum += match[k];      //修正匹配的高斯分布函数k的match_sum值
        double learning_rate_weight = (1.0/((double)window_current + 1.0)); //increased by one since sum

		//修正每个高斯分布的权值
        g_point->g_values[k].weight = g_point->g_values[k].weight +
            (learning_rate_weight*((double)match[k] - g_point->g_values[k].weight));

        //如果存在匹配的高斯分布函数k（当前像素为背景像素）,则需要再对第k个高斯分布函数的均值mean和方差variance更新
        if( g_point->g_values[k].match_sum > 0 && match[k] )
        {
            double learning_rate_gaussian = (double)match[k]/((double)g_point->g_values[k].match_sum);//参数学习率p(p=a/w)
            for( m = 0; m < nChannels; m++ )
            {
				//参数更新公式：u=(1-p)*u0+p*x; o*o=(1-p)*o*o+p*tmpDiff*tmpDiff
				//当前像素的通道m的值与原始模型值之差
                const double tmpDiff = src_pixel[m] - g_point->g_values[k].mean[m];
                g_point->g_values[k].mean[m] = g_point->g_values[k].mean[m] +
                    (learning_rate_gaussian*tmpDiff);
                g_point->g_values[k].variance[m] = g_point->g_values[k].variance[m]+
                    (learning_rate_gaussian*((tmpDiff*tmpDiff) - g_point->g_values[k].variance[m]));
            }
        }
    }
}
/*************************************************************************
**	函数名称：
**		UpdateFullNoMatch()
**	参数
**      unsigned char* gm_image                        当前要处理的图像
**		int nChannels                                  图像的通道数
**      int p                                          图像像素点的位置
**      int* match                                     当前处理的图像像素值的匹配值指针
**		GaussBGPoint* g_point                          高斯背景点指针
**		const GaussBGStatModelParams *bg_model_params  高斯背景模型参数
**	返回值：
**		无
**	说明:
**		当所有的高斯函数均不匹配时，说明有新的分布出现，需要将原高斯函数中sort_key最小的
**      替换为新的高斯函数(权值小，方差大),其余的高斯函数对应的只需更新权值
************************************************************************/
void UpdateFullNoMatch( unsigned char* gm_image, int nChannels, int p, int* match,
								 GaussBGPoint* g_point,
								 const GaussBGStatModelParams *bg_model_params)
{
    int k, m;
    double alpha;
    int match_sum_total = 0;
	
    //最后一个高斯分布设置为新的高斯函数
    g_point->g_values[bg_model_params->n_gauss - 1].match_sum = 1;  //将新的高斯分布函数的match_sum置为1

    
    //计算各个高斯分布match_sum 之和
    for( k = 0; k < bg_model_params->n_gauss ; k++ )
        match_sum_total += g_point->g_values[k].match_sum;
    //给新的高斯分布函数赋一个较小的权值
    g_point->g_values[bg_model_params->n_gauss - 1].weight = 1./(double)match_sum_total;

	//将新的高斯分布函数的variance[m]全部置为variance_init；mean[m]的值置为当前像素各通道的值
    for( m = 0; m < nChannels ; m++ )
    {
        g_point->g_values[bg_model_params->n_gauss - 1].variance[m] = bg_model_params->variance_init;
        g_point->g_values[bg_model_params->n_gauss - 1].mean[m] = (unsigned char)gm_image[p + m];
    }

    //对其他的高斯分布函数做权值更新：w=(1-a)*w+a*m (a:模型学习率，m是匹配，匹配就是1，不匹配就是0)
    alpha = 1.0 - (1.0/bg_model_params->win_size);
    for( k = 0; k < bg_model_params->n_gauss - 1; k++ )
    {
        g_point->g_values[k].weight *= alpha;   
        if( match[k] )                          //如果匹配的话 w=(1-a)*w+a*m 加上后半段
            g_point->g_values[k].weight += alpha;
    }
}
/*************************************************************************
**	函数名称：
**		UpdatePartialNoMatch()
**	参数
**      double *pixel                                  当前要处理的图像像素值
**		int nChannels                                  图像的通道数
**      int* match                                     当前处理的图像像素值的匹配值指针
**		GaussBGPoint* g_point                          高斯背景点指针
**		const GaussBGStatModelParams *bg_model_params  高斯背景模型参数
**	返回值：
**		无
**	说明:
**		当所有的高斯函数均不匹配时，需要将原高斯函数中sort_key最小的
**      替换为新的高斯函数(权值小，方差大),然后归一化所有的高斯分布权值
************************************************************************/
void UpdatePartialNoMatch(double *pixel,
                        int nChannels,
                        int* /*match*/,
                        GaussBGPoint* g_point,
                        const GaussBGStatModelParams *bg_model_params)
{
    int k, m;
    //将最后高斯分布函数的match_sum置为1
    g_point->g_values[bg_model_params->n_gauss - 1].match_sum = 1;
    
    //计算各个高斯分布match_sum 之和
    int match_sum_total = 0;
    for(k = 0; k < bg_model_params->n_gauss ; k++)
        match_sum_total += g_point->g_values[k].match_sum;

	//将最后一个高斯分布函数的variance[m]全部置为variance_init；mean[m]的值置为当前像素各通道的值
    for(m = 0; m < nChannels; m++)
    {
        g_point->g_values[bg_model_params->n_gauss - 1].variance[m] = bg_model_params->variance_init;
        g_point->g_values[bg_model_params->n_gauss - 1].mean[m] = pixel[m];
    }
	//归一化所有的高斯分布权值
    for(k = 0; k < bg_model_params->n_gauss; k++)
    {
        g_point->g_values[k].weight = (double)g_point->g_values[k].match_sum /
            (double)match_sum_total;
    }
}
/*************************************************************************
**	函数名称：
**		GetSortKey()
**	参数
**      const int nChannels                            图像的通道数
**		double* sort_key                               weight/sqrt(variance_sum)的值的指针
**		const GaussBGPoint* g_point                    高斯背景点指针
**		const GaussBGStatModelParams *bg_model_params  高斯背景模型参数
**	返回值：
**		无
**	说明:
**		计算各个高斯分布weight/sqrt(variance_sum)的值，
**      后面将对该值进行排序（该值越大则表示背景的可能性就越大）
************************************************************************/
void GetSortKey( const int nChannels, double* sort_key, const GaussBGPoint* g_point,
						  const GaussBGStatModelParams *bg_model_params )
{
    int k, m;
    for( k = 0; k < bg_model_params->n_gauss; k++ )
    {
        // 避免被零除
        if( g_point->g_values[k].match_sum > 0 )
        {
            // 假设高斯成分各个部分是独立的
            double variance_sum = 0.0;
            for( m = 0; m < nChannels; m++ )                              //计算各通道方差和
                variance_sum += g_point->g_values[k].variance[m];
            
            sort_key[k] = g_point->g_values[k].weight/sqrt(variance_sum); //计算高斯分布weight/sqrt(variance_sum)
        }
        else
            sort_key[k]= 0.0;
    }
}
/*************************************************************************
**	函数名称：
**		InsertionSortGaussians()
**	参数
**		const GaussBGPoint* g_point                    高斯背景点指针
**		double* sort_key                               weight/sqrt(variance_sum)的值的指针
**		GaussBGStatModelParams *bg_model_params  高斯背景模型参数
**	返回值：
**		无
**	说明:
**		对各个高斯分布weight/sqrt(variance_sum)的值进行插入排序，降序排列
************************************************************************/
void InsertionSortGaussians( GaussBGPoint* g_point, double* sort_key, GaussBGStatModelParams *bg_model_params )
{
    int i, j;
    for( i = 1; i < bg_model_params->n_gauss; i++ )
    {
        double index = sort_key[i];
        for( j = i; j > 0 && sort_key[j-1] < index; j-- )           //降序排列
        {
            double temp_sort_key = sort_key[j];                     //交换两个sortkey的值
            sort_key[j] = sort_key[j-1];
            sort_key[j-1] = temp_sort_key;
            
            GaussBGValues temp_gauss_values = g_point->g_values[j]; //交换sortkey对应的高斯分布的值
            g_point->g_values[j] = g_point->g_values[j-1];
            g_point->g_values[j-1] = temp_gauss_values;
        }
    }
}
/*************************************************************************
**	函数名称：
**		BackgroundTest()
**	参数
**      const int nChannels     图像的通道数
**		int n                   高斯背景点位置 i×width+j
**		int p                   图像像素的位置 n×nChannels
**      int *match              图像像素值的匹配值指针
**		GaussBGModel* bg_model  高斯背景模型
**	返回值：
**		无
**	说明:
**		检测当前处理的图像像素值为前景值还是背景值
************************************************************************/
void BackgroundTest( const int nChannels, int n, int p, int *match, GaussBGModel* bg_model )
{
    int m, b;
    unsigned char pixelValue = (unsigned char)255; // 初始值设为255，找到高斯匹配函数，再将其设置为0
    double weight_sum = 0.0;
    GaussBGPoint* g_point = bg_model->g_point;
    
    for( m = 0; m < nChannels; m++)                 //得到高斯模型的估计背景
        bg_model->background[p+m]   = (unsigned char)(g_point[n].g_values[0].mean[m]+0.5);
    
    for( b = 0; b < bg_model->params.n_gauss; b++)  //与多个高斯函数进行匹配，检测该像素是否匹配
    {
        weight_sum += g_point[n].g_values[b].weight;
        if( match[b] )                              //如果为真，说明该像素已与某高斯函数匹配，该像素为背景
            pixelValue = 0;
        if( weight_sum > bg_model->params.bg_threshold )
            break;                                  //如果if语句为真，则前b个高斯分布被选为描述背景的函数
    }
    
    bg_model->foreground[p/nChannels] = pixelValue; //设置前景值
}
/************************************************************************/
/*                    B.基于贝叶斯模型的背景建模的运动目标检测          */
/************************************************************************/ 

/*************************************************************************
**	函数名称：
**		CreateFGDStatModel()
**	参数
**      unsigned char*Image   初始帧图像
**		int Width             图像的宽度
**		int Height            图像的高度
**      int nChannels         图像的通道数
**		FGDStatModelParams* parameters  贝叶斯背景模型参数
**	返回值：
**		贝叶斯背景统计模型 
**	说明:
**		初始时创建贝叶斯背景统计模型相关参数
************************************************************************/
FGDStatModel* CreateFGDStatModel( unsigned char*Image,int Width,int Height,int nChannels, FGDStatModelParams* parameters )
{
    FGDStatModel* p_model = 0;
    
    int i,pixel_count, buf_size;
    FGDStatModelParams params;
	
    //init parameters
    if( parameters == NULL )
    {
        params.Lc = LC;
        params.N1c = N1C;
        params.N2c = N2C;
        params.Lcc = LCC;
        params.N1cc = N1CC;
        params.N2cc = N2CC;
        params.delta = DELTA;
        params.alpha1 = ALPHA_1;
        params.alpha2 = ALPHA_2;
        params.alpha3 = ALPHA_3;
        params.T = FGD_T;
    }
    else
    {
        params = *parameters;
    }
    p_model = (FGDStatModel*)malloc( sizeof(*p_model) );
    memset( p_model, 0, sizeof(*p_model) );
    p_model->params = params;
	p_model->width = Width;
	p_model->height = Height;
    //init storages
    pixel_count = Width*Height;
    
    buf_size = pixel_count*sizeof(p_model->pixel_stat[0]);
    p_model->pixel_stat = (BGPixelStat*)malloc(buf_size) ;
    memset( p_model->pixel_stat, 0, buf_size );
    
    buf_size = params.N2c*sizeof(p_model->pixel_stat[0].ctable[0]);
	for (i=0;i<pixel_count;i++)
	{
		p_model->pixel_stat[i].ctable =(BGPixelCStatTable*)malloc(buf_size);
		memset( p_model->pixel_stat[i].ctable, 0, buf_size );
	}
    buf_size = params.N2cc*sizeof(p_model->pixel_stat[0].cctable[0]);
	for (i=0;i<pixel_count;i++)
	{
		p_model->pixel_stat[i].cctable =(BGPixelCCStatTable*)malloc(buf_size);
		memset( p_model->pixel_stat[i].cctable, 0, buf_size );
	}
	//init temporary images
	p_model->Ftd = (unsigned char*)malloc(sizeof(unsigned char)*(Width*Height*1));
	p_model->Fbd = (unsigned char*)malloc(sizeof(unsigned char)*(Width*Height*1));
	p_model->foreground = (unsigned char*)malloc(sizeof(unsigned char)*(Width*Height*1));
	
	p_model->background = (unsigned char*)malloc(sizeof(unsigned char)*(Width*Height*nChannels));
	p_model->prev_frame = (unsigned char*)malloc(sizeof(unsigned char)*(Width*Height*nChannels));
	memcpy(p_model->background,Image,sizeof(unsigned char)*(Width*Height*nChannels));
	memcpy(p_model->prev_frame,Image,sizeof(unsigned char)*(Width*Height*nChannels));
	return p_model;
}

/*************************************************************************
**	函数名称：
**		_max_element()
**	参数
**      double* start          序列的头指针
**		double* end            序列的尾指针
**	返回值：
**		最大值指针
**	说明:
**		找出序列的最大值指针
************************************************************************/
static double* _max_element( double* start, double* end )
{
    double* p = start++;
    for( ; start != end; start++ )
        if( *p < *start )
            p = start;
		return p;
}

/*************************************************************************
**	函数名称：
**		ChangeDetection()
**	参数
**      unsigned char * prev_Image            先前帧图像
**		unsigned char * curr_Image            当前要处理的图像
**      unsigned char * changeMask_Image      前景二值图
**      int Width                             图像的宽度
**		int Height                            图像的高度
**      int nChannels                         图像的通道数
**	返回值：
**		无 
**	说明:
**		通过简单的帧差法和简单自适应分割得到前景二值图
************************************************************************/
int ChangeDetection( unsigned char * prev_Image,
				  unsigned char * curr_Image,
				  unsigned char * changeMask_Image,int Width,int Height,int nChannels)
{
	int WidthStep = Width*nChannels;

    int i, j, b, x, y, thres;
    const int PIXELRANGE=256;
	
  //  if(nChannels!=3) return 0;
	
	memset(changeMask_Image,0,sizeof(unsigned char)*Width*Height);
    // All operations per colour
    for (b=0 ; b<nChannels ; b++) {
        // create histogram
		
        long HISTOGRAM[PIXELRANGE]; 
        for (i=0 ; i<PIXELRANGE; i++) HISTOGRAM[i]=0;
        
        for (y=0 ; y<Height ; y++)
        {
            unsigned char* rowStart1 = (unsigned char*)curr_Image + y * WidthStep + b;
            unsigned char* rowStart2 = (unsigned char*)prev_Image + y * WidthStep + b;
            for (x=0 ; x<Width; x++, rowStart1+=nChannels, rowStart2+=nChannels) {
                int diff = abs( int(*rowStart1) - int(*rowStart2) );
                HISTOGRAM[diff]++;
            }
        }
		
        double relativeVariance[PIXELRANGE];
        for (i=0 ; i<PIXELRANGE; i++) relativeVariance[i]=0;
		
        for (thres=PIXELRANGE-2; thres>=0 ; thres--)
        {
            double sum=0;
            double sqsum=0;
            int count=0;
            for (j=thres ; j<PIXELRANGE ; j++) {
                sum   += double(j)*double(HISTOGRAM[j]);
                sqsum += double(j*j)*double(HISTOGRAM[j]);
                count += HISTOGRAM[j];
            }
            count = count == 0 ? 1 : count;
            double my = sum / count;
            double sigma = sqrt( sqsum/count - my*my);
            relativeVariance[thres] = sigma;
        }
        // find maximum
        unsigned char bestThres = 0;
		
        double* pBestThres = _max_element(relativeVariance, relativeVariance+PIXELRANGE);
        bestThres = (unsigned char)(*pBestThres); if (bestThres <10) bestThres=10;
		
        for (y=0 ; y<Height ; y++)
        {
            unsigned char* rowStart1 = (unsigned char*)(curr_Image) + y * WidthStep + b;
            unsigned char* rowStart2 = (unsigned char*)(prev_Image) + y * WidthStep + b;
            unsigned char* rowStart3 = (unsigned char*)(changeMask_Image) + y * Width;
            for (x = 0; x < Width; x++, rowStart1+=nChannels,
                rowStart2+=nChannels, rowStart3+=1) {
                // OR between different color channels
                int diff = abs( int(*rowStart1) - int(*rowStart2) );
                if ( diff > bestThres)
                    *rowStart3 |=255;
            }
        }
    }
	
    return 1;
}

#define MIN_PV 1E-10
#define V_C(k,l) ctable[k].v[l]
#define PV_C(k) ctable[k].Pv
#define PVB_C(k) ctable[k].Pvb
#define V_CC(k,l) cctable[k].v[l]
#define PV_CC(k) cctable[k].Pv
#define PVB_CC(k) cctable[k].Pvb
#define  BGFG_FGD_BG_UPDATE_TRESH 0.5f
/*************************************************************************
**	函数名称：
**		UpdateFGDStatModel()
**	参数
**      unsigned char*Image   当前要处理的图像
**		int Width             图像的宽度
**		int Height            图像的高度
**      int nChannels         图像的通道数
**		FGDStatModel* _model  贝叶斯状态模型
**	返回值：
**		无 
**	说明:
**		通过新输入图像即时的更新当前贝叶斯模型
************************************************************************/
void UpdateFGDStatModel( unsigned char*Image,int Width,int Height,int nChannels, FGDStatModel*  model )
{
	int WidthStep = Width*nChannels;

    int            mask_step = Width;
    unsigned char* prev_frame = model->prev_frame;
    int            FG_pixels_count = 0;
    int            deltaC  = int(model->params.delta * 256 / model->params.Lc);
    int            deltaCC = int(model->params.delta * 256 / model->params.Lcc);
    int            i, j, k, l;

	memset(model->foreground,0,sizeof(unsigned char)*Width*Height);
    //form FG pixels candidates using image differencing with adaptive threshold [P.Rosin, Thresholding for change detection, ICCV, 1998 ]
    ChangeDetection( prev_frame, Image, model->Ftd , Width, Height, nChannels);
    ChangeDetection( model->background, Image, model->Fbd , Width, Height, nChannels);

    for( i = 0; i < Height; i++ )
    {
        for( j = 0; j < Width; j++ )
        {
            if( ((unsigned char*)model->Fbd)[i*mask_step+j] || ((unsigned char*)model->Ftd)[i*mask_step+j] )
            {
                float Pb=0, Pv=0, Pvb=0;
                BGPixelStat* stat = model->pixel_stat + i * Width + j;
                BGPixelCStatTable* ctable = stat->ctable;
                BGPixelCCStatTable* cctable = stat->cctable;
    
                unsigned char* curr_data = (unsigned char*)(Image)+i*WidthStep+j*nChannels;
                unsigned char* prev_data = (unsigned char*)(prev_frame)+i*WidthStep+j*nChannels;

                int val = 0;
                // is it a motion pixel?
                if( ((unsigned char*)model->Ftd)[i*mask_step+j] )
                {
                    if( !stat->is_trained_dyn_model ) val = 1; //没有检测到运动背景之前
                    else
                    {
                        //compare with stored CCt vectors
                        for( k = 0; PV_CC(k) > model->params.alpha2 && k < model->params.N1cc; k++ )
                        {
//                             if ( abs( V_CC(k,0) - prev_data[0]) <=  deltaCC &&
//                                  abs( V_CC(k,1) - prev_data[1]) <=  deltaCC &&
//                                  abs( V_CC(k,2) - prev_data[2]) <=  deltaCC &&
//                                  abs( V_CC(k,3) - curr_data[0]) <=  deltaCC &&
//                                  abs( V_CC(k,4) - curr_data[1]) <=  deltaCC &&
//                                  abs( V_CC(k,5) - curr_data[2]) <=  deltaCC)
							bool condition = true;
							for (int n=0;n<nChannels;n++)
							{
								if (abs( V_CC(k,n) - prev_data[n]) >  deltaCC ||
									abs( V_CC(k,n+3) - curr_data[n]) >  deltaCC)
									condition = false;
							}
							if (condition)
                            {
                                Pv += PV_CC(k);
                                Pvb += PVB_CC(k);
                            }
                        }
                        Pb = stat->Pbcc;
                        if( 2 * Pvb * Pb <= Pv ) val = 1;
                    }
                }
                else if( stat->is_trained_st_model ) //检测到运动背景后
                {
                    //compare with stored Ct vectors
                    for( k = 0; PV_C(k) > model->params.alpha2 && k < model->params.N1c; k++ )
                    {
//                         if ( abs( V_C(k,0) - curr_data[0]) <=  deltaC &&
//                              abs( V_C(k,1) - curr_data[1]) <=  deltaC &&
//                              abs( V_C(k,2) - curr_data[2]) <=  deltaC )
						bool condition = true;
						for (int n=0;n<nChannels;n++)
						{
							if ( abs(V_C(k,n) - curr_data[n] >  deltaC))
								condition = false;
						}
						if (condition)
                        {
                            Pv += PV_C(k);
                            Pvb += PVB_C(k);
                        }
                    }
                    Pb = stat->Pbc;
                    if( 2 * Pvb * Pb <= Pv ) val = 1;
                }
                //update FG
                ((unsigned char*)model->foreground)[i*mask_step+j] = (unsigned char)(val*255);
                FG_pixels_count += val;
            }// end if( change detection...
        }//for j...
    } //for i...
    //end BG/FG classification

    //check ALL BG update condition
    if( ((float)FG_pixels_count/(Width*Height)) > BGFG_FGD_BG_UPDATE_TRESH )
    {
         for( i = 0; i < Height; i++ )
             for( j = 0; j < Width; j++ )
             {
                 BGPixelStat* stat = model->pixel_stat + i * Width + j;
                 stat->is_trained_st_model = stat->is_trained_dyn_model = 1;
             }
    }

    //update BG model 
	//1.更新特征统计表
    for( i = 0; i < Height; i++ )
    {
        for( j = 0; j < Width; j++ )
        {
            BGPixelStat* stat = model->pixel_stat + i * Width + j;
            BGPixelCStatTable* ctable = stat->ctable;
            BGPixelCCStatTable* cctable = stat->cctable;

            unsigned char *curr_data = (unsigned char*)(Image)+i*WidthStep+j*nChannels;
            unsigned char *prev_data = (unsigned char*)(prev_frame)+i*WidthStep+j*nChannels;

            if( ((unsigned char*)model->Ftd)[i*mask_step+j] || !stat->is_trained_dyn_model )
            {
                float alpha = stat->is_trained_dyn_model ? model->params.alpha2 : model->params.alpha3;
                float diff = 0;
                int dist, min_dist = 2147483647, indx = -1;

                //update Pb
                stat->Pbcc *= (1.f-alpha);
                if( !((unsigned char*)model->foreground)[i*mask_step+j] )
                {
                    stat->Pbcc += alpha;
                }

                // find best Vi match
                for(k = 0; PV_CC(k) && k < model->params.N2cc; k++ )
                {
                    // Exponential decay of memory
                    PV_CC(k)  *= (1-alpha);
                    PVB_CC(k) *= (1-alpha);
                    if( PV_CC(k) < MIN_PV )
                    {
                        PV_CC(k) = 0;
                        PVB_CC(k) = 0;
                        continue;
                    }

                    dist = 0;
                    for( l = 0; l < nChannels; l++ )
                    {
                        int val = abs( V_CC(k,l) - prev_data[l] );
                        if( val > deltaCC ) break;
                        dist += val;
                        val = abs( V_CC(k,l+3) - curr_data[l] );
                        if( val > deltaCC) break;
                        dist += val;
                    }
                    if( l == nChannels && dist < min_dist )
                    {
                        min_dist = dist;
                        indx = k;
                    }
                }


                if( indx < 0 )
                {
					//N2th elem in the table is replaced by a new features
                    indx = model->params.N2cc - 1;
                    PV_CC(indx) = alpha;
                    PVB_CC(indx) = alpha;
                    //udate Vt
                    for( l = 0; l < nChannels; l++ )
                    {
                        V_CC(indx,l) = prev_data[l];
                        V_CC(indx,l+3) = curr_data[l];
                    }
                }
                else                        //如果在动态表中找到相似像素，使它相应的概率值增加
                {//update
                    PV_CC(indx) += alpha;
                    if( !((unsigned char*)model->foreground)[i*mask_step+j] )
                    {
                        PVB_CC(indx) += alpha;
                    }
                }

                //re-sort CCt table by Pv  //对动态点表进行排序
                for( k = 0; k < indx; k++ )
                {
                    if( PV_CC(k) <= PV_CC(indx) )
                    {
                        //shift elements
                        BGPixelCCStatTable tmp1, tmp2 = cctable[indx];
                        for( l = k; l <= indx; l++ )
                        {
                            tmp1 = cctable[l];
                            cctable[l] = tmp2;
                            tmp2 = tmp1;
                        }
                        break;
                    }
                }

                float sum1=0, sum2=0;
                //check "once-off" changes
                for(k = 0; PV_CC(k) && k < model->params.N1cc; k++ )
                {
                    sum1 += PV_CC(k);
                    sum2 += PVB_CC(k);
                }
                if( sum1 > model->params.T ) stat->is_trained_dyn_model = 1;
                
                diff = sum1 - stat->Pbcc * sum2;
                //update stat table
                if( diff >  model->params.T )
                {
                    //printf("once off change at motion mode\n");
                    //new BG features are discovered
                    for( k = 0; PV_CC(k) && k < model->params.N1cc; k++ )
                    {
                        PVB_CC(k) =
                            (PV_CC(k)-stat->Pbcc*PVB_CC(k))/(1-stat->Pbcc);
                    }
                    assert(stat->Pbcc<=1 && stat->Pbcc>=0);
                }
            }

            //case of stational pixel
            if( !((unsigned char*)model->Ftd)[i*mask_step+j] )
            {
                float alpha = stat->is_trained_st_model ? model->params.alpha2 : model->params.alpha3;
                float diff = 0;
                int dist, min_dist = 2147483647, indx = -1;

                //update Pb
                stat->Pbc *= (1.f-alpha);
                if( !((unsigned char*)model->foreground)[i*mask_step+j] )
                {
                    stat->Pbc += alpha;
                }

                //find best Vi match
                for( k = 0; k < model->params.N2c; k++ )
                {
                    // Exponential decay of memory
                    PV_C(k) *= (1-alpha);
                    PVB_C(k) *= (1-alpha);
                    if( PV_C(k) < MIN_PV )
                    {
                        PV_C(k) = 0;
                        PVB_C(k) = 0;
                        continue;
                    }
                    
                    dist = 0;
                    for( l = 0; l < nChannels; l++ )
                    {
                        int val = abs( V_C(k,l) - curr_data[l] );
                        if( val > deltaC ) break;
                        dist += val;
                    }
                    if( l == nChannels && dist < min_dist )
                    {
                        min_dist = dist;
                        indx = k;
                    }
                }

                if( indx < 0 )
                {//N2th elem in the table is replaced by a new features
                    indx = model->params.N2c - 1;
                    PV_C(indx) = alpha;
                    PVB_C(indx) = alpha;
                    //udate Vt
                    for( l = 0; l < nChannels; l++ )
                    {
                        V_C(indx,l) = curr_data[l];
                    }
                } else
                {//update
                    PV_C(indx) += alpha;
                    if( !((unsigned char*)model->foreground)[i*mask_step+j] )
                    {
                        PVB_C(indx) += alpha;
                    }
                }

                //re-sort Ct table by Pv
                for( k = 0; k < indx; k++ )
                {
                    if( PV_C(k) <= PV_C(indx) )
                    {
                        //shift elements
                        BGPixelCStatTable tmp1, tmp2 = ctable[indx];
                        for( l = k; l <= indx; l++ )
                        {
                            tmp1 = ctable[l];
                            ctable[l] = tmp2;
                            tmp2 = tmp1;
                        }
                        break;
                    }
                }

                //check "once-off" changes
                float sum1=0, sum2=0;
                for( k = 0; PV_C(k) && k < model->params.N1c; k++ )
                {
                    sum1 += PV_C(k);
                    sum2 += PVB_C(k);
                }
                diff = sum1 - stat->Pbc * sum2;
                if( sum1 > model->params.T ) stat->is_trained_st_model = 1;

                //update stat table
                if( diff >  model->params.T )
                {
                    //printf("once off change at stat mode\n");
                    //new BG features are discovered
                    for( k = 0; PV_C(k) && k < model->params.N1c; k++ )
                    {
                        PVB_C(k) = (PV_C(k)-stat->Pbc*PVB_C(k))/(1-stat->Pbc);
                    }
                    stat->Pbc = 1 - stat->Pbc;
                }
            }//if !(change detection) at pixel (i,j)

            //update the reference BG image
			//2. 更新参考背景图像
            if( !((unsigned char*)model->foreground)[i*mask_step+j])//对背景点进行更新
            {
                unsigned char* ptr = ((unsigned char*)model->background) + i*WidthStep+j*nChannels;
                
                if( !((unsigned char*)model->Ftd)[i*mask_step+j] &&
                    !((unsigned char*)model->Fbd)[i*mask_step+j] )//对缓慢逐渐变化的背景点进行更新
                {
                    //apply IIR filter
                    for( l = 0; l < nChannels; l++ )
                    {
                        int a = int(ptr[l]*(1 - model->params.alpha1) + model->params.alpha1*curr_data[l]);
                        ptr[l] = (unsigned char)a;
                    }
                }
                else                                              //对运动的背景点（如摇曳的树叶背景、闪烁的水纹）进行更新
                {
                    for( l = 0; l < nChannels; l++ )
                    {
                        ptr[l] = curr_data[l];
                    }
                }
            }
        }//j
    }//i

    //keep prev frame
	memcpy(model->prev_frame,Image,sizeof(unsigned char)*(Width*Height*nChannels));

}
/*************************************************************************
**	函数名称：
**		ReleaseFGDStatModel()
**	参数
**		FGDStatModel* _model  贝叶斯状态模型
**	返回值：
**		无 
**	说明:
**		检测结束释放贝叶斯模型相关资源
************************************************************************/
void ReleaseFGDStatModel( FGDStatModel* _model )
{
    if( _model )
    {
        FGDStatModel* model = _model;
        if( model->pixel_stat )
        {
			for (int i=0;i<model->width*model->height;i++)
			{
				free( model->pixel_stat[i].ctable );
				free( model->pixel_stat[i].cctable );
			}
            free( model->pixel_stat );
        }
		
        free( model->Ftd );
        free( model->Fbd );
        free( model->foreground );
        free( model->background );
        free( model->prev_frame );
        free( _model );
    }
}
/************************************************************************/
/*                    C.基于简单背景建模的运动目标检测                  */
/************************************************************************/

/*************************************************************************
**	函数名称：
**		SetBgHistgram()
**	参数
**		short int *pHistgram  保存统计的背景直方图
**      BYTE *sGray           当前进入的图像
**		int nWidth            图像的宽度
**		int nHeight           图像的高度
**	返回值：
**		无 
**	说明:
**		建立背景时域统计直方图
************************************************************************/
void SetBgHistgram(short int *pHistgram, BYTE *sGray,int nWidth, int nHeight)
{
	for(int j=0;j<nHeight;j++)
	{
		for(int i=0;i<nWidth;i++)
		{
			int p=*(sGray+j*nWidth+i)+256*(j*nWidth+i);
			*(pHistgram+p)=*(pHistgram+p)+1;
		}
	}	
	
}
/*************************************************************************
**	函数名称：
**		GetPeakBgByHistgram()
**	参数
**		short int *pHistgram  统计的背景直方图
**      BYTE *sBg             估计的背景图像
**		int nWidth            图像的宽度
**		int nHeight           图像的高度
**	返回值：
**		无 
**	说明:
**		通过时域直方图建立峰值背景
************************************************************************/
void GetPeakBgByHistgram(short int *pHistgram,BYTE *sBg,int nWidth,int nHeight)
{
	int j,i,k;
	for(j=0;j<nHeight;j++)
	{
		for(i=0;i<nWidth;i++)
		{
			int c=0; 
			int color=0;
			for(k=0;k<256;k++)
			{
				if(*(pHistgram+256*(j*nWidth+i)+k)>c)
				{
					c=*(pHistgram+256*(j*nWidth+i)+k);					
					color=k;
				}
			}
			
			*(sBg+j*nWidth+i)=color;
			
		}
	}		
}
/*************************************************************************
**	函数名称：
**		GetMedianBgByHistgram()
**	参数
**		short int *pHistgram  统计的背景直方图
**      BYTE *sBg             估计的背景图像
**		int nWidth            图像的宽度
**		int nHeight           图像的高度
**      int countFrames       已经处理的帧数
**	返回值：
**		无 
**	说明:
**		通过时域直方图建立中值背景
************************************************************************/
void GetMedianBgByHistgram(short int *pHistgram,BYTE *sBg,int nWidth,int nHeight,int countFrames)
{
	int j,i,k;
	int MedianValueIndex = (countFrames+1) /2;
	for(j=0;j<nHeight;j++)
	{
		for(i=0;i<nWidth;i++)
		{
			int count=0; 
			int MedianValue=0;
			for(k=0;k<256;k++)
			{
				count +=*(pHistgram+256*(j*nWidth+i)+k);
				if(count>MedianValueIndex)
				{					
					MedianValue=k;
					break;
				}
			}
			
			*(sBg+j*nWidth+i)=MedianValue;
			
		}
	}		
}
/*************************************************************************
**	函数名称：
**		GetRuningAvgBg()
**	参数
**		short int *pHistgram  统计的背景直方图
**      BYTE *sBg             估计的背景图像
**		int nWidth            图像的宽度
**		int nHeight           图像的高度
**      int countFrames       已经处理的帧数
**	返回值：
**		无 
**	说明:
**		通过时域直方图建立中值背景
************************************************************************/
void GetRuningAvgBg(unsigned char* Bg,unsigned char*Image,int Width,int Height,int countFrames)
{
	for(int j=0;j<Height;j++)
	{
		for(int i=0;i<Width;i++)
		{
			double p = double (*(Bg+j*Width+i))*countFrames+double (*(Image+j*Width+i));
			*(Bg+j*Width+i) = (unsigned char)(p/(countFrames+1));
		}
	}
}
/*************************************************************************
**	函数名称：
**		CreateSimpleBGModel()
**	参数
**      unsigned char*Image   初始帧图像
**		int Width             图像的宽度
**		int Height            图像的高度
**      int Thresh            检测二值化阈值
**		int SimBGType         简单背景方式1: 均值背景 2：中值背景 3：直方图峰值背景
**	返回值：
**		简单背景模型 
**	说明:
**		初始时创建简单背景统计模型相关参数
************************************************************************/
SimpleBGModel* CreateSimpleBGModel( unsigned char*Image,int Width,int Height,int Thresh,int SimBGType)
{
	//给SimpleBGModel分配内存
	SimpleBGModel* p_model = 0;
    p_model = (SimpleBGModel*)malloc( sizeof(*p_model) );
    memset( p_model, 0, sizeof(*p_model) );

	//填充SimpleBGModel
	p_model->SimBGType = SimBGType;
    p_model->thresh = Thresh;
	p_model->width = Width;
	p_model->height = Height;
	p_model->countFrames = 1;
	p_model->foreground = (unsigned char*)malloc(sizeof(unsigned char)*(Width*Height*1));
	p_model->background = (unsigned char*)malloc(sizeof(unsigned char)*(Width*Height*1));
	if (p_model->SimBGType==SIMPLEBG_HISTPEAK || p_model->SimBGType==SIMPLEBG_MEDIAN)
	{
		p_model->pBackHist = (short int *)malloc(sizeof(short int)*(Width*Height*256));
		memset(p_model->pBackHist,0,sizeof(short int)*Width*Height*256);
		SetBgHistgram(p_model->pBackHist,Image,Width,Height);
	}
	memcpy(p_model->background,Image,sizeof(unsigned char)*(Width*Height*1));

	return p_model;
}  
/*************************************************************************
**	函数名称：
**		UpdateSimpleModel()
**	参数
**      unsigned char*Image    当前要处理的图像
**		int Width              图像的宽度
**		int Height             图像的高度
**		SimpleBGModel* _model  简单背景模型
**	返回值：
**		无 
**	说明:
**		通过新输入图像即时的更新当前简单背景模型
************************************************************************/
void UpdateSimpleModel( unsigned char*Image,int Width,int Height,SimpleBGModel* model )
{
	if (model->SimBGType==SIMPLEBG_AVG)                        //背景时域均值建立背景
	{
		GetRuningAvgBg(model->background,Image,Width,Height,model->countFrames);
		InfraredSubtract(model->foreground,Image,model->background,Width,Height,model->thresh);
		++model->countFrames;
	}
	else if (model->SimBGType==SIMPLEBG_MEDIAN)                    //背景时域中值建立背景
	{
		++model->countFrames;
		SetBgHistgram(model->pBackHist,Image,Width,Height);
		GetMedianBgByHistgram(model->pBackHist,model->background,Width,Height,model->countFrames);
		InfraredSubtract(model->foreground,Image,model->background,Width,Height,model->thresh);
	}
	else                                             //背景直方图峰值建立背景
	{
		SetBgHistgram(model->pBackHist,Image,Width,Height);
		GetPeakBgByHistgram(model->pBackHist,model->background,Width,Height);
		InfraredSubtract(model->foreground,Image,model->background,Width,Height,model->thresh);
	}
}
/*************************************************************************
**	函数名称：
**		ReleaseSimpleModel()
**	参数
**		SimpleBGModel* _model  简单背景模型
**	返回值：
**		无 
**	说明:
**		检测结束释放简单背景模型相关资源
************************************************************************/
void ReleaseSimpleModel(SimpleBGModel* _model )
{
	if( _model )
    {
        SimpleBGModel* model = _model;
        free( model->foreground );
        free( model->background );
		if (model->SimBGType==SIMPLEBG_HISTPEAK)
			free(model->pBackHist);
        free( _model );
    }
}
/************************************************************************/
/*                                 结束                                 */
/************************************************************************/ 
/************************************************************************/
/*                    D.基于核密度的背景建模的运动目标检测          */
/************************************************************************/ 




/*************************************************************************
**	函数名称：
**		excle()
**	参数
**       i   值为i的像素 
**	
**	返回值：
**		无 
**	说明:
**		建立查找表加快计算速度
************************************************************************/

void excle(double* exclee) 
{
	if(!exclee)
		return;
	
	int k;
	double cova=2/0.68*sqrt(2);	
	double pi=3.14159265;
	
	for(k=0;k<256;k++)
	{
		exclee[k]=exp((-k*k/(2*cova*cova))/(cova*sqrt(2*pi)));
	}
}
/*************************************************************************
**	函数名称：
**		createkdemodel()
**	参数

SampleImage   前10帧为样本图像最后一帧为当前处理的图像，共11帧
height        每一帧图像的行数
lineByte      每帧帧图像的列数
N_L           样本所含的帧数+当前帧共N_L帧
exclee        索引表
fore          最后的目标提取结果
thre          手动设置阈值
	
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
					double thre)					   
{  
	int i,j,k,index;
	int SampleNum=N_L-1; // N_L取11，SampleNum取10


    unsigned char* foreground=new unsigned char[height*lineByte];
	memset( foreground, 0, sizeof(unsigned char)*height*lineByte);
	
    memcpy(foreground, SampleImage[SampleNum], sizeof(unsigned char)*height*lineByte);//将当前帧复制到数组foreground中，备份


	double *ab=new double [lineByte*height]; 
 	memset( ab, 0, sizeof(double)*height*lineByte);


	for (i=0;i<height;i++)
    {
        for (j=0;j<lineByte;j++)
        {
            k = i*lineByte + j;
            for (index=0;index<SampleNum;index++)
            {
                ab[k] += exclee[abs(SampleImage[SampleNum][k]-SampleImage[index][k])]; 
            }
        }
    }
	

	//*****************阈值**********************************
	
	
	for (j=0;j<lineByte*height;j++)
	{
		if (ab[j]<thre)
		{
			foreground[j]=255;
		}
		
	}
 
	//二值化
	for (j=0;j<lineByte*height;j++)
	{
		if (foreground[j]<255)
		{
			fore[j]=0;
		}
		else
			fore[j] = 255;
	}


	


 	//调用去除噪声点的函数对fore去除噪声点

	delete []foreground;
	foreground = NULL;

 
	delete []ab;
	ab = NULL;

 }






/************************************************************************/
/*                                 结束                                 */
/************************************************************************/ 

