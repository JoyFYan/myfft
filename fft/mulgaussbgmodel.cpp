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
/*                    �ġ�����Ŀ����ģ��                              */
/*                                                                      */
/************************************************************************/

/************************************************************************/
/*                    4.2 ��չĿ������غ���(��)                      */
/************************************************************************/ 

/************************************************************************/
/*                    A.���ڶ��˹������ģ���˶�Ŀ����                */
/************************************************************************/ 

/*************************************************************************
**	�������ƣ�
**		CreateGaussianBGModel()
**	����
**      unsigned char*Image   ��ʼ֡ͼ��
**		int Width             ͼ��Ŀ��
**		int Height            ͼ��ĸ߶�
**      int nChannels         ͼ���ͨ����
**		GaussBGStatModelParams* parameters  ��˹����ģ�Ͳ���
**	����ֵ��
**		��˹����ͳ��ģ�� 
**	˵��:
**		��ʼʱ������˹����ͳ��ģ����ز���
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
        params.win_size = BGFG_MOG_WINDOW_SIZE;                        // ��ʼ���׶ε�֡�����û��Զ���ģ��ѧϰ��a=1/win_size;
        params.bg_threshold = BGFG_MOG_BACKGROUND_THRESHOLD;
        params.std_threshold = BGFG_MOG_STD_THRESHOLD;
        params.weight_init = BGFG_MOG_WEIGHT_INIT;
        params.variance_init = BGFG_MOG_SIGMA_INIT*BGFG_MOG_SIGMA_INIT;//����
        params.n_gauss = BGFG_MOG_NGAUSSIANS;                          //��˹�ֲ������ĸ���
    }
    else
    {
        params = *parameters;                                           //�û��Զ������
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
    var_init = 2 * params.std_threshold * params.std_threshold;            //��ʼ������
   
    for( i = 0, p = 0, n = 0; i < Height; i++ )        //��
    {
        for( j = 0; j < Width; j++, n++ )               //��
        {
			//���¼����ǶԵ�һ����˹��������ʼ��
            bg_model->g_point[n].g_values[0].weight = 1;     //Ȩֵ��Ϊ1
            bg_model->g_point[n].g_values[0].match_sum = 1;  //��˹������ƥ��Ĵ���
            for( m = 0; m < nChannels; m++)
            {
                bg_model->g_point[n].g_values[0].variance[m] = var_init;
				//��ֵ��Ϊ��ǰ���ص�ֵ
                bg_model->g_point[n].g_values[0].mean[m] = (unsigned char)Image[p + m];
            }

			//����һ����ĸ�˹�ֲ������ĳ�ʼ��(��ֵ��Ȩֵ��ƥ�����������)
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
    }//g_point[]:���أ�g_values[]:��˹�ֲ�������mean[]:ͨ��

    bg_model->countFrames = 0; //���ó�ʼ֡Ϊ0
	return bg_model;
}

/*************************************************************************
**	�������ƣ�
**		UpdateGaussianBGModel()
**	����
**      unsigned char*Image      ��ǰҪ�����ͼ��
**		int Width                ͼ��Ŀ��
**		int Height               ͼ��ĸ߶�
**      int nChannels            ͼ���ͨ����
**		GaussBGModel*  bg_model  ��˹����ģ��
**	����ֵ��
**		��
**	˵��:
**		ͨ��������ͼ��ʱ�ĸ��µ�ǰ��˹����ģ�ͣ�����ģ�͵ĸ��£�
**      ����Ҫ���¸�˹�ֲ������Ĳ�������Ҫ���¸���˹������Ȩ��
************************************************************************/
void UpdateGaussianBGModel(  unsigned char*Image,int Width,int Height,int nChannels, GaussBGModel*  bg_model )
{
    int i, j, k;
    bg_model->countFrames++;
    
    for( i = 0; i < Height; i++ )
    {
        for( j = 0; j < Width; j++ )
        {
            int match[BGFG_MOG_MAX_NGAUSSIANS];                //�Ը�˹��������ǣ�match[m]=1��ʾ����mΪƥ��ĸ�˹�ֲ�����
            double sort_key[BGFG_MOG_MAX_NGAUSSIANS];          //���������ÿ����˹�����ľ�ֵ�뷽���ֵ
            const int n = i*Width+j;
            const int p = n*nChannels;
            
            // A few short cuts
            GaussBGPoint* g_point = &bg_model->g_point[n];
            const GaussBGStatModelParams bg_model_params = bg_model->params;
            double pixel[4];                                 // pixel[]������ǰ���صĸ�ͨ��RGBֵ
            int no_match;
            
            for( k = 0; k < nChannels; k++ )
                pixel[k] = (unsigned char)Image[p+k];

            //����Ƿ����뵱ǰ����ƥ��ĸ�˹����
            no_match = MatchTest( pixel, nChannels, match, g_point, &bg_model_params );
            if( bg_model->countFrames == bg_model->params.win_size )   //����ѧϰ�����Զ���                 
            {
                UpdateFullWindow( pixel, nChannels, match, g_point, &bg_model->params );//����ƥ�亯��ʱ����ظ�˹�����ĸ���
                if( no_match == -1)        //û�и�˹����ƥ��ʱ�������µĸ�˹����
                    UpdateFullNoMatch( Image, nChannels, p, match, g_point, &bg_model_params );
            }
            else                                                       //����ѧϰ����ȡ���ڸ�˹�ֲ���ƥ�����
            {
                UpdatePartialWindow( pixel, nChannels, match, g_point, &bg_model_params );//����ƥ�亯��ʱ����ظ�˹�����ĸ���
                if( no_match == -1)       //û�и�˹����ƥ��ʱ�������µĸ�˹����
                    UpdatePartialNoMatch( pixel, nChannels, match, g_point, &bg_model_params );
            }
			//�����˹�ֲ��ľ�ֵ�ͷ����ֵ
            GetSortKey( nChannels, sort_key, g_point, &bg_model_params );
			//�Ը�˹�ֲ��ľ�ֵ�ͷ����ֵ���н�������
            InsertionSortGaussians( g_point, sort_key, (GaussBGStatModelParams *)&bg_model_params );
			//������⣬ǰ���ĸ���
            BackgroundTest( nChannels, n, p, match, bg_model );
        }
    }
}
/*************************************************************************
**	�������ƣ�
**		ReleaseGaussianBGModel()
**	����
**		GaussBGModel** _bg_model  ��˹����ģ��
**	����ֵ��
**		�� 
**	˵��:
**		�������ͷŸ�˹����ģ�������Դ
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
**	�������ƣ�
**		MatchTest()
**	����
**      double* src_pixel                              ����ǰ֡�и����ص�RGBֵ
**		int nChannels                                  ͼ���ͨ����
**      int* match                                     ��ǰ�����ͼ������ֵ��ƥ��ֵָ��
**		GaussBGPoint* g_point                          ��˹������ָ��
**		const GaussBGStatModelParams *bg_model_params  ��˹����ģ�Ͳ���
**	����ֵ��
**		Ϊ-1��û��ƥ���λ�ã����򷵻ظ�˹ƥ�亯����λ��
**	˵��:
**		�õ�ǰ���ص�ֵ���Ѵ��ڵĸ�˹�ֲ������Ƚϣ������Ƿ����ƥ���
**      �ĸ�˹�ֲ�������������򷵻� kֵ(��˹�ֲ����������)
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
            double d = g_point->g_values[k].mean[m]- src_pixel[m];//ͨ��m��ԭʼģ��ֵ�뵱ǰ���ص�ֵ֮��
            sum_d2 += (d*d);
            var_threshold += g_point->g_values[k].variance[m];
        }  //difference < STD_LIMIT*STD_LIMIT or difference**2 < STD_LIMIT*STD_LIMIT*VAR

		//��ǰsum_d2Ϊd0,d1,d2��ƽ���ͣ�var_threshold��ֵΪ���ظ�ͨ������֮��

        var_threshold = bg_model_params->std_threshold*bg_model_params->std_threshold*var_threshold;

        if(sum_d2 < var_threshold)    //�鿴�Ƿ������ĳ��˹�ֲ�ƥ��
		{
            match[k] = 1;
            matchPosition = k;
            break;                   //����͵�k����˹����ƥ�䣬����ֹ�����������ƥ��
        }
    }
    
    return matchPosition;
}
/*************************************************************************
**	�������ƣ�
**		UpdateFullWindow()
**	����
**      double* src_pixel                              ����ǰ֡�и����ص�RGBֵ
**		int nChannels                                  ͼ���ͨ����
**      int* match                                     ��ǰ�����ͼ������ֵ��ƥ��ֵָ��
**		GaussBGPoint* g_point                          ��˹������ָ��
**		const GaussBGStatModelParams *bg_model_params  ��˹����ģ�Ͳ���
**	����ֵ��
**		��
**	˵��:
**		����ÿ����˹�ֲ���Ȩֵ(��ƥ��ĸ�˹����k�Ӵ�Ȩֵ����������СȨֵ)�����ǰ��Ľ����
**      ����ƥ��ĸ�˹�ֲ�����k,����Ҫ�ٶԵ�k����˹�ֲ������ľ�ֵmean�ͷ���variance������
************************************************************************/
void UpdateFullWindow( double* src_pixel, int nChannels, int* match,
					  GaussBGPoint* g_point,
					  const GaussBGStatModelParams *bg_model_params )
{
	//�û��Զ���ģ��ѧϰ��a
    const double learning_rate_weight = (1.0/(double)bg_model_params->win_size);

    for(int k = 0; k < bg_model_params->n_gauss; k++)
	{
		//��ÿ����˹�ֲ���Ȩֵ��������w=(1-a)w+a*m (a:ģ��ѧϰ�ʣ�m��ƥ�䣬ƥ�����1����ƥ�����0)
        g_point->g_values[k].weight = g_point->g_values[k].weight +
            (learning_rate_weight*((double)match[k] -g_point->g_values[k].weight));

		//�������ƥ��ĸ�˹�ֲ�����k����ǰ����Ϊ�������أ�,����Ҫ�ٶԵ�k����˹�ֲ������ľ�ֵmean�ͷ���variance����
        if(match[k])
		{
            double learning_rate_gaussian = (double)match[k]/(g_point->g_values[k].weight*
                (double)bg_model_params->win_size);                 //����ѧϰ��p(p=a/w)
            for(int m = 0; m < nChannels; m++)
			{
				//�������¹�ʽ��u=(1-p)*u0+p*x; o*o=(1-p)*o*o+p*tmpDiff*tmpDiff
				//��ǰ���ص�ͨ��m��ֵ��ԭʼģ��ֵ֮��
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
**	�������ƣ�
**		UpdatePartialWindow()
**	����
**      double* src_pixel                              ����ǰ֡�и����ص�RGBֵ
**		int nChannels                                  ͼ���ͨ����
**      int* match                                     ��ǰ�����ͼ������ֵ��ƥ��ֵָ��
**		GaussBGPoint* g_point                          ��˹������ָ��
**		const GaussBGStatModelParams *bg_model_params  ��˹����ģ�Ͳ���
**	����ֵ��
**		��
**	˵��:
**		�����еĸ�˹�ֲ�����������.����ÿ����˹�ֲ���Ȩֵ�������������ǰ��Ľ���д���ƥ��ĸ�˹�ֲ�����k,
**      ����Ҫ�ٶԵ�k����˹�ֲ�������match_sum�޸ģ����ն���Щƥ��ĸ�˹�ֲ�����k�Ĳ���match_sum>0������ֵmean�ͷ���variance����
************************************************************************/
void UpdatePartialWindow( double* src_pixel, int nChannels, int* match, GaussBGPoint* g_point, const GaussBGStatModelParams *bg_model_params )
{
    int k, m;
    int window_current = 0;
    
    for( k = 0; k < bg_model_params->n_gauss; k++ )
        window_current += g_point->g_values[k].match_sum;//window_currentΪk����˹�ֲ�������match_sumֵ֮��
    
    for( k = 0; k < bg_model_params->n_gauss; k++ )
    {
        g_point->g_values[k].match_sum += match[k];      //����ƥ��ĸ�˹�ֲ�����k��match_sumֵ
        double learning_rate_weight = (1.0/((double)window_current + 1.0)); //increased by one since sum

		//����ÿ����˹�ֲ���Ȩֵ
        g_point->g_values[k].weight = g_point->g_values[k].weight +
            (learning_rate_weight*((double)match[k] - g_point->g_values[k].weight));

        //�������ƥ��ĸ�˹�ֲ�����k����ǰ����Ϊ�������أ�,����Ҫ�ٶԵ�k����˹�ֲ������ľ�ֵmean�ͷ���variance����
        if( g_point->g_values[k].match_sum > 0 && match[k] )
        {
            double learning_rate_gaussian = (double)match[k]/((double)g_point->g_values[k].match_sum);//����ѧϰ��p(p=a/w)
            for( m = 0; m < nChannels; m++ )
            {
				//�������¹�ʽ��u=(1-p)*u0+p*x; o*o=(1-p)*o*o+p*tmpDiff*tmpDiff
				//��ǰ���ص�ͨ��m��ֵ��ԭʼģ��ֵ֮��
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
**	�������ƣ�
**		UpdateFullNoMatch()
**	����
**      unsigned char* gm_image                        ��ǰҪ�����ͼ��
**		int nChannels                                  ͼ���ͨ����
**      int p                                          ͼ�����ص��λ��
**      int* match                                     ��ǰ�����ͼ������ֵ��ƥ��ֵָ��
**		GaussBGPoint* g_point                          ��˹������ָ��
**		const GaussBGStatModelParams *bg_model_params  ��˹����ģ�Ͳ���
**	����ֵ��
**		��
**	˵��:
**		�����еĸ�˹��������ƥ��ʱ��˵�����µķֲ����֣���Ҫ��ԭ��˹������sort_key��С��
**      �滻Ϊ�µĸ�˹����(ȨֵС�������),����ĸ�˹������Ӧ��ֻ�����Ȩֵ
************************************************************************/
void UpdateFullNoMatch( unsigned char* gm_image, int nChannels, int p, int* match,
								 GaussBGPoint* g_point,
								 const GaussBGStatModelParams *bg_model_params)
{
    int k, m;
    double alpha;
    int match_sum_total = 0;
	
    //���һ����˹�ֲ�����Ϊ�µĸ�˹����
    g_point->g_values[bg_model_params->n_gauss - 1].match_sum = 1;  //���µĸ�˹�ֲ�������match_sum��Ϊ1

    
    //���������˹�ֲ�match_sum ֮��
    for( k = 0; k < bg_model_params->n_gauss ; k++ )
        match_sum_total += g_point->g_values[k].match_sum;
    //���µĸ�˹�ֲ�������һ����С��Ȩֵ
    g_point->g_values[bg_model_params->n_gauss - 1].weight = 1./(double)match_sum_total;

	//���µĸ�˹�ֲ�������variance[m]ȫ����Ϊvariance_init��mean[m]��ֵ��Ϊ��ǰ���ظ�ͨ����ֵ
    for( m = 0; m < nChannels ; m++ )
    {
        g_point->g_values[bg_model_params->n_gauss - 1].variance[m] = bg_model_params->variance_init;
        g_point->g_values[bg_model_params->n_gauss - 1].mean[m] = (unsigned char)gm_image[p + m];
    }

    //�������ĸ�˹�ֲ�������Ȩֵ���£�w=(1-a)*w+a*m (a:ģ��ѧϰ�ʣ�m��ƥ�䣬ƥ�����1����ƥ�����0)
    alpha = 1.0 - (1.0/bg_model_params->win_size);
    for( k = 0; k < bg_model_params->n_gauss - 1; k++ )
    {
        g_point->g_values[k].weight *= alpha;   
        if( match[k] )                          //���ƥ��Ļ� w=(1-a)*w+a*m ���Ϻ���
            g_point->g_values[k].weight += alpha;
    }
}
/*************************************************************************
**	�������ƣ�
**		UpdatePartialNoMatch()
**	����
**      double *pixel                                  ��ǰҪ�����ͼ������ֵ
**		int nChannels                                  ͼ���ͨ����
**      int* match                                     ��ǰ�����ͼ������ֵ��ƥ��ֵָ��
**		GaussBGPoint* g_point                          ��˹������ָ��
**		const GaussBGStatModelParams *bg_model_params  ��˹����ģ�Ͳ���
**	����ֵ��
**		��
**	˵��:
**		�����еĸ�˹��������ƥ��ʱ����Ҫ��ԭ��˹������sort_key��С��
**      �滻Ϊ�µĸ�˹����(ȨֵС�������),Ȼ���һ�����еĸ�˹�ֲ�Ȩֵ
************************************************************************/
void UpdatePartialNoMatch(double *pixel,
                        int nChannels,
                        int* /*match*/,
                        GaussBGPoint* g_point,
                        const GaussBGStatModelParams *bg_model_params)
{
    int k, m;
    //������˹�ֲ�������match_sum��Ϊ1
    g_point->g_values[bg_model_params->n_gauss - 1].match_sum = 1;
    
    //���������˹�ֲ�match_sum ֮��
    int match_sum_total = 0;
    for(k = 0; k < bg_model_params->n_gauss ; k++)
        match_sum_total += g_point->g_values[k].match_sum;

	//�����һ����˹�ֲ�������variance[m]ȫ����Ϊvariance_init��mean[m]��ֵ��Ϊ��ǰ���ظ�ͨ����ֵ
    for(m = 0; m < nChannels; m++)
    {
        g_point->g_values[bg_model_params->n_gauss - 1].variance[m] = bg_model_params->variance_init;
        g_point->g_values[bg_model_params->n_gauss - 1].mean[m] = pixel[m];
    }
	//��һ�����еĸ�˹�ֲ�Ȩֵ
    for(k = 0; k < bg_model_params->n_gauss; k++)
    {
        g_point->g_values[k].weight = (double)g_point->g_values[k].match_sum /
            (double)match_sum_total;
    }
}
/*************************************************************************
**	�������ƣ�
**		GetSortKey()
**	����
**      const int nChannels                            ͼ���ͨ����
**		double* sort_key                               weight/sqrt(variance_sum)��ֵ��ָ��
**		const GaussBGPoint* g_point                    ��˹������ָ��
**		const GaussBGStatModelParams *bg_model_params  ��˹����ģ�Ͳ���
**	����ֵ��
**		��
**	˵��:
**		���������˹�ֲ�weight/sqrt(variance_sum)��ֵ��
**      ���潫�Ը�ֵ�������򣨸�ֵԽ�����ʾ�����Ŀ����Ծ�Խ��
************************************************************************/
void GetSortKey( const int nChannels, double* sort_key, const GaussBGPoint* g_point,
						  const GaussBGStatModelParams *bg_model_params )
{
    int k, m;
    for( k = 0; k < bg_model_params->n_gauss; k++ )
    {
        // ���ⱻ���
        if( g_point->g_values[k].match_sum > 0 )
        {
            // �����˹�ɷָ��������Ƕ�����
            double variance_sum = 0.0;
            for( m = 0; m < nChannels; m++ )                              //�����ͨ�������
                variance_sum += g_point->g_values[k].variance[m];
            
            sort_key[k] = g_point->g_values[k].weight/sqrt(variance_sum); //�����˹�ֲ�weight/sqrt(variance_sum)
        }
        else
            sort_key[k]= 0.0;
    }
}
/*************************************************************************
**	�������ƣ�
**		InsertionSortGaussians()
**	����
**		const GaussBGPoint* g_point                    ��˹������ָ��
**		double* sort_key                               weight/sqrt(variance_sum)��ֵ��ָ��
**		GaussBGStatModelParams *bg_model_params  ��˹����ģ�Ͳ���
**	����ֵ��
**		��
**	˵��:
**		�Ը�����˹�ֲ�weight/sqrt(variance_sum)��ֵ���в������򣬽�������
************************************************************************/
void InsertionSortGaussians( GaussBGPoint* g_point, double* sort_key, GaussBGStatModelParams *bg_model_params )
{
    int i, j;
    for( i = 1; i < bg_model_params->n_gauss; i++ )
    {
        double index = sort_key[i];
        for( j = i; j > 0 && sort_key[j-1] < index; j-- )           //��������
        {
            double temp_sort_key = sort_key[j];                     //��������sortkey��ֵ
            sort_key[j] = sort_key[j-1];
            sort_key[j-1] = temp_sort_key;
            
            GaussBGValues temp_gauss_values = g_point->g_values[j]; //����sortkey��Ӧ�ĸ�˹�ֲ���ֵ
            g_point->g_values[j] = g_point->g_values[j-1];
            g_point->g_values[j-1] = temp_gauss_values;
        }
    }
}
/*************************************************************************
**	�������ƣ�
**		BackgroundTest()
**	����
**      const int nChannels     ͼ���ͨ����
**		int n                   ��˹������λ�� i��width+j
**		int p                   ͼ�����ص�λ�� n��nChannels
**      int *match              ͼ������ֵ��ƥ��ֵָ��
**		GaussBGModel* bg_model  ��˹����ģ��
**	����ֵ��
**		��
**	˵��:
**		��⵱ǰ�����ͼ������ֵΪǰ��ֵ���Ǳ���ֵ
************************************************************************/
void BackgroundTest( const int nChannels, int n, int p, int *match, GaussBGModel* bg_model )
{
    int m, b;
    unsigned char pixelValue = (unsigned char)255; // ��ʼֵ��Ϊ255���ҵ���˹ƥ�亯�����ٽ�������Ϊ0
    double weight_sum = 0.0;
    GaussBGPoint* g_point = bg_model->g_point;
    
    for( m = 0; m < nChannels; m++)                 //�õ���˹ģ�͵Ĺ��Ʊ���
        bg_model->background[p+m]   = (unsigned char)(g_point[n].g_values[0].mean[m]+0.5);
    
    for( b = 0; b < bg_model->params.n_gauss; b++)  //������˹��������ƥ�䣬���������Ƿ�ƥ��
    {
        weight_sum += g_point[n].g_values[b].weight;
        if( match[b] )                              //���Ϊ�棬˵������������ĳ��˹����ƥ�䣬������Ϊ����
            pixelValue = 0;
        if( weight_sum > bg_model->params.bg_threshold )
            break;                                  //���if���Ϊ�棬��ǰb����˹�ֲ���ѡΪ���������ĺ���
    }
    
    bg_model->foreground[p/nChannels] = pixelValue; //����ǰ��ֵ
}
/************************************************************************/
/*                    B.���ڱ�Ҷ˹ģ�͵ı�����ģ���˶�Ŀ����          */
/************************************************************************/ 

/*************************************************************************
**	�������ƣ�
**		CreateFGDStatModel()
**	����
**      unsigned char*Image   ��ʼ֡ͼ��
**		int Width             ͼ��Ŀ��
**		int Height            ͼ��ĸ߶�
**      int nChannels         ͼ���ͨ����
**		FGDStatModelParams* parameters  ��Ҷ˹����ģ�Ͳ���
**	����ֵ��
**		��Ҷ˹����ͳ��ģ�� 
**	˵��:
**		��ʼʱ������Ҷ˹����ͳ��ģ����ز���
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
**	�������ƣ�
**		_max_element()
**	����
**      double* start          ���е�ͷָ��
**		double* end            ���е�βָ��
**	����ֵ��
**		���ֵָ��
**	˵��:
**		�ҳ����е����ֵָ��
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
**	�������ƣ�
**		ChangeDetection()
**	����
**      unsigned char * prev_Image            ��ǰ֡ͼ��
**		unsigned char * curr_Image            ��ǰҪ�����ͼ��
**      unsigned char * changeMask_Image      ǰ����ֵͼ
**      int Width                             ͼ��Ŀ��
**		int Height                            ͼ��ĸ߶�
**      int nChannels                         ͼ���ͨ����
**	����ֵ��
**		�� 
**	˵��:
**		ͨ���򵥵�֡��ͼ�����Ӧ�ָ�õ�ǰ����ֵͼ
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
**	�������ƣ�
**		UpdateFGDStatModel()
**	����
**      unsigned char*Image   ��ǰҪ�����ͼ��
**		int Width             ͼ��Ŀ��
**		int Height            ͼ��ĸ߶�
**      int nChannels         ͼ���ͨ����
**		FGDStatModel* _model  ��Ҷ˹״̬ģ��
**	����ֵ��
**		�� 
**	˵��:
**		ͨ��������ͼ��ʱ�ĸ��µ�ǰ��Ҷ˹ģ��
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
                    if( !stat->is_trained_dyn_model ) val = 1; //û�м�⵽�˶�����֮ǰ
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
                else if( stat->is_trained_st_model ) //��⵽�˶�������
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
	//1.��������ͳ�Ʊ�
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
                else                        //����ڶ�̬�����ҵ��������أ�ʹ����Ӧ�ĸ���ֵ����
                {//update
                    PV_CC(indx) += alpha;
                    if( !((unsigned char*)model->foreground)[i*mask_step+j] )
                    {
                        PVB_CC(indx) += alpha;
                    }
                }

                //re-sort CCt table by Pv  //�Զ�̬����������
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
			//2. ���²ο�����ͼ��
            if( !((unsigned char*)model->foreground)[i*mask_step+j])//�Ա�������и���
            {
                unsigned char* ptr = ((unsigned char*)model->background) + i*WidthStep+j*nChannels;
                
                if( !((unsigned char*)model->Ftd)[i*mask_step+j] &&
                    !((unsigned char*)model->Fbd)[i*mask_step+j] )//�Ի����𽥱仯�ı�������и���
                {
                    //apply IIR filter
                    for( l = 0; l < nChannels; l++ )
                    {
                        int a = int(ptr[l]*(1 - model->params.alpha1) + model->params.alpha1*curr_data[l]);
                        ptr[l] = (unsigned char)a;
                    }
                }
                else                                              //���˶��ı����㣨��ҡҷ����Ҷ��������˸��ˮ�ƣ����и���
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
**	�������ƣ�
**		ReleaseFGDStatModel()
**	����
**		FGDStatModel* _model  ��Ҷ˹״̬ģ��
**	����ֵ��
**		�� 
**	˵��:
**		�������ͷű�Ҷ˹ģ�������Դ
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
/*                    C.���ڼ򵥱�����ģ���˶�Ŀ����                  */
/************************************************************************/

/*************************************************************************
**	�������ƣ�
**		SetBgHistgram()
**	����
**		short int *pHistgram  ����ͳ�Ƶı���ֱ��ͼ
**      BYTE *sGray           ��ǰ�����ͼ��
**		int nWidth            ͼ��Ŀ��
**		int nHeight           ͼ��ĸ߶�
**	����ֵ��
**		�� 
**	˵��:
**		��������ʱ��ͳ��ֱ��ͼ
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
**	�������ƣ�
**		GetPeakBgByHistgram()
**	����
**		short int *pHistgram  ͳ�Ƶı���ֱ��ͼ
**      BYTE *sBg             ���Ƶı���ͼ��
**		int nWidth            ͼ��Ŀ��
**		int nHeight           ͼ��ĸ߶�
**	����ֵ��
**		�� 
**	˵��:
**		ͨ��ʱ��ֱ��ͼ������ֵ����
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
**	�������ƣ�
**		GetMedianBgByHistgram()
**	����
**		short int *pHistgram  ͳ�Ƶı���ֱ��ͼ
**      BYTE *sBg             ���Ƶı���ͼ��
**		int nWidth            ͼ��Ŀ��
**		int nHeight           ͼ��ĸ߶�
**      int countFrames       �Ѿ������֡��
**	����ֵ��
**		�� 
**	˵��:
**		ͨ��ʱ��ֱ��ͼ������ֵ����
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
**	�������ƣ�
**		GetRuningAvgBg()
**	����
**		short int *pHistgram  ͳ�Ƶı���ֱ��ͼ
**      BYTE *sBg             ���Ƶı���ͼ��
**		int nWidth            ͼ��Ŀ��
**		int nHeight           ͼ��ĸ߶�
**      int countFrames       �Ѿ������֡��
**	����ֵ��
**		�� 
**	˵��:
**		ͨ��ʱ��ֱ��ͼ������ֵ����
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
**	�������ƣ�
**		CreateSimpleBGModel()
**	����
**      unsigned char*Image   ��ʼ֡ͼ��
**		int Width             ͼ��Ŀ��
**		int Height            ͼ��ĸ߶�
**      int Thresh            ����ֵ����ֵ
**		int SimBGType         �򵥱�����ʽ1: ��ֵ���� 2����ֵ���� 3��ֱ��ͼ��ֵ����
**	����ֵ��
**		�򵥱���ģ�� 
**	˵��:
**		��ʼʱ�����򵥱���ͳ��ģ����ز���
************************************************************************/
SimpleBGModel* CreateSimpleBGModel( unsigned char*Image,int Width,int Height,int Thresh,int SimBGType)
{
	//��SimpleBGModel�����ڴ�
	SimpleBGModel* p_model = 0;
    p_model = (SimpleBGModel*)malloc( sizeof(*p_model) );
    memset( p_model, 0, sizeof(*p_model) );

	//���SimpleBGModel
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
**	�������ƣ�
**		UpdateSimpleModel()
**	����
**      unsigned char*Image    ��ǰҪ�����ͼ��
**		int Width              ͼ��Ŀ��
**		int Height             ͼ��ĸ߶�
**		SimpleBGModel* _model  �򵥱���ģ��
**	����ֵ��
**		�� 
**	˵��:
**		ͨ��������ͼ��ʱ�ĸ��µ�ǰ�򵥱���ģ��
************************************************************************/
void UpdateSimpleModel( unsigned char*Image,int Width,int Height,SimpleBGModel* model )
{
	if (model->SimBGType==SIMPLEBG_AVG)                        //����ʱ���ֵ��������
	{
		GetRuningAvgBg(model->background,Image,Width,Height,model->countFrames);
		InfraredSubtract(model->foreground,Image,model->background,Width,Height,model->thresh);
		++model->countFrames;
	}
	else if (model->SimBGType==SIMPLEBG_MEDIAN)                    //����ʱ����ֵ��������
	{
		++model->countFrames;
		SetBgHistgram(model->pBackHist,Image,Width,Height);
		GetMedianBgByHistgram(model->pBackHist,model->background,Width,Height,model->countFrames);
		InfraredSubtract(model->foreground,Image,model->background,Width,Height,model->thresh);
	}
	else                                             //����ֱ��ͼ��ֵ��������
	{
		SetBgHistgram(model->pBackHist,Image,Width,Height);
		GetPeakBgByHistgram(model->pBackHist,model->background,Width,Height);
		InfraredSubtract(model->foreground,Image,model->background,Width,Height,model->thresh);
	}
}
/*************************************************************************
**	�������ƣ�
**		ReleaseSimpleModel()
**	����
**		SimpleBGModel* _model  �򵥱���ģ��
**	����ֵ��
**		�� 
**	˵��:
**		�������ͷż򵥱���ģ�������Դ
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
/*                                 ����                                 */
/************************************************************************/ 
/************************************************************************/
/*                    D.���ں��ܶȵı�����ģ���˶�Ŀ����          */
/************************************************************************/ 




/*************************************************************************
**	�������ƣ�
**		excle()
**	����
**       i   ֵΪi������ 
**	
**	����ֵ��
**		�� 
**	˵��:
**		�������ұ�ӿ�����ٶ�
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
**	�������ƣ�
**		createkdemodel()
**	����

SampleImage   ǰ10֡Ϊ����ͼ�����һ֡Ϊ��ǰ�����ͼ�񣬹�11֡
height        ÿһ֡ͼ�������
lineByte      ÿ֡֡ͼ�������
N_L           ����������֡��+��ǰ֡��N_L֡
exclee        ������
fore          ����Ŀ����ȡ���
thre          �ֶ�������ֵ
	
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
					double thre)					   
{  
	int i,j,k,index;
	int SampleNum=N_L-1; // N_Lȡ11��SampleNumȡ10


    unsigned char* foreground=new unsigned char[height*lineByte];
	memset( foreground, 0, sizeof(unsigned char)*height*lineByte);
	
    memcpy(foreground, SampleImage[SampleNum], sizeof(unsigned char)*height*lineByte);//����ǰ֡���Ƶ�����foreground�У�����


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
	

	//*****************��ֵ**********************************
	
	
	for (j=0;j<lineByte*height;j++)
	{
		if (ab[j]<thre)
		{
			foreground[j]=255;
		}
		
	}
 
	//��ֵ��
	for (j=0;j<lineByte*height;j++)
	{
		if (foreground[j]<255)
		{
			fore[j]=0;
		}
		else
			fore[j] = 255;
	}


	


 	//����ȥ��������ĺ�����foreȥ��������

	delete []foreground;
	foreground = NULL;

 
	delete []ab;
	ab = NULL;

 }






/************************************************************************/
/*                                 ����                                 */
/************************************************************************/ 

