/************************************************************************
* ��Ȩ���� (C)2011, ���ӿƼ���ѧ.ͼ��ͼ�����źŴ���Ӧ���о���
* 
* �ļ����ƣ�  preproc.c
* �ļ���ʶ��  
* ����ժҪ��  ͼ��Ԥ����ģ��Դ�ļ�
* ����˵����  
* ��ǰ�汾��  
* ��    �ߣ�  
* ������ڣ�  
* 
************************************************************************/
#include "../preproc/preproc.h"
//#include "../preproc/edge_canny.h"


/**************************************************************************
*                             ����                                        *
**************************************************************************/


/**************************************************************************
*                             ��                                          *
**************************************************************************/


/**************************************************************************
*                          ��������                                       *
**************************************************************************/


/**************************************************************************
*                           ȫ�ֱ���                                      *
**************************************************************************/

/**************************************************************************
*                           �ֲ�����ԭ��                                  *
**************************************************************************/

 
/**************************************************************************
*                       ȫ�ֺ���ʵ��                                      *
**************************************************************************/
/**************************************************************************
* �������ƣ� LocationProc
* ���������� location msg ����
* ��������� 
* ��������� 
* �� �� ֵ�� 
* ����˵���� 
* �޸�����    �汾��     �޸���      �޸�����
* -----------------------------------------------
* 2011/09/01   V1.0       XXXX          XXXX
**************************************************************************/

/************************************************************************/
/*                                                                      */
/*                    һ������ͼ��Ǿ���У��ģ��                        */
/*                                                                      */
/************************************************************************/
/************************************************************************
**  �������ƣ�
**      TwopointNUC()
**	����
**	    unsigned char *image    ָ��ͼ���ָ��
**		int width               ͼ��Ŀ��
**		int height              ͼ��ĸ߶�
**      double * G              �������
**      double * O              ƫ�þ���
**	����ֵ��
**	    ��
**	˵����
**      ����Ǿ���У���㷨      
*************************************************************************/
void TwopointNUC(unsigned char * image,int width,int height,double * G , double * O )
{
	int x,y;
	long lLBytes = (width+3 )/4*4;
	unsigned char * temp = new unsigned char [lLBytes*height];
	memcpy(temp,image,sizeof(unsigned char )*lLBytes*height);
	for (y = 0;y<height;y++)
	{
		for (x = 0;x<width;x++)
		{
			temp[lLBytes*y+x] = G[y*width+x]*temp[lLBytes*y+x]+O[y*width+x];
		}
	}
	double maxT = 1.0;
	double offset = image[0] - temp[0];
	for( y=0;y<height;y++)
	{
		for(int x=0;x<width;x++)
		{   
			if ((temp[lLBytes*y+x]+offset)>maxT)
				maxT=temp[lLBytes*y+x]+offset;       // ����ͼ����������ֵ
		}
	}
	
	for(y=0;y<height;y++)
	{
		for(int x=0;x<width;x++)
		{   
			image[lLBytes*y+x]=(int)((temp[lLBytes*y+x]+offset)*255/maxT);//��У��ͼ����й�һ��
		}
	}
	delete [] temp;
	
}

/************************************************************************
**  �������ƣ�
**      HighpassNUC()
**	����
**	    unsigned char *image    ָ��ͼ���ָ��
**		int width               ͼ��Ŀ��
**		int height              ͼ��ĸ߶�
**      double *F               ָ�򾭹���ͨ�˲�������������,2ά˫�������飬width��height
**      double M                ʱ�䳣��M
**	����ֵ��
**	    ��
**	˵����
**      ʱ���ͨ�˲��Ǿ���У���㷨      
*************************************************************************/
void HighpassNUC(unsigned char* image , int width , int height, double * F,double M)
{
	long    lLBytes = (width*8+31)/32*4;
	unsigned char * Temp = new unsigned char [lLBytes* height];
	memcpy(Temp, image, sizeof(unsigned char) *lLBytes*height);

	////��ǰ֡�ĻҶ�ƽ�����������Ȳ���
	double Avg = 0;
	for (int y = 0;y<height;y++)
	{
		for (int x = 0;x<width;x++)
		{
			
			Avg = Avg+Temp[width*y+x];
		}
	}
	Avg = Avg/(width*height);
	
	
	for ( y =0;y<height;y++ )	
	{
		for (int x = 0;x<width;x++)
		{
			F[y* width +x] = (1/M)*Temp[y*lLBytes +x] + (1-1/M)*F[y* width +x];
			Temp[y*lLBytes +x] = (int) (1.1*Temp[y*lLBytes +x] - F[y* width +x] );	
		}
	}
	////�Խ��ͼ����й�һ��������Աȶ�
	double maxT = 1.0;
	
	for( y=0;y<height;y++)
	{
		for(int x=0;x<width;x++)
		{   
			if (Temp[lLBytes*y+x]>maxT)
				maxT=Temp[lLBytes*y+x];                         // ����ͼ����������ֵ
		}
	}
	
	for(y=0;y<height;y++)
	{
		for(int x=0;x<width;x++)
		{   
			image[lLBytes*y+x]=(int)(Temp[lLBytes*y+x]*255/maxT);//�Բ�ֵͼ����й�һ��
		}
	}
	delete [] Temp;	
}

/***********************************************************************
**  �������ƣ�
**      ANNNUC()
**  ����
**      unsigned char * image    ָ��ͼ���ָ��
**	    int width                ͼ��Ŀ��
**	    int height               ͼ��ĸ߶�
**      float *G                 ָ�����ڵ�������������,2ά˫��������width��height
**      float *O                 ָ�����ڵ�����ƫ������,2ά˫��������width��height
**  ����ֵ��
**      ��
**  ˵����
**      ����������ķǾ���У���㷨                                               
*************************************************************************/
void ANNNUC(unsigned char * image, int width ,int height,float * G,float *O)
{
	long   lLBytes = (width*8+31)/32*4;
	float u1 = 0.000001f;
	float u2 = 0.000002f;
	
	float f = 0;
	unsigned char * Temp = new unsigned char[lLBytes*height];
	unsigned char *TempIm = new unsigned char[lLBytes*height];
	memcpy(Temp,image,sizeof(unsigned char)*lLBytes*height);
	memcpy(TempIm,image,sizeof(unsigned char)*lLBytes*height);
	for (int y=0;y<height;y++)
	{
		for (int x =0;x<width;x++)
		{
			Temp[y*lLBytes+x] = (unsigned char) (G[y*width+x]*TempIm[y*lLBytes+x]+O[y*width+x]);
		}
	}

	for ( y=1;y<height-1;y++)
	{
		for(int x=1;x<width-1;x++)
		{
			f=(Temp[(y-1)*lLBytes+x]+Temp[(y+1)*lLBytes+x]+Temp[y*lLBytes+(x+1)]+Temp[y*lLBytes+(x-1)]+3*Temp[lLBytes*y+x])/7;
			G[y*width+x] = G[y*width+x] + 2*u1*TempIm[lLBytes*y+x]*(f-Temp[lLBytes*y+x]);
			O[y*width+x] = O[y*width+x] + 2*u2*(f-Temp[lLBytes*y+x]);
		}
	}
	float maxT = 1.0;
	
	for( y=0;y<height;y++)
	{
		for(int x=0;x<width;x++)
		{   
			if (Temp[lLBytes*y+x]>maxT)
				maxT=Temp[lLBytes*y+x];       // ����ͼ����������ֵ
		}
	}
	
	for(y=0;y<height;y++)
	{
		for(int x=0;x<width;x++)
		{   
			image[lLBytes*y+x]=(unsigned char)(Temp[lLBytes*y+x]*255/maxT);//��ͼ����й�һ��
			
		}
	}
	delete [] Temp;
	delete [] TempIm;
}

/************************************************************************
**  �������ƣ�
**      MatrixMutiply()
**  ����
**      double* pOutput  ָ���������2ά˫��������l��n
**      double* pInA     ָ���������A,2ά˫��������l��m
**      double* pInB     ָ���������B,2ά˫��������m��n
**      int l            �������A ����
**      int m            �������A���� (�������B����)
**      int n            �������B����
**  ����ֵ��
**      ��   
**  ˵����
**      ʵ��l��m ����A��m��n����B�����                        
*************************************************************************/
void MatrixMultiply(double* pOutput, double* pInA, double* pInB, int l, int m, int n)
{
	memset(pOutput, 0, sizeof(double) * l * n);
	int i, j, k;
	for(i = 0; i < l; i++)
		for(j = 0; j < n; j++)
		{
			for(k = 0; k < m; k++)
			{
				pOutput[i * n + j] += pInA[i * m + k] * pInB[k*n + j];
			}
		}
}

/************************************************************************
**  �������ƣ�
**      MatrixAdd()
**  ����
**     double * pOutput  ָ���������2ά˫��������m��n
**     double * pInA     ָ���������A��2ά˫��������m��n
**     double * pInB     ָ���������B��2ά˫��������m��n
**     int m             �������A��B������ 
**     int n             �������A��B������
**  ����ֵ��
**     ��
**  ˵����
**     ʵ��2��m��n�ľ������                                                                       
**************************************************************************/
void MatrixAdd(double *pOutput, double *pInA,double *pInB, int m, int n)
{
	//	memset(pOutput,0,sizeof(float )*m*n);
	for (int i =0 ; i <m ; i++)
	{
		for(int j= 0 ; j<n ; j++)
		{
			pOutput[i*n+j] = pInA[i*n+j] + pInB[i*n+j];
		}
	}
	
}

/************************************************************************
**   ��������
**       Static_KalmanNUC()
**   ����
**       unsigned char * image   ָ������ͼ��
**       int width               ͼ����
**       int height              ͼ��߶�
**       double **X              ָ�����ڵ�������ģ�Ͳ����������ָ�룬
**                               2ά˫��������width��height 
**   ����ֵ��
**       ��
**   ˵����
**      ������̬kalman�˲��ķǾ���У���㷨                                                                   
*************************************************************************/
void Static_KalmanNUC(unsigned char * image , int width ,int height , double ** X)
{
	long    lLBytes = ((width*8)+31)/32*4;
	double * TempX = new double [2];
	unsigned char * Temp = new unsigned char [lLBytes*height];
	memcpy(Temp,image,sizeof(unsigned char)*lLBytes*height);
	//������ʼ��
    double alpha = 0.9;
	double beta = 0.9;
	double sigmaA = 0.25;
	double sigmaB =0.1;
	
	double *A = new double[4];
	A[0] = alpha;
	A[1] = 0;
	A[2] = 0;
	A[3] = beta;
	double R =5.5;
	double *D = new double[4];
	D[0] = 1-alpha;
	D[1] = 0;
	D[2] = 0;
	D[3] = 1-beta;
	double *x0 = new double[2];
	x0[0] = 1;
	x0[1] = 0;
    double * M = new double[4];
	MatrixMultiply(M,D,x0,2,2,1);
	double b =( 1-beta*beta)*(sigmaB-R);
	double delta = b*b+4*R*( 1-beta*beta)*sigmaB ;
	double d  = (b+sqrt(delta))/2;
	
	double * K = new double[2];
	K[0] = 0;
	K[1] = d/(d+R);
    double T =5;
	double sigmaT =5.5;
	double deltaV = 4.5;
	double w;
	double bk;
	double * TempK = new double [2];
	memset(TempK,0,sizeof(double)*2);
    for(int y = 0;y<height;y++)
	{
		for(int x = 0;x<width;x++)
		{
			memcpy(TempX,X[y*width + x],sizeof(double)*2);
			MatrixMultiply(TempX,A,TempX,2,2,1);
			MatrixAdd(TempX,TempX,M,2,1);
			
			double TempY = image[lLBytes*y+x] -TempX[1];
			TempK[0] =0;
			TempK[1] = TempY*K[1];
			MatrixAdd(TempX,TempX,TempK,2,1);
			w = TempX[0]*sigmaT/(TempX[0]+deltaV);
			bk = T-w*(TempX[0]*T +TempX[1]);
            Temp[lLBytes*y+x] = w*Temp[lLBytes*y+x] + bk;
			memcpy(X[y*width + x],TempX,sizeof(double)*2);
		}
	}
	double maxT = 1.0;
	
	for( y=0;y<height;y++)
	{
		for(int x=0;x<width;x++)
		{   
			if (Temp[lLBytes*y+x]>maxT)
				maxT=Temp[lLBytes*y+x];       // ����ͼ����������ֵ
		}
	}
	
	for(y=0;y<height;y++)
	{
		for(int x=0;x<width;x++)
		{   
			image[lLBytes*y+x]= (unsigned char)(Temp[lLBytes*y+x]*255/maxT);//�Բ�ֵͼ����й�һ��
		}
	}
	delete [] Temp;
	delete [] TempX;
	delete [] TempK;
	delete [] K;
	delete [] M;
	delete [] x0;
	delete [] D;
	delete [] A;
	
}

/************************************************************************/
/*                                                                      */
/*                    ��������ͼ��������������ǿģ��                    */
/*                                                                      */
/************************************************************************/

/*************************************************************************
**	�������ƣ�
**		Medianfilter_3x3()
**	����
**		unsigned char* image  ָ��ͼ���ָ��
**		int width             ͼ��Ŀ��
**		int height            ͼ��ĸ߶�
**	����ֵ��
**		�� 
**	˵��:
**		ͼ����ֵ�˲���ȥ����������
************************************************************************/
void Medianfilter_3x3(unsigned char* image,int width,int height)
{
	int i,j,k,l,count;
	int ii,mid;
	BYTE mean[9];
	int ImageWidth, ImageHeight;
	ImageWidth = width;
	ImageHeight = height;
	int lineByte=(ImageWidth+3)/4*4;
	unsigned char* pBuf=new unsigned char[lineByte*ImageHeight]; 
	memcpy(pBuf,image,sizeof(unsigned char)*lineByte*ImageHeight);
	for(i=1;i<ImageHeight-1;i++)
	{
		for(j=1;j<ImageWidth-1;j++)
		{
			count = 0;	
			for(k=-1;k<=1;k++)
			{
				for(l=-1;l<=1;l++)
				{
					mean[count]=image[(k+i)*lineByte+l+j];
					count++;
				}
			}
			for(k=0;k<9;k++)
			{
				mid=mean[k];
				ii=k;
				for(l=k+1;l<9;l++)
				{
					if(mid>mean[l])
					{
						mid=mean[l];
						ii=l;
					}
				}
				mean[ii]=mean[k];
				mean[k]=mid;
			}
			pBuf[i*lineByte+j]= mean[4];
		}
	}
	memcpy(image,pBuf,sizeof(unsigned char)*lineByte*ImageHeight);
	delete []pBuf;
}
/*************************************************************************
**	�������ƣ�
**		Near_Mean()
**	����
**		unsigned char* image  ָ��ͼ���ָ��
**		int width             ͼ��Ŀ��
**		int height            ͼ��ĸ߶�
**      int Ksize             ģ��Ĵ�С
**	����ֵ��
**		�� 
**	˵��:
**		ͼ���ֵ�˲���ȥ�����������ƽ��ͼ��
************************************************************************/
void Near_Mean(unsigned char* image,int width,int height, int Ksize)
{
	int lLBytes = ((width*8)+31)/32*4;
	int k2 = (Ksize-1)/2;
	int i,j,x,y;
	float mean;
	
	unsigned char* pBuf=new unsigned char[lLBytes* height]; 
	memcpy(pBuf,image,sizeof(unsigned char)*lLBytes* height);
	for(y=k2; y<height-k2; y++)
	{
		for(x=k2; x<width-k2; x++)
		{
			mean=0.0;
			for(j=-k2; j<=k2; j++)
			{
				for(i=-k2; i<=k2; i++)
					mean += (unsigned char)image[(y+j)*lLBytes+x+i];
			}
			mean /=(float)(Ksize*Ksize);
			pBuf[y*lLBytes+x] = (unsigned char)(mean+0.5);
		}
	}
	memcpy(image,pBuf,sizeof(unsigned char)*lLBytes* height);
	delete [] pBuf;
}
/*************************************************************************
**	�������ƣ�
**		Gaussian()
**	����
**		unsigned char* image  ָ��ͼ���ָ��
**		int width             ͼ��Ŀ��
**		int height            ͼ��ĸ߶�
**      int Ksize             ģ��Ĵ�С
**      int sigma             ��˹�˲���������
**	����ֵ��
**		�� 
**	˵��:
**		ͼ���˹�˲���ȥ����˹�����������ƽ��ͼ��
************************************************************************/
void Gaussian(unsigned char* image,int width,int height, int  Ksize , float sigma)
{
	/* ������˹�� */
	int i,j,x,y,k2;
	k2    = (Ksize-1)/2;
	//int lLBytes = ((width*8)+31)/32*4;
	float *kernel= (float*)malloc(Ksize*Ksize*sizeof(float));
	/* �����˹�˲��� */
	for (i=0; i<Ksize; i++) 
	{
		for (j=0;j < Ksize;j++)
		{
			int   tmp1 = (i-k2)*(i-k2) + (j-k2)*(j-k2);
			float tmp2 = (-1*tmp1) / (2*sigma*sigma);
			kernel[i*Ksize+j] = (float)exp(tmp2);
		}
	}
	/* ��һ����˹�˲���*/
	float sum = 0;
	for (i=0;i < Ksize;i++)
	{
		for (j=0;j < Ksize;j++)
			sum = sum + kernel[i*Ksize+j];
	}
	for (i=0;i < Ksize;i++)
	{	
		for (j=0;j < Ksize;j++)
			kernel[i*Ksize+j] = kernel[i*Ksize+j] / sum;
	}
	
	//���и�˹�˲�
	/************************************************************************/
	/* lBytes=width                                                                     */
	/************************************************************************/
	unsigned char* pBuf=new unsigned char[width* height]; 
	memcpy(pBuf,image,sizeof(unsigned char)*width* height);
	float fsum;
	for(y=k2; y<height-k2; y++)
	{
		for(x=k2; x<width-k2; x++)
		{
			fsum=0.0;
			for(j=-k2; j<=k2; j++)
			{
				for(i=-k2; i<=k2; i++)
					fsum += image[(y+j)*width+x+i]*kernel[(j+k2)*Ksize+i+k2];
			}
			if(fsum > 255) 
				pBuf[y*width+x] = (unsigned char)255;
			else
				pBuf[y*width+x] = (unsigned char)(fsum+0.5);
		}
	}
	memcpy(image,pBuf,sizeof(unsigned char)*width* height);
	free(kernel);
	delete [] pBuf;
}
/*************************************************************************
**	�������ƣ�
**		DifGaussian()
**	����
**		unsigned char* image  ָ��ͼ���ָ��
**		int width             ͼ��Ŀ��
**		int height            ͼ��ĸ߶�
**      int Ksize             ģ��Ĵ�С
**	����ֵ��
**		�� 
**	˵��:
**		ͼ��������Ը�˹�˲�(�Ľ��˸�˹�˲�)��ȥ����˹�����������ƽ��ͼ��
************************************************************************/
void DifGaussian(unsigned char* image,int width,int height, int  Ksize)
{
	int     i,j,x,y;
	float	fsum;
	int     k2 = (Ksize-1)/2;
	int     lLBytes = ((width*8)+31)/32*4; 
    
	float divx,divy;
	float sita;
	float sum1,sum2;
	float sigmax2,DS;
	float sigma;
	float k=20;      //��������
    float average=0;
    float sigmax=0;
    float sigmay=0;
	
	unsigned char* pBuf=new unsigned char[lLBytes* height]; 
	float *kernel=new float[Ksize*Ksize];
	for( y=k2; y<height-k2; y++) 
	{
		for( x=k2; x<width-k2; x++)
		{			
			//1.����ÿһ���Ӧ��ģ��
			sigmax2=float(255)/(unsigned char)(image[y*lLBytes+x]);  
			float sum=0;
			for ( j=-k2; j<=k2; j++)
			{
				for ( i=-k2; i<=k2; i++) 
				{
					
					sum+= (unsigned char)image[(y+j)*lLBytes+x+i];
				}
			}
			//����������
			average=sum/float(Ksize*Ksize);
			sum1=0;
			for ( j=-k2; j<=k2; j++)
			{
				for ( i=-k2; i<=k2; i++) 
				{    
					sum1+= ((unsigned char)image[(y+j)*lLBytes+x+i]-average)*
						((unsigned char)image[(y+j)*lLBytes+x+i]-average);
				}
			}
			DS=sum1/float(Ksize*Ksize);
			sigmax=(float)sqrt(sigmax2);
			sigmay=k/(k+DS)*sigmax;
			
			//����ÿ��������Ӧ�ĽǶ�
			divx=float(image[y*lLBytes+x+1]-image[y*lLBytes+x]);
			divy=float(image[(y+1)*lLBytes+x]-image[y*lLBytes+x]);
			float l=divy/divx;
			sita=float(atan(l)+PI/2);
			sigma=20;
			for (i=0; i<Ksize; i++)
			{
				for (j=0;j < Ksize;j++) 
				{
					int   tmp1 = (i-k2)*(i-k2) + (j-k2)*(j-k2);
					float tmp2 = (-1*tmp1) / (2*sigma*sigma);
					kernel[i*Ksize+j] = (float)exp(tmp2);
				}
			}
			
			//2.���й�һ��
			sum2=0;
			for(j=0;j<Ksize;j++)
			{
				for(i=0;i<Ksize;i++)
					sum2+=kernel[j*Ksize+i];
			}
			for(j=0;j<Ksize;j++)
			{
				for(i=0;i<Ksize;i++)
					kernel[j*Ksize+i]=kernel[j*Ksize+i]/sum2;
			}
			
			//3.���о��
            fsum = 0;
			for ( j=-k2; j<=k2; j++)
			{
				for ( i=-k2; i<=k2; i++) 
				{
					
					fsum += image[(y+j)*lLBytes+x+i] * kernel[(j+k2)*Ksize+i+k2];
				}
			}
			if(fsum > 255) 
				pBuf[y*lLBytes+x] = (unsigned char)255;
			else
				pBuf[y*lLBytes+x] = (unsigned char)(fsum+0.5);
		}
	}
	memcpy(image,pBuf,sizeof(unsigned char)*lLBytes* height);
	delete [] kernel;
	delete [] pBuf;
}

int  level=9;					//ͼ�����ֽ�߶�					

struct  array2_2_N_N			//�ṹ�壬����״̬ת�Ƹ���
{
	double x1,x2,x3,x4;
};

struct array2_N_N				//�ṹ�壬��������״̬�ĸ���
{
	double x1,x2;
};


/*************************************************************************
**	�������ƣ�
**		hmtmodel()       
**	����
**		double si[][level]			��������ͼ��ķ���ϵ��ģ��	
**		double es[][2][level]		��������ͼ��ĸ���ת�ƾ���ģ��
**	����ֵ��
**		�� 
**	˵��:
**		��������ͼ���ģ��ϵ��	
************************************************************************/
void  hmtmodel(double **si, double ***es)
{
	//�������ָ�˹�����Ĳ���
	double 	alpha_big = 2.5;			//�󷽲��˥������	
	double 	C1_big = 13;
	double 	alpha_sm = 2.5;             //С�����˥������
	double 	C1_sm = 7;
	double	beta = 1;
	
	//���������ϵ�����㲻ͬ�߶��µķ���ϵ��
	int n1,n2,n3;
	for (n2=0;n2<level;n2++)
	{
		si[0][n2]=(double)(pow(2,C1_sm)*pow(2,-alpha_sm*(n2+1)));
		si[1][n2]=(double)(pow(2,C1_big)*pow(2,-alpha_big*(n2+1)));
	}
	
	//���������ϵ�����㲻ͬ�߶��µĳ�ʼ���ʷֲ�
	double *p00=new double[level];
	double *p10=new double[level];
	double *p11=new double[level];
	memset(p00,0,level*sizeof(double));
	memset(p10,0,level*sizeof(double));
	memset(p11,0,level*sizeof(double));
	
	for (n1=0;n1<3;n1++)
	{
		p00[n1]=1;
		p10[n1]=0;
		p11[n1]=1;
	}
	
	for (n1=3;n1<level;n1++)
	{
		p00[n1]=(double)(0.8+0.2*(1-pow(2,-(beta*(n1-3)))));
		p10[n1]=(double)(1-p00[n1]);
		p11[n1]=(double)(0.9-0.4*(1-pow(2,-(beta*(n1-3)))));
	}
	
	//���㲻ͬ�߶��µĳ�ʼת�Ƹ���
	for (n3=0;n3<level;n3++)
	{
		es[0][0][n3]=1-p10[n3];
		es[0][1][n3]=1-p11[n3];
		es[1][0][n3]=p10[n3];
		es[1][1][n3]=p11[n3];
	}

	delete[]p00;
	delete[]p10;
	delete[]p11;
}

/*************************************************************************
**	�������ƣ�
**		 CompareGray()       
**	����
**		const void *a		����Ƚϲ���1
**		const void *b		����Ƚϲ���2
**	����ֵ��
**		���� 1���� -1
**	˵��:
**		��ֵ����ıȽϺ�����������ȡС��ϵ������ֵ
************************************************************************/ 
int CompareGray(const void *a, const void *b)
{
	return  *(double*)a>*(double*)b ?1:-1;
}

/*************************************************************************
**	�������ƣ�
**		 posthh()       
**	����
**		double			 *wavelete	С���ֽ�ϵ��ָ��	
**		array2_2_N_N	 *ES		ת�Ƹ���ģ��ϵ��ָ��
**		array2_N_N		 *SI		����ģ��ϵ��ָ��
**		array2_N_N		 *P1		�洢������ĺ������ָ��
**	����ֵ��
**		��
**	˵��:
**		����Խ�����С��ϵ���ĺ�����ʣ���rice��ѧ�ṩ��MATLAB����ı�
************************************************************************/
void  posthh(double *wavelete,array2_2_N_N *ES,array2_N_N *SI,array2_N_N *P1,int N)
{
	array2_N_N *BE=(array2_N_N *)malloc(N*N*sizeof(array2_N_N));
	memset(BE,0,N*N*sizeof(array2_N_N));
	array2_N_N *BEP=(array2_N_N *)malloc(N*N*sizeof(array2_N_N));
	memset(BEP,0,N*N*sizeof(array2_N_N));
	array2_N_N *BER=(array2_N_N *)malloc(N*N*sizeof(array2_N_N));
	memset(BER,0,N*N*sizeof(array2_N_N));
	array2_N_N *wtmp=(array2_N_N *)malloc(N*N*sizeof(array2_N_N));
	memset(wtmp,0,N*N*sizeof(array2_N_N));
	array2_N_N *scale=(array2_N_N *)malloc(N*N*sizeof(array2_N_N));
	memset(scale,0,N*N*sizeof(array2_N_N));
	array2_N_N *AL=(array2_N_N *)malloc(N*N*sizeof(array2_N_N));
	memset(AL,0,N*N*sizeof(array2_N_N));

	for ( int n1=0;n1<N*N;n1++)
	{
		double temp1=1.0/sqrtf(2*PI*SI[n1].x1)*exp(-pow(wavelete[n1],2)/(2*SI[n1].x1))+0.0000000000001;
		wtmp[n1].x1=(double)temp1;
		double temp2=1.0/sqrtf(2*PI*SI[n1].x2)*exp(-pow(wavelete[n1],2)/(2*SI[n1].x2))+0.0000000000001;
		wtmp[n1].x2=(double)temp2;
		scale[n1].x1=scale[n1].x2=(double)((temp1+temp2)/2);
	}

	for (int si=(int)(pow(2,level-1));si<N;si++)
	{
		for (int sj=(int)(pow(2,level-1));sj<N;sj++)
		{
			int k=si*N+sj;
			BE[k].x1=wtmp[k].x1/scale[k].x1;
			BE[k].x2=wtmp[k].x2/scale[k].x2;
		}
	}

	//���ϲ���
	for (int k=level;k>1;k--)
	{
		int J=(int)(pow(2,k-1));
		int J2=J*J;
		
		double k1=ES[J*N+J].x1;
		double k2=ES[J*N+J].x2;		
		double k3=ES[J*N+J].x3;
		double k4=ES[J*N+J].x4;
		
		for (si=J;si<2*J;si++)
		{
			for (int sj=J;sj<2*J;sj++)
			{				
				int k=si*N+sj;	
				BEP[k].x1=k1*BE[k].x1+k3*BE[k].x2;
				BEP[k].x2=k2*BE[k].x1+k4*BE[k].x2;
			}
		}
		//����beta���Ӿ����BE����
		array2_N_N *BCtmp=(array2_N_N*)malloc(J*J/4*sizeof(array2_N_N));
		memset(BCtmp,0,J*J/4*sizeof(array2_N_N));

		//һ��ѭ����ֻ���±ꡣ
		for (si=J;si<2*J;si+=2)
		{
			for (int sj=J;sj<2*J;sj+=2)
			{				
				int k1=si*N+sj;	
				int k2=k1+1;	
				int k3=k1+N;	
				int k4=k3+1;	
				int m1=si/2*N+sj/2;			
			
				double temp1=BEP[k1].x1*BEP[k2].x1*
							BEP[k3].x1*BEP[k4].x1;	
			 	double temp2=BEP[k1].x2*BEP[k2].x2*
 					    	BEP[k3].x2*BEP[k4].x2;	
				BE[m1].x1=temp1*wtmp[m1].x1/(scale[m1].x1*((temp1+temp2)/2));
				BE[m1].x2=temp2*wtmp[m1].x2/(scale[m1].x2*((temp1+temp2)/2));
		
				BER[k1].x1=BE[m1].x1/BEP[k1].x1;
				BER[k2].x1=BE[m1].x1/BEP[k2].x1;
				BER[k3].x1=BE[m1].x1/BEP[k3].x1;
				BER[k4].x1=BE[m1].x1/BEP[k4].x1;

				BER[k1].x2=BE[m1].x2/BEP[k1].x2;
				BER[k2].x2=BE[m1].x2/BEP[k2].x2;
				BER[k3].x2=BE[m1].x2/BEP[k3].x2;
				BER[k4].x2=BE[m1].x2/BEP[k4].x2;
			}
		}
 		free(BCtmp);
	}

	//DOWN_step
	//��AL����ֵ
	AL[1*N+1].x1=22.627416997970;
	AL[1*N+1].x2=1448.154687870049;

	for (k=2;k<=level;k++)
	{
		int J=int(pow(2,k-1));
		int J2=J*J;

		array2_N_N *Atmp1=(array2_N_N*)malloc(J*J*sizeof(array2_N_N));
		memset(Atmp1,0,J*J*sizeof(array2_N_N));

		for (int k1=J/2;k1<J;k1++)
		{
			for (int k2=J/2;k2<J;k2++)
			{
				int m1=k1*N+k2;
				int n1=2*((k1-J/2)*J+k2-J/2);
				int n2=n1+1;
				int n3=n1+J;
				int n4=n3+1;
				Atmp1[n1].x1=Atmp1[n2].x1=Atmp1[n3].x1=Atmp1[n4].x1=AL[m1].x1;
				Atmp1[n1].x2=Atmp1[n2].x2=Atmp1[n3].x2=Atmp1[n4].x2=AL[m1].x2;
			}
		}

		for (si=0;si<J;si++)
		{
			for (int sj=0;sj<J;sj++)
			{
				int k1=(si+J)*N+sj+J;
				int k2=si*J+sj;
				Atmp1[k2].x1*=BER[k1].x1;
				Atmp1[k2].x2*=BER[k1].x2;				
			}
		}

		array2_N_N *Atmp=(array2_N_N*)malloc(2*J*J*sizeof(array2_N_N));
		memset(Atmp,0,2*J*J*sizeof(array2_N_N));
		array2_N_N *EStmp=(array2_N_N*)malloc(2*J*J*sizeof(array2_N_N));
		memset(EStmp,0,2*J*J*sizeof(array2_N_N));		
		array2_N_N *ALtmp=(array2_N_N*)malloc(J*J*sizeof(array2_N_N));
		memset(ALtmp,0,J*J*sizeof(array2_N_N));

		for (int i1=0;i1<J;i1++)
		{
			for (int i2=0;i2<J;i2++)
			{
				int k1=i1*J+i2;
				int k2=2*k1;
				int k3=k2+1;
				int k4=(i1+J)*N+(i2+J);

				Atmp[k2].x1=Atmp[k2].x2=Atmp1[k1].x1;
				Atmp[k3].x1=Atmp[k3].x2=Atmp1[k1].x2;

				EStmp[k2].x1=ES[k4].x1;
				EStmp[k2].x2=ES[k4].x3;		
				EStmp[k3].x1=ES[k4].x2;
				EStmp[k3].x2=ES[k4].x4;
			}
		}
		//����EStmp*Atmp=Altmp
		for (si=0;si<J;si++)
		{
			for (int sj=0;sj<2*J;sj+=2)
			{
				int k1=si*2*J+sj;
				int k2=k1+1;
				int k3=k1/2;

 				ALtmp[k3].x1=Atmp[k1].x1*EStmp[k1].x1+Atmp[k2].x1*EStmp[k2].x1;
 				ALtmp[k3].x2=Atmp[k1].x2*EStmp[k1].x2+Atmp[k2].x2*EStmp[k2].x2;
			}
		}

		for (si=0;si<J;si++)
		{
			for (int sj=0;sj<J;sj++)
			{
				int k1=(si+J)*N+sj+J;
				int k2=si*J+sj;
				AL[k1]=ALtmp[k2];
			}
		}
		free(Atmp1);
		free(Atmp);
		free(EStmp);
		free(ALtmp);
	}
//��Ҷ˹���ʼ���
	for (k=2;k<=level;k++)
	{
		int J=int(pow(2,k-1));
		for (si=J;si<2*J;si++)
		{
			for (int sj=J;sj<2*J;sj++)
			{
				int k1=si*N+sj;
				double  temp1=AL[k1].x1*BE[k1].x1;
				double  temp2=AL[k1].x2*BE[k1].x2;
				P1[k1].x1=temp1/(temp1+temp2);
				P1[k1].x2=temp2/(temp1+temp2);
			}
		}
	}
	int m1=1*N+1;
	double temp1=AL[m1].x1*BE[m1].x1;
	double temp2=AL[m1].x2*BE[m1].x2;
	P1[m1].x1=temp1/(temp1+temp2);
	P1[m1].x2=temp2/(temp1+temp2);

	free(BE);
	free(BEP);
	free(BER);
	free(wtmp);
	free(scale);
	free(AL);
}

/*************************************************************************
**	�������ƣ�
**		 posthl()       
**	����
**		double			 *wavelete	С���ֽ�ϵ��ָ��	
**		array2_2_N_N	 *ES		ת�Ƹ���ģ��ϵ��ָ��
**		array2_N_N		 *SI		����ģ��ϵ��ָ��
**		array2_N_N		 *P1		�洢������ĺ������ָ��
**	����ֵ��
**		��
**	˵��:
**		���㴹ֱ������С��ϵ���ĺ�����ʣ���rice��ѧ�ṩ��MATLAB����ı�
************************************************************************/
void  posthl(double *wavelete,array2_2_N_N *ES,array2_N_N *SI,array2_N_N *P1,int N)
{
	array2_N_N *BE=(array2_N_N *)malloc(N*N*sizeof(array2_N_N));
	memset(BE,0,N*N*sizeof(array2_N_N));
	array2_N_N *BEP=(array2_N_N *)malloc(N*N*sizeof(array2_N_N));
	memset(BEP,0,N*N*sizeof(array2_N_N));
	array2_N_N *BER=(array2_N_N *)malloc(N*N*sizeof(array2_N_N));
	memset(BER,0,N*N*sizeof(array2_N_N));
	array2_N_N *wtmp=(array2_N_N *)malloc(N*N*sizeof(array2_N_N));
	memset(wtmp,0,N*N*sizeof(array2_N_N));
	array2_N_N *scale=(array2_N_N *)malloc(N*N*sizeof(array2_N_N));
	memset(scale,0,N*N*sizeof(array2_N_N));
	array2_N_N *AL=(array2_N_N *)malloc(N*N*sizeof(array2_N_N));
	memset(AL,0,N*N*sizeof(array2_N_N));

	//�����˹����
	for ( int n1=0;n1<N*N;n1++)
	{
		double temp1=1.0/sqrtf(2*PI*SI[n1].x1)*exp(-pow(wavelete[n1],2)/(2*SI[n1].x1))+0.0000000000001;
		wtmp[n1].x1=(double)temp1;
		double temp2=1.0/sqrtf(2*PI*SI[n1].x2)*exp(-pow(wavelete[n1],2)/(2*SI[n1].x2))+0.0000000000001;
		wtmp[n1].x2=(double)temp2;
		scale[n1].x1=scale[n1].x2=(double)((temp1+temp2)/2);
	}

	for (int si=0;si<(int)(pow(2,level-1));si++)
	{
		for (int sj=(int)(pow(2,level-1));sj<N;sj++)
		{
			int k=si*N+sj;
			BE[k].x1=wtmp[k].x1/scale[k].x1;
			BE[k].x2=wtmp[k].x2/scale[k].x2;
		}
	}

	//���ϲ���
	for (int k=level;k>1;k--)
	{
		int J=(int)(pow(2,k-1));
		int J2=J*J;
		
		double k1=ES[J*N+J].x1;
		double k2=ES[J*N+J].x2;		
		double k3=ES[J*N+J].x3;
		double k4=ES[J*N+J].x4;
	
		for (si=0;si<J;si++)
		{
			for (int sj=J;sj<2*J;sj++)
			{				
				int k=si*N+sj;					
				BEP[k].x1=k1*BE[k].x1+k3*BE[k].x2;
				BEP[k].x2=k2*BE[k].x1+k4*BE[k].x2;
			}
		}		
		//����beta���Ӿ����BE����
		array2_N_N *BCtmp=(array2_N_N*)malloc(J*J/4*sizeof(array2_N_N));
		memset(BCtmp,0,J*J/4*sizeof(array2_N_N));	
		for (si=0;si<J;si+=2)
		{
			for (int sj=J;sj<2*J;sj+=2)
			{				
				int k1=si*N+sj;	
				int k2=k1+1;	
				int k3=k1+N;	
				int k4=k3+1;
				
				int m1=si/2*N+sj/2;			

				double temp1=BEP[k1].x1*BEP[k2].x1*
							BEP[k3].x1*BEP[k4].x1;	
			 	double temp2=BEP[k1].x2*BEP[k2].x2*
 					    	BEP[k3].x2*BEP[k4].x2;	
				BE[m1].x1=temp1*wtmp[m1].x1/(scale[m1].x1*((temp1+temp2)/2));
				BE[m1].x2=temp2*wtmp[m1].x2/(scale[m1].x2*((temp1+temp2)/2));
				
				BER[k1].x1=BE[m1].x1/BEP[k1].x1;
				BER[k2].x1=BE[m1].x1/BEP[k2].x1;
				BER[k3].x1=BE[m1].x1/BEP[k3].x1;
				BER[k4].x1=BE[m1].x1/BEP[k4].x1;

				BER[k1].x2=BE[m1].x2/BEP[k1].x2;
				BER[k2].x2=BE[m1].x2/BEP[k2].x2;
				BER[k3].x2=BE[m1].x2/BEP[k3].x2;
				BER[k4].x2=BE[m1].x2/BEP[k4].x2;
			}
		}
 		free(BCtmp);
	}
	//DOWN_step
	//��AL����ֵ
	AL[1].x1=(double)22.627416997970;
	AL[1].x2=(double)1448.154687870049;

	for (k=2;k<=level;k++)
	{
		int J=int(pow(2,k-1));
		int J2=J*J;

		array2_N_N *Atmp1=(array2_N_N*)malloc(J*J*sizeof(array2_N_N));
		memset(Atmp1,0,J*J*sizeof(array2_N_N));

		for (int k1=0;k1<J/2;k1++)
		{
			for (int k2=J/2;k2<J;k2++)
			{
				int m1=k1*N+k2;
				int n1=2*(k1*J+k2-J/2);
				int n2=n1+1;
				int n3=n1+J;
				int n4=n3+1;
				Atmp1[n1].x1=Atmp1[n2].x1=Atmp1[n3].x1=Atmp1[n4].x1=AL[m1].x1;
				Atmp1[n1].x2=Atmp1[n2].x2=Atmp1[n3].x2=Atmp1[n4].x2=AL[m1].x2;
			}
		}
		for (si=0;si<J;si++)
		{
			for (int sj=0;sj<J;sj++)
			{
				int k1=si*N+sj+J;
				int k2=si*J+sj;
				Atmp1[k2].x1*=BER[k1].x1;
				Atmp1[k2].x2*=BER[k1].x2;				
			}
		}

		array2_N_N *Atmp=(array2_N_N*)malloc(2*J*J*sizeof(array2_N_N));
		memset(Atmp,0,2*J*J*sizeof(array2_N_N));	
		array2_N_N *EStmp=(array2_N_N*)malloc(2*J*J*sizeof(array2_N_N));
		memset(EStmp,0,2*J*J*sizeof(array2_N_N));		
		array2_N_N *ALtmp=(array2_N_N*)malloc(J*J*sizeof(array2_N_N));
		memset(ALtmp,0,J*J*sizeof(array2_N_N));

		for (int i1=0;i1<J;i1++)
		{
			for (int i2=0;i2<J;i2++)
			{
				int k1=i1*J+i2;
				int k2=2*k1;
				int k3=k2+1;
				int k4=i1*N+(i2+J);

				Atmp[k2].x1=Atmp[k2].x2=Atmp1[k1].x1;
				Atmp[k3].x1=Atmp[k3].x2=Atmp1[k1].x2;
			
				EStmp[k2].x1=ES[k4].x1;
				EStmp[k2].x2=ES[k4].x3;
				EStmp[k3].x1=ES[k4].x2;
				EStmp[k3].x2=ES[k4].x4;
			}
		}
		//����EStmp*Atmp=Altmp
		for (si=0;si<J;si++)
		{
			for (int sj=0;sj<2*J;sj+=2)
			{
				int k1=si*2*J+sj;
				int k2=k1+1;
				int k3=k1/2;

 				ALtmp[k3].x1=Atmp[k1].x1*EStmp[k1].x1+Atmp[k2].x1*EStmp[k2].x1;
 				ALtmp[k3].x2=Atmp[k1].x2*EStmp[k1].x2+Atmp[k2].x2*EStmp[k2].x2;
			}
		}
		for (si=0;si<J;si++)
		{
			for (int sj=0;sj<J;sj++)
			{
				int k1=si*N+sj+J;
				int k2=si*J+sj;
					AL[k1]=ALtmp[k2];
			}
		}
		free(Atmp1);
		free(Atmp);
		free(EStmp);
		free(ALtmp);
	}
	//��Ҷ˹���ʼ���
	for (k=2;k<=level;k++)
	{
		int J=int(pow(2,k-1));
		for (si=0;si<J;si++)
		{
			for (int sj=J;sj<2*J;sj++)
			{
				int k1=si*N+sj;
				double  temp1=AL[k1].x1*BE[k1].x1;
				double  temp2=AL[k1].x2*BE[k1].x2;
				P1[k1].x1=temp1/(temp1+temp2);
				P1[k1].x2=temp2/(temp1+temp2);
			}
		}
	}
	
	int m1=1;
	double temp1=AL[m1].x1*BE[m1].x1;
	double temp2=AL[m1].x2*BE[m1].x2;
	P1[m1].x1=temp1/(temp1+temp2);
	P1[m1].x2=temp2/(temp1+temp2);

	free(BE);
	free(BEP);
	free(BER);
	free(wtmp);
	free(scale);
	free(AL);
}

/*************************************************************************
**	�������ƣ�
**		 postlh()       
**	����
**		double			 *wavelete	С���ֽ�ϵ��ָ��	
**		array2_2_N_N	 *ES		ת�Ƹ���ģ��ϵ��ָ��
**		array2_N_N		 *SI		����ģ��ϵ��ָ��
**		array2_N_N		 *P1		�洢������ĺ������ָ��
**	����ֵ��
**		��
**	˵��:
**		����ˮƽ������С��ϵ���ĺ�����ʣ���rice��ѧ�ṩ��MATLAB����ı�
************************************************************************/
void  postlh(double *wavelete,array2_2_N_N *ES,array2_N_N *SI,array2_N_N *P1,int N)
{
	array2_N_N *BE=(array2_N_N *)malloc(N*N*sizeof(array2_N_N));
	memset(BE,0,N*N*sizeof(array2_N_N));
	array2_N_N *BEP=(array2_N_N *)malloc(N*N*sizeof(array2_N_N));
	memset(BEP,0,N*N*sizeof(array2_N_N));
	array2_N_N *BER=(array2_N_N *)malloc(N*N*sizeof(array2_N_N));
	memset(BER,0,N*N*sizeof(array2_N_N));
	array2_N_N *wtmp=(array2_N_N *)malloc(N*N*sizeof(array2_N_N));
	memset(wtmp,0,N*N*sizeof(array2_N_N));
	array2_N_N *scale=(array2_N_N *)malloc(N*N*sizeof(array2_N_N));
	memset(scale,0,N*N*sizeof(array2_N_N));
	array2_N_N *AL=(array2_N_N *)malloc(N*N*sizeof(array2_N_N));
	memset(AL,0,N*N*sizeof(array2_N_N));

	//�����˹����
	for ( int n1=0;n1<N*N;n1++)
	{
		double temp1=1.0/sqrtf(2*PI*SI[n1].x1)*exp(-pow(wavelete[n1],2)/(2*SI[n1].x1))+0.0000000000001;
		wtmp[n1].x1=(double)temp1;
		double temp2=1.0/sqrtf(2*PI*SI[n1].x2)*exp(-pow(wavelete[n1],2)/(2*SI[n1].x2))+0.0000000000001;
		wtmp[n1].x2=(double)temp2;
		scale[n1].x1=scale[n1].x2=(double)((temp1+temp2)/2);
	}
	for (int si=(int)(pow(2,level-1));si<N;si++)
	{
		for (int sj=0;sj<(int)(pow(2,level-1));sj++)
		{
			int k=si*N+sj;
			BE[k].x1=wtmp[k].x1/scale[k].x1;
			BE[k].x2=wtmp[k].x2/scale[k].x2;
		}
	}
	//���ϲ���
	for (int k=level;k>1;k--)
	{
		int J=(int)(pow(2,k-1));
		int J2=J*J;
		
		double k1=ES[J*N+J].x1;
		double k2=ES[J*N+J].x2;		
		double k3=ES[J*N+J].x3;
		double k4=ES[J*N+J].x4;
		
		for (si=J;si<2*J;si++)
		{
			for (int sj=0;sj<J;sj++)
			{				
				int k=si*N+sj;				
				BEP[k].x1=k1*BE[k].x1+k3*BE[k].x2;
				BEP[k].x2=k2*BE[k].x1+k4*BE[k].x2;
			}
		}
		//����beta���Ӿ����BE����
		array2_N_N *BCtmp=(array2_N_N*)malloc(J*J/4*sizeof(array2_N_N));
		memset(BCtmp,0,J*J/4*sizeof(array2_N_N));		
		//һ��ѭ����ֻ���±ꡣ
		for (si=J;si<2*J;si+=2)
		{
			for (int sj=0;sj<J;sj+=2)
			{				
				int k1=si*N+sj;	
				int k2=k1+1;	
				int k3=k1+N;	
				int k4=k3+1;
				
				int m1=si/2*N+sj/2;			

				double temp1=BEP[k1].x1*BEP[k2].x1*
							BEP[k3].x1*BEP[k4].x1;	
			 	double temp2=BEP[k1].x2*BEP[k2].x2*
 					    	BEP[k3].x2*BEP[k4].x2;	
				BE[m1].x1=temp1*wtmp[m1].x1/(scale[m1].x1*((temp1+temp2)/2));
				BE[m1].x2=temp2*wtmp[m1].x2/(scale[m1].x2*((temp1+temp2)/2));
			
				BER[k1].x1=BE[m1].x1/BEP[k1].x1;
				BER[k2].x1=BE[m1].x1/BEP[k2].x1;
				BER[k3].x1=BE[m1].x1/BEP[k3].x1;
				BER[k4].x1=BE[m1].x1/BEP[k4].x1;

				BER[k1].x2=BE[m1].x2/BEP[k1].x2;
				BER[k2].x2=BE[m1].x2/BEP[k2].x2;
				BER[k3].x2=BE[m1].x2/BEP[k3].x2;
				BER[k4].x2=BE[m1].x2/BEP[k4].x2;
			}
		}
 		free(BCtmp);
	}
	//DOWN_step
	//	��AL����ֵ
	AL[1*N].x1=(double)22.627416997970;
	AL[1*N].x2=(double)1448.154687870049;

	for (k=2;k<=level;k++)
	{
		int J=int(pow(2,k-1));
		int J2=J*J;

		array2_N_N *Atmp1=(array2_N_N*)malloc(J*J*sizeof(array2_N_N));
		memset(Atmp1,0,J*J*sizeof(array2_N_N));

		for (int k1=J/2;k1<J;k1++)
		{
			for (int k2=0;k2<J/2;k2++)
			{
				int m1=k1*N+k2;
				int n1=2*((k1-J/2)*J+k2);
				int n2=n1+1;
				int n3=n1+J;
				int n4=n3+1;
				Atmp1[n1].x1=Atmp1[n2].x1=Atmp1[n3].x1=Atmp1[n4].x1=AL[m1].x1;
				Atmp1[n1].x2=Atmp1[n2].x2=Atmp1[n3].x2=Atmp1[n4].x2=AL[m1].x2;
			}
		}		
		for (si=0;si<J;si++)
		{
			for (int sj=0;sj<J;sj++)
			{
				int k1=(si+J)*N+sj;
				int k2=si*J+sj;
				Atmp1[k2].x1*=BER[k1].x1;
				Atmp1[k2].x2*=BER[k1].x2;				
			}
		}

		array2_N_N *Atmp=(array2_N_N*)malloc(2*J*J*sizeof(array2_N_N));
		memset(Atmp,0,2*J*J*sizeof(array2_N_N));
	
		array2_N_N *EStmp=(array2_N_N*)malloc(2*J*J*sizeof(array2_N_N));
		memset(EStmp,0,2*J*J*sizeof(array2_N_N));
		
		array2_N_N *ALtmp=(array2_N_N*)malloc(J*J*sizeof(array2_N_N));
		memset(ALtmp,0,J*J*sizeof(array2_N_N));

		for (int i1=0;i1<J;i1++)
		{
			for (int i2=0;i2<J;i2++)
			{
				int k1=i1*J+i2;
				int k2=2*k1;
				int k3=k2+1;
				int k4=(i1+J)*N+i2;

				Atmp[k2].x1=Atmp[k2].x2=Atmp1[k1].x1;
				Atmp[k3].x1=Atmp[k3].x2=Atmp1[k1].x2;
				
				EStmp[k2].x1=ES[k4].x1;
				EStmp[k2].x2=ES[k4].x3;
				
				EStmp[k3].x1=ES[k4].x2;
				EStmp[k3].x2=ES[k4].x4;
			}
		}

		//����EStmp*Atmp=Altmp
		for (si=0;si<J;si++)
		{
			for (int sj=0;sj<2*J;sj+=2)
			{
				int k1=si*2*J+sj;
				int k2=k1+1;
				int k3=k1/2;

 				ALtmp[k3].x1=Atmp[k1].x1*EStmp[k1].x1+Atmp[k2].x1*EStmp[k2].x1;
 				ALtmp[k3].x2=Atmp[k1].x2*EStmp[k1].x2+Atmp[k2].x2*EStmp[k2].x2;
			}
		}
		for (si=0;si<J;si++)
		{
			for (int sj=0;sj<J;sj++)
			{
				int k1=(si+J)*N+sj;
				int k2=si*J+sj;
					AL[k1]=ALtmp[k2];
			}
		}
		free(Atmp1);
		free(Atmp);
		free(EStmp);
		free(ALtmp);
	}

	//��Ҷ˹���ʾ���
	for (k=2;k<=level;k++)
	{
		int J=int(pow(2,k-1));
		for (si=J;si<2*J;si++)
		{
			for (int sj=0;sj<J;sj++)
			{
				int k1=si*N+sj;
				double  temp1=AL[k1].x1*BE[k1].x1;
				double  temp2=AL[k1].x2*BE[k1].x2;
				P1[k1].x1=temp1/(temp1+temp2);
				P1[k1].x2=temp2/(temp1+temp2);
			}
		}
	}
		
	int m1=1*N;
	double temp1=AL[m1].x1*BE[m1].x1;
	double temp2=AL[m1].x2*BE[m1].x2;
	P1[m1].x1=temp1/(temp1+temp2);
	P1[m1].x2=temp2/(temp1+temp2);

	free(BE);
	free(BEP);
	free(BER);
	free(wtmp);
	free(scale);
	free(AL);
}

/*************************************************************************
**	�������ƣ�
**		 HMT_denoise()       
**	����
**		unsigned char *image		����ͼ���ָ��
**		int width                   ����ͼ��Ŀ��
**		int height					����ͼ��ĸ߶�
**		int wavelet_style           С����������
**		int layer					�ֽ���ع���С������
**	����ֵ��
**		��
**	˵��:
**		С����������Ʒ�ȥ�룬�Լ����ͼ�����С���ֽⲢ������HMT�ĺ�����ʣ�
**		�Եõ���ȥ��ģ�ͽ���С�����任���Ϳ��Եõ�ȥ���ͼ��
************************************************************************/
void HMT_denoise(unsigned char *image,int width,int height,int wavelet_style,int layer)
{
	int N= width;
	level = log10(N)/log10(2);
	if (level<layer)
		layer =level;
	double **si = new double* [2];
	for (int i=0;i<2;i++)
	{
		si[i]=new double [level];
	}
	for (i=0;i<2;i++)
	{
		for (int j=0;j<level;j++)
				si[i][j] =0;
	}
	double ***es=new double** [2];
	for (i=0;i<2;i++)
	{
		es[i]=new double*[2];
		for (int j=0;j<2;j++)
			es[i][j]=new double [level];
	}

	for (i=0;i<2;i++)
	{
		for (int j=0;j<2;j++)
			for (int k=0;k<level;k++)
				es[i][j][k] =0;
	}

//	double	si[2][level]={0};					//��ʼ��HMTģ�͵ķ���
//	double   es[2][2][level]={0};				//��ʼ��HMTģ�͵�ת�Ƹ��ʾ���
	hmtmodel(si,es);							//����HMTģ�͵Ĳ���

	//��ES��ʼ��ת�ƾ��󸳳�ֵ
	array2_2_N_N *ES=(array2_2_N_N *)malloc(N*N*sizeof(array2_2_N_N));
	memset(ES,0,N*N*sizeof(array2_2_N_N));	

	//�Գ�ʼ�ķ���ֲ�SI����ֵ
	array2_N_N *SI=(array2_N_N *)malloc(N*N*sizeof(array2_N_N));
	memset(SI,0,N*N*sizeof(array2_N_N));
	
	//����ES��ʼ��ת�ƾ���SI��ʼ�ķ���ֲ������ð󶨼�����
	for (int ii=level-1;ii>=0;ii--)
	{
		int k1=(int)(pow(2,ii+1));	
		for (int lrow=0;lrow<k1;lrow++)
		{
			for (int lcol=0;lcol<k1;lcol++)
			{
				int c1=lrow*N+lcol;		
				ES[c1].x1=es[0][0][ii];
				ES[c1].x2=es[0][1][ii];					
				ES[c1].x3=es[1][0][ii];			
				ES[c1].x4=es[1][1][ii];
				SI[c1].x1=si[0][ii];				
				SI[c1].x2=si[1][ii];
			}
		}
	}

	ES[0].x1=0;
	ES[0].x2=0;
	ES[0].x3=0;
	ES[0].x4=0;
	SI[0].x1=1;
	SI[0].x2=1;

	//�Ժ���ͼ�����С���ֽ⡣
	double *wavelet_coefficient=new double[N*N];     //����С��ϵ����ָ��
	memset(wavelet_coefficient,0,N*N*sizeof(double));

	//С�������ͼ���ϵ��
	//������HarrС��
	double HaarCoeffs [] = { 1.0/sqrt(2), 1.0/sqrt(2) };
	//DBС����Ŀǰ������DBС����Ч���ܺ�
	double Daub4Coeffs [] = { 0.4829629131445341,  0.8365163037378077,
		          0.2241438680420134, -0.1294095225512603 };
	double Daub6Coeffs [] = { 0.3326705529500825,  0.8068915093110924,
		0.4598775021184914, -0.1350110200102546,
		-0.0854412738820267,  0.0352262918857095 };

	double Coif1Coeffs[] = {
		-0.015656,	-0.072733,	0.384865,	0.852572,  0.337898,-0.072733  };

	//Symlets wavelet
	double Sym2Coeffs[] = {	-0.129410,	0.224144,	0.836516,	0.482963};

	double Sym3Coeffs[] = {	0.035226,	-0.085441,	-0.135011,	0.459878,	
			0.806892,	0.332671};


	double Daub20Coeffs [] = {
		-0.000013,	
			0.000094,	
			-0.000116,	
			-0.000686,	
			0.001992,	
			0.001395,	
			-0.010733,	
			0.003607,	
			0.033213,	
			-0.029458,	
			-0.071394,	
			0.093057,	
			0.127369,	
			-0.195946,	
			-0.249846,	
			0.281172,	
			0.688459,	
			0.527201,	
			0.188177,	
			0.026670	
};

	filterset choosedfilset;			//����һ��filterset�Ķ���

	//ѡ��С�������༰��ϵ��
	switch (wavelet_style)
	{
	case 0:
		{
			filterset filset(FALSE, HaarCoeffs,         2, 0);
			choosedfilset.copy(filset);
		}
		break;

	case 1:
		{
			filterset filset(FALSE, Daub4Coeffs,         4, 0);
			choosedfilset.copy(filset);
		}
		break;

	case 2:
		{
			filterset filset(FALSE, Daub6Coeffs,         6, 0);
			choosedfilset.copy(filset);
		}
		break;

	case 3:
		{
			filterset filset(FALSE, Coif1Coeffs,         6, 0);
			choosedfilset.copy(filset);
		}
		break;

	case 4:
		{
			filterset filset(FALSE, Sym2Coeffs,         4, 0);
			choosedfilset.copy(filset);
		}
		break;

	case 5:
		{
			filterset filset(FALSE, Sym3Coeffs,         6, 0);
			choosedfilset.copy(filset);
		}
		break;
	case 6:
		{
			filterset filset(FALSE, Daub20Coeffs,        20, 0);
			choosedfilset.copy(filset);
		}
		break;

	default:
		break; 
	}

	wavelett wlet(&choosedfilset);
	
//	int i;
	double *wavelet1=new double[N*N];
	memset(wavelet1,0,N*N*sizeof(double));
	for (i=0;i<N*N;i++)
	{
		wavelet1[i]=(double)((double)image[i]/255.0);	//��ͼ���һ��0~1
	}

	//����С���ֽ�
	wlet.transform2d(wavelet1,wavelet_coefficient,N, N,layer,-1);
	//�ֽ��С��ϵ��ת��һ�¡�
	double *temp=new double[N*N];
	memcpy(temp,wavelet_coefficient,N*N*sizeof(double));
	for(int i1=0;i1<N;i1++)
	{
		for (int i2=0;i2<N;i2++)
		{
			int k1=i1*N+i2;
			int k2=i2*N+i1;
			wavelet_coefficient[k1]=temp[k2];    //С��ϵ��ת��һ��
		}
	}
	
	//��ȡ��һ��С��ϵ������ֵ��
	//����һ�ηֽ�ĶԽǷ���С��ϵ������ȡ������������ֵ������ȡ��ֵ;
	double *median_image=new double[N*N/4];
	memset(median_image,0,N*N/4*sizeof(double));
	for (int p1=0;p1<N/2;p1++)
	{
		for (int p2=0;p2<N/2;p2++)
		{  
		 median_image[p1*N/2+p2]=
				      (double)fabs(wavelet_coefficient[(p1+N/2)*N+p2+N/2]);
		}
 	}
	//���Խ�ϵ����������
	qsort(median_image,N*N/4,sizeof(double),CompareGray);

	double median=(double)(median_image[N*N/8]);  //�õ���ֵ
	//�����˹��������
	double  sigma2=(double)pow(median/0.67,2);

	//��p1��������ʾ�������ռ�,p1��һ��2*N*N�����顣
	array2_N_N   *P1=(array2_N_N *)malloc(N*N*sizeof(array2_N_N));
	memset(P1,0,N*N*sizeof(array2_N_N));
	//SI=SI+sigma;
	for ( i=0;i<N*N;i++)
	{
		SI[i].x1+=sigma2;
		SI[i].x2+=sigma2;
	}

	posthh(wavelet_coefficient,ES,SI,P1,N);  //�ԽǷ����Ͻ���HMT������ʼ��㣬����ֵ����P1�С�

	postlh(wavelet_coefficient,ES,SI,P1,N);  //ˮƽ�����Ͻ���HMT������ʼ��㣬����ֵ����P1�С�

	posthl(wavelet_coefficient,ES,SI,P1,N);  //��ֱ�����Ͻ���HMT������ʼ��㣬����ֵ����P1�С�
	//SI=SI-sigma;
	for ( i=0;i<N;i++)
	{
		for (int j=0;j<N;j++)
		{
			int k1=i*N+j;
			SI[k1].x1=SI[k1].x1-sigma2;
			SI[k1].x2=SI[k1].x2-sigma2;
		}
	}
	
	//���ñ�Ҷ˹����ģ�ͼ�����Ʋ���������С��ϵ��
	double *wavelet_HMT=new double[N*N];
	memset(wavelet_HMT,0,N*N*sizeof(double));

	for (i=0;i<N*N;i++)
	{
		double temp1=SI[i].x1/(SI[i].x1+sigma2)*P1[i].x1*wavelet_coefficient[i];
		double temp2=SI[i].x2/(SI[i].x2+sigma2)*P1[i].x2*wavelet_coefficient[i];
		wavelet_HMT[i]=(double)(temp1+temp2);
	}

	//���任������ϽǸ�ֵ
	for (i=0;i<4;i++)
	{
		for (int j=0;j<4;j++)
		{
			int k=i*N+j;
			wavelet_HMT[k]=wavelet_coefficient[k];			
		}
	}
	memset(temp,0,N*N*sizeof(double));

	//��������ת��һ�µ���
	for (i=0;i<N;i++)
	{
		for (int j=0;j<N;j++)
		{
			int k1=i*N+j;
			int k2=j*N+i;
			temp[k1]=wavelet_HMT[k2];
		}
	}

	memcpy(wavelet_HMT,temp,N*N*sizeof(double));

	//����С�����任���������´�ص�image�С�
	memset(temp,0,N*N*sizeof(double));
	wlet.invert2d(wavelet_HMT,temp,N,N,layer,-1);

	for (i=0;i<N*N;i++)
	{
		int temp1=(int)(temp[i]*255+1);
		if (temp1<0)
			temp1=-temp1;
		if(temp1>255)
			temp1=255;
		image[i]=(BYTE)temp1;
	}
	
	free(ES);
	free(SI);
	free(P1);
	free(wavelet_coefficient);

	delete[]median_image;
	delete[]temp;
	delete[]wavelet1;
	delete [] wavelet_HMT;
	for (i=0;i<2;i++)
		delete [] si[i];
	delete [] si;
	for (i=0;i<2;i++)
	{
		for (int j=0;j<2;j++)
		{
			delete [] es[i][j];
		}
		delete [] es[i];
	}
	delete [] es;
}
/*************************************************************************
**	�������ƣ�
**		Double_Dimension_histEqual()
**	����
**		unsigned char* image  ָ��ͼ���ָ��
**		int width    ͼ��Ŀ��
**		int height   ͼ��ĸ߶�
**	����ֵ��
**		�� 
**	˵��:
**		˫��ֱ��ͼ��������ͼ�����ǿ������Աȶ�
**		��ֱ��ͼ���⴦��õ��ĵȼ����еȼ�����У��õ�����ͼ��
************************************************************************/
void Double_Dimension_histEqual(unsigned char* image,int width,int height)
{
	//1.����ͼ���ֱ��ͼ
	///////////////////////////////////////////////////////////////////
	int *m_histArray=new int [256];
	int i,x,y;
	
	//ֱ��ͼ������0
	for(i=0;i<256;i++)
		m_histArray[i]=0;
	
	//ͳ�ƻҶ�ֱ��ͼ
	for(y=0;y<height;y++)
	{
		for(x=0;x<width;x++)
		{
			int temp=image[width*y+x];
			m_histArray[temp]++;
		}
	}
	//���⻯����
	int sum=0;
	float Mape[256];            //ӳ���
	float mape_new[256];
	//����ӳ���
	for(i=0;i<256;i++)
	{
		sum+=m_histArray[i];
		Mape[i]=(float)(sum*255.0/(width*height)+0.5);  //ӳ������в�����ͬ��ֵ
	}
	
	//ͼ��ͶӰ��ӳ�����
	for (y=0;y<height;y++)
	{
		for (x=0;x<width;x++)
		{
			int temp=(int)Mape[image[width*y+x]];   //ӳ�����ȡֵ   
			image[width*y+x]=temp;
		}
	}
	//ֱ��ͼ���е�����ֵ��ʧ����ʣ�������ֵƽ�����зֲ���
	//����ͳ��ֱ��ͼ
	for(i=0;i<256;i++)
	{
		m_histArray[i]=0;
		Mape[i]=0;
		mape_new[i]=0;
	}
	
	for(y=0;y<height;y++)
	{
		for(x=0;x<width;x++)
		{
			int temp=image[width*y+x];
			m_histArray[temp]++;
		}
	}
	//ͳ�Ʒ���ĻҶȵȼ�������
	int num=1;
	for(i=0;i<256;i++)
	{
		if (m_histArray[i])
		{
			Mape[num-1]=(float)i;
			num++;
		}
	}
	
	//���½���ӳ���
	int count1=1;
	for(i=0;i<256;i++)
	{
		if (m_histArray[i])
		{
			mape_new[count1-1]=float(count1*255.0/(num-1));
			count1++;
		}
	}
	//�������ȼ����ȷֲ�
	for(y=0;y<height;y++)
	{
		for(x=0;x<width;x++)
		{
			for (i=0;i<num-1;i++)
			{
				if(image[width*y+x]==Mape[i])
				{
					image[width*y+x]=(int)mape_new[i];
				}
			}
		}
	}
	delete [] m_histArray;
}
/*************************************************************************
**	�������ƣ�
**		LateralInhibition()
**	����
**		unsigned char* image  ָ��ͼ���ָ��
**		int width    ͼ��Ŀ��
**		int height   ͼ��ĸ߶�
**      int Ksize    ģ��Ĵ�С
**	����ֵ��
**		�� 
**	˵��:
**		��������������ͼ����ǿ����ǿ�߿�ͻ������������
**		�������ġ�������������ͼ����ǿ�е�Ӧ�á������͵�
************************************************************************/
void LateralInhibition(unsigned char* image, int width, int height, int Ksize)
{
	int x,y,i,j;
	int LineBytes = ((width*8)+31)/32*4;  //һ����ռ���ֽ���
	int k2 = (Ksize-1)/2;
	unsigned char * NewImage= new unsigned char[height*LineBytes];
	memcpy(NewImage, image, LineBytes*height);

    unsigned char *MeanImage= new unsigned char[height*LineBytes];
	memcpy(MeanImage, image, LineBytes*height);//����ԭʼͼ�����ڴ��������3*3�����ڵľ�ֵ

	float *gg= new float[height*LineBytes];
    memset(gg,0,height*LineBytes*sizeof(float));
	float *EE = new float[height*LineBytes];
    memset(EE,0,height*LineBytes*sizeof(float));

	float *kernel=new float[Ksize*Ksize];

	//��ģ������,���õ��ǵ����˹�ֲ���ȷ������ϵ��
	for ( j=-k2; j<=k2; j++) 
	{
		for ( i=-k2; i<=k2; i++)
		{
			float dx=abs(j);
			float dy=abs(i);
			float sigma=10;
			float d=sqrt(dx*dx+dy*dy);
			float c=sqrt(2*PI);
			float beta=20;
			kernel[(j+k2)*Ksize+(i+k2)]=-1/(beta*c*sigma)*exp(-d*d/(2*sigma*sigma));
			kernel[0]=0;
			
		}
	}
   //�Ľ���3����ģ��
	
	for( y=k2; y<height-k2; y++) 
   {	
	   for( x=k2; x<width-k2; x++)
	   { 
		   //ˮƽģ��
		   float g1=0,g2=0;
		   for ( j=-k2; j<=k2; j++)
		   {  	  
			g1+=NewImage[(y+j)*LineBytes+x+1];
			g2+=NewImage[(y+j)*LineBytes+x-1];
		   }
		   g1=g1/3;
		   g2=g2/3;
		   
		   //��ֱģ��
		   float g3=0,g4=0;
		   for ( i=-k2; i<=k2; i++)
		   {  	  
			g3+=NewImage[(y+1)*LineBytes+x+i];
			g4+=NewImage[(y-1)*LineBytes+x-i];
		   }
		   g3=g3/3;
		   g4=g4/3;
		  

		   //��бģ��45
		   float g5=0,g6=0;
		   g5=(NewImage[(y-1)*LineBytes+x-1]+NewImage[(y-1)*LineBytes+x]+NewImage[y*LineBytes+x-1])/3;//���Ͻ������ؾ�ֵ
		   g6=(NewImage[(y+1)*LineBytes+x]+NewImage[(y+1)*LineBytes+x+1]+NewImage[y*LineBytes+x+1])/3;//���½������ؾ�ֵ
		
           //��бģ��
		   float g7=0,g8=0;
		   g7=(NewImage[y*LineBytes+x-1]+NewImage[(y+1)*LineBytes+x-1]+NewImage[(y+1)*LineBytes+x])/3;//���Ͻ������ؾ�ֵ
		   g8=(NewImage[(y-1)*LineBytes+x]+NewImage[(y-1)*LineBytes+x+1]+NewImage[y*LineBytes+x+1])/3;//���Ͻ������ؾ�ֵ
		  
		   
		   //�ҳ��ĸ���������ֵΪgmax
		   unsigned char g[4];
		  
           g[0]=abs(g1-g2);
		   g[1]=abs(g3-g4);
		   g[2]=abs(g5-g6);
		   g[3]=abs(g7-g8);
		   for (int k =0;k<4;k++)
		   {
			   if(g[k]>gg[y*LineBytes+x])
				    gg[y*LineBytes+x]=g[k];
		   }

	   }
   }
   //�ҳ�gg[]������һ��
   unsigned char max;
   max=0;
   for( y=k2; y<height-k2; y++) 
   {
	   for( x=k2; x<width-k2; x++)
	   {
		   if(gg[y*LineBytes+x]>max)
			   max=gg[y*LineBytes+x];
	   }
   }

   //����E[y*LineBytes+x]
   for( y=k2; y<height-k2; y++) 
   {
	   for( x=k2; x<width-k2; x++)
	   {
		   EE[y*LineBytes+x]= gg[y*LineBytes+x]/max;
	   }
   }

   //���빫ʽ
   for( y=k2; y<height-k2; y++) 
   {
	   for( x=k2; x<width-k2; x++)
	   {
		   float sum = 0;
		   float sum1 = 0;
		   
		   for ( j=-k2; j<=k2; j++) 
		   { 
			   for ( i=-k2; i<=k2; i++)
			   {
				   sum += MeanImage[(y+j)*LineBytes+x+i];
				   sum1 += NewImage[(y+j)*LineBytes+x+i] * kernel[(j+k2)*Ksize+(i+k2)];

			   }
		   }
          MeanImage[y*LineBytes+x]=(int)sum/(Ksize*Ksize);
		  image[y*LineBytes+x]=MeanImage[y*LineBytes+x]+sum1+EE[y*LineBytes+x]*(MeanImage[y*LineBytes+x]-NewImage[y*LineBytes+x]);
	   }
   }

   //����߽� 
   for( y=0; y<=height-1; y++) 
   {
	   for( x=0; x<=width-1; x++)
	   {
		   //�����ϱ߽�
		   if(y==0 || y==1)
		   {
			   for(int x1=0;x1<width-1;x1++)
			   {
				 image[y*LineBytes+x1]=image[2*LineBytes+x1];
			   }
		   }
		   //�����±߽�
		   if(y==height-1 || y== height-2)
		   {
			   for(int x1=0;x1<width-1;x1++)
			   {
				 image[y*LineBytes+x1]=image[(height-3)*LineBytes+x1];
			   }
		   }
		   //������߽�
		  if(x==0 || x==1)
		   {
			   for(int y1=1;y1<height-2;y1++)
			   {
				 image[y1*LineBytes+x]=image[y1*LineBytes+2];
			   }
		   }
		  //�����ұ߽�
		  if(x==width-1 || x==width-2)
		   {
			   for(int y1=1;y1<height-2;y1++)
			   {
				 image[y1*LineBytes+x]=image[y1*LineBytes+width-3];
			   }
		   }
	   }
   }

   int Min=255;
   int Max=1;

   for(y=0; y<=height-1; y++)
   {
	   for(x=0; x<=width-1; x++)
	   {
		   if(image[y*LineBytes+x]>Max)
			   Max = image[y*LineBytes+x];
		   if(image[y*LineBytes+x]<Min)
			   Min = image[y*LineBytes+x];
	   }
   }

   for(y=0; y<=height-1; y++)
   {
	   for(x=0; x<=width-1; x++)
	   {
		   float s=float(image[y*LineBytes+x]-Min);
		   float u=float(Max-Min);
		   image[y*LineBytes+x] = int(s*255/u);
	   }
   }
   delete [] gg ;
   delete [] EE;
   delete [] NewImage;
   delete [] MeanImage;
}




/**************************************************************************
*                           �ֲ�����ʵ��                                  *
**************************************************************************/

void  Image_Reserve(unsigned char *image,int width,int height)
{
	int x,y;
	long lLBytes = (width+3 )/4*4;
	unsigned char * temp = new unsigned char [lLBytes*height];
	memset(temp,0,sizeof(unsigned char )*lLBytes*height);
	for (y = 0;y<height;y++)
	{
		for (x = 0;x<width;x++)
		{
			int k=lLBytes*y+x;
			temp[k]=255-image[k];
		}
	}
	memcpy(image,temp,sizeof(unsigned char )*lLBytes*height);
	delete []temp;	
	temp = NULL;
}

void Roberts(unsigned char *image,int width,int height)
{
	//ÿ��������ռ�ֽ��������ͼ��������ͼ����ͬ
	int lineByte=(width+3)/4*4;
	//�������ͼ�񻺳���
	unsigned char *temp=new unsigned char[lineByte*height];
	memset(temp,0,lineByte*height*sizeof(unsigned char));
	//ѭ��������ͼ�������
	int i,j;
	//�м����
	int x, y, t;
	//Roberts����
	for(i=1;i<height-1;i++)
	{
		for(j=1;j<width-1;j++)
		{
			//x�����ݶ�
			x=*(image+i*lineByte+j)
				-*(image+(i+1)*lineByte+j);
			//y�����ݶ�
			y=*(image+i*lineByte+j)
				-*(image+i*lineByte+(j+1));
			t=(int)(sqrt(x*x+y*y)+0.5);
			if(t>255)
				t=255;
			*(temp+i*lineByte+j)=t;
		}
	}
	memcpy(image,temp,lineByte*height*sizeof(unsigned char));
	delete[]temp;
	temp = NULL;
}


void Sobel(unsigned char *image,int width,int height)
{
	//ÿ��������ռ�ֽ��������ͼ��������ͼ����ͬ
	int lineByte=(width+3)/4*4;
	//�������ͼ�񻺳���
	unsigned char *temp=new unsigned char[lineByte*height];
	memset(temp,0,lineByte*height*sizeof(unsigned char));
	//ѭ��������ͼ�������
	int i,j;
	//�м����
	int x, y, t;
	//Sobel����
	for(i=1;i<height-1;i++)
	{
		for(j=1;j<width-1;j++)
		{
			//x�����ݶ�
			x= *(image+(i-1)*lineByte+(j+1))
				+ 2 * *(image+i*lineByte+(j+1))
				+ *(image+(i+1)*lineByte+(j+1))
				- *(image+(i-1)*lineByte+(j-1))
				- 2 * *(image+i*lineByte+(j-1))
				- *(image+(i+1)*lineByte+(j-1));
			
			//y�����ݶ�
			y= *(image+(i-1)*lineByte+(j-1))
				+ 2 * *(image+(i-1)*lineByte+j)
				+ *(image+(i-1)*lineByte+(j+1))
				- *(image+(i+1)*lineByte+(j-1))
				- 2 * *(image+(i+1)*lineByte+j)
				- *(image+(i+1)*lineByte+(j+1));
			
// 			xx=*(image+(i-1)*lineByte+(j-1))
// 				+ 2 * *(image+(i-1)*lineByte+j)
// 				+ *(image+(i-1)*lineByte+(j+1))
// 				- *(image+(i+1)*lineByte+(j-1))
// 				- 2 * *(image+(i+1)*lineByte+j)
// 				- *(image+(i+1)*lineByte+(j+1));
// 		
// 			yy=*(image+(i-1)*lineByte+(j-1))
// 				+ 2 * *(image+(i-1)*lineByte+j)
// 				+ *(image+(i-1)*lineByte+(j+1))
// 				- *(image+(i+1)*lineByte+(j-1))
// 				- 2 * *(image+(i+1)*lineByte+j)
// 				- *(image+(i+1)*lineByte+(j+1));

			t=(int)(sqrt(x*x+y*y)+0.5);
			if(t>255)
				t=255;
			*(temp+i*lineByte+j)=t;
			}
		}
	memcpy(image,temp,lineByte*height*sizeof(unsigned char));
	delete[]temp;


	
}	

void Prewitt(unsigned char *image,int width,int height)
{
	//ÿ��������ռ�ֽ��������ͼ��������ͼ����ͬ
	int lineByte=(width+3)/4*4;
	//�������ͼ�񻺳���
	unsigned char *temp=new unsigned char[lineByte*height];
	memset(temp,0,lineByte*height*sizeof(unsigned char));
	//ѭ��������ͼ�������
	int i,j;
	//�м����
	int x, y, t;
	//Prewitt����
	for(i=1;i<height-1;i++)
	{
		for(j=1;j<width-1;j++)
		{
			//x�����ݶ�
			x= *(image+(i-1)*lineByte+(j+1))
				+ *(image+i*lineByte+(j+1))
				+ *(image+(i+1)*lineByte+(j+1))
				- *(image+(i-1)*lineByte+(j-1))
				- *(image+i*lineByte+(j-1))
				- *(image+(i+1)*lineByte+(j-1));
			
			//y�����ݶ�
			y= *(image+(i-1)*lineByte+(j-1))
				+ *(image+(i-1)*lineByte+j)
				+ *(image+(i-1)*lineByte+(j+1))
				- *(image+(i+1)*lineByte+(j-1))
				- *(image+(i+1)*lineByte+j)
				- *(image+(i+1)*lineByte+(j+1));
			
			t=(int)(sqrt(x*x+y*y)+0.5);
			if(t>255)
				t=255;
			*(temp+i*lineByte+j)=t;
		}
	}
	memcpy(image,temp,lineByte*height*sizeof(unsigned char));
	delete[]temp;
}

void Laplacian(unsigned char *image,int width,int height)
{
	//ÿ��������ռ�ֽ��������ͼ��������ͼ����ͬ
	int lineByte=(width+3)/4*4;
	//�������ͼ�񻺳���
	unsigned char *temp=new unsigned char[lineByte*height];
	memset(temp,0,lineByte*height*sizeof(unsigned char));
	//ѭ��������ͼ�������
	int i,j;
	//�м����
	int t;
	//Laplacian����
	for(i=1;i<height-1;i++)
	{
		for(j=1;j<width-1;j++)
		{
			t= 4 * *(image+i*lineByte+j)
				- *(image+(i-1)*lineByte+j)
				- *(image+(i+1)*lineByte+j)
				- *(image+i*lineByte+(j-1))
				- *(image+i*lineByte+(j+1));
			t=(int)(abs(t)+0.5);
			if(t>255)
				t=255;
			*(temp+i*lineByte+j)=t;
		}
	}
	memcpy(image,temp,lineByte*height*sizeof(unsigned char));
	delete[]temp;
}


////////////////////////////////////////////////////////////////////////////////
double* GaussKernel( long Ksize=3, double sigma = 1.0 )
{
	long  i,j,k2;
	k2    = (Ksize-1)/2;
	
	double *kernel= (double*)malloc(Ksize*Ksize*sizeof(double));
	
	/* Create Filter */
	for (i=0; i<Ksize; i++) {
		for (j=0;j < Ksize;j++) {
			long   tmp1 = (i-k2)*(i-k2) + (j-k2)*(j-k2);
			double tmp2 = (-1*tmp1) / (2*sigma*sigma);
			kernel[i*Ksize+j] = exp(tmp2);
		}
	}
	/* Normalize Filter */
	double sum = 0;
	for (i=0;i < Ksize;i++)
		for (j=0;j < Ksize;j++)
			sum = sum + kernel[i*Ksize+j];
		for (i=0;i < Ksize;i++)
			for (j=0;j < Ksize;j++)
				kernel[i*Ksize+j] = kernel[i*Ksize+j] / sum;	
			
			return kernel;
			//free(kernel);
}

bool Conv2d(LPSTR lpBits, long lWidth, long lHeight, double *kernel, 
			long Ksize = 3, double fcoef=1)
{
	HGLOBAL         hMemy;
	LPSTR           lpNewBits;
	
	long    i,j,x,y;
	double	fsum;
	long    k2 = (Ksize-1)/2;
	long    lLBytes = ((lWidth*8)+31)/32*4;
	
	if(!(hMemy = ::GlobalAlloc(GHND,lLBytes*lHeight)))
		return FALSE;
	
	lpNewBits = (LPSTR)::GlobalLock(hMemy);
	memcpy(lpNewBits, lpBits, lLBytes*lHeight);
	
	for( y=k2; y<lHeight-k2; y++) {
		for( x=k2; x<lWidth-k2; x++) {			
			fsum = 0;
			for ( j=-k2; j<=k2; j++) {
				for ( i=-k2; i<=k2; i++) {
					fsum += (unsigned char)lpNewBits[(y+j)*lLBytes+x+i] * kernel[(j+k2)*Ksize+i+k2];
				}
			}
			fsum *= fcoef;
			if(fsum > 255) 
				lpBits[y*lLBytes+x] = (unsigned char)255;
			else
				lpBits[y*lLBytes+x] = (unsigned char)(fsum+0.5);
		}
	}
	::GlobalUnlock(hMemy);	
	::GlobalFree(hMemy);
	return TRUE;
}
/*******************************************************
  Program Steps:
  1) Gaussian kernel for noise removal
  2) Use Sobel masks to approximate gradient
  3) Find orientation of gradient
  4) Implement nonmaximum suppression to assign edges
  5) Hysteresis thresholding (2 thresholds)
********************************************************/
#define  M_PI   3.141592653
BOOL Canny(LPSTR lpBits,long lWidth, long lHeight, long lThresh = 40, long hThresh = 120 )
{
	long    	x, y, i, j;
	float		sumX, sumY, SUM;
	int 		leftPixel, rightPixel;
	int 		P1, P2, P3, P4, P5, P6, P7, P8;
	float		ORIENT;
	int			edgeDirection;	
	HGLOBAL     hMemG;
	LPSTR       imgG;
	
	long lLBytes = ((lWidth*8)+31)/32*4;
	if(!(hMemG = ::GlobalAlloc(GHND,lLBytes*lHeight))) 
		return FALSE;
	imgG = (LPSTR)::GlobalLock(hMemG);

	/* 3x3 GX Sobel mask. */
	double dx[] = {-1, 0, 1,-2, 0, 2,-1, 0, 1};
	
	/* 3x3 GY Sobel mask. */
	double dy[] = { 1, 2, 1, 0, 0, 0,-1,-2,-1};
	
	/***********************************************
	*	3X3 GAUSSIAN MASK ALGORITHM STARTS HERE    *
	***********************************************/

	memcpy(imgG, lpBits, lLBytes*lHeight);
    Conv2d(imgG, lWidth, lHeight, GaussKernel(3,1.0), 3);

	/***********************************************
	*	SOBEL GRADIENT APPROXIMATION
	***********************************************/
	for( y=1; y < lHeight-1; y++)  {
		for( x=1; x<lWidth-1; x++)  {
			
			sumX = 0;
			sumY = 0;		   
			/***********************************
			* X gradient approximation
			***********************************/
			for(j=-1; j<=1; j++) {
				for(i=-1; i<=1; i++) {
					sumX = sumX + (long)((BYTE)imgG[(y+j)*lLBytes+x+i] * dx[(j+1)*3+i+1]);
				}
			}
			/**************************
			* Y gradient approximation
			**************************/
			for(i=-1; i<=1; i++)  {
				for(j=-1; j<=1; j++)  {
					sumY = sumY + (long)((BYTE)imgG[(y+j)*lLBytes+x+i] * dy[(j+1)*3+i+1]);
				}
			}
			/***********************************************
			* GRADIENT MAGNITUDE APPROXIMATION 
			***********************************************/
			SUM =(float) (fabs(sumX) + fabs(sumY));
			if(SUM>255)    SUM=255;
			if(SUM<0)      SUM=0;

			/***************************
			* Magnitude orientation
			***************************/
			/* Cannot divide by zero*/
			if(sumX == 0) {
				if(sumY==0) ORIENT = 0.0;
				else if (sumY<0) {
					sumY = -sumY;
					ORIENT = 90.0;
				}
				else ORIENT = 90.0;
			}
			/* Can't take invtan of angle in 2nd Quad */
			else if(sumX<0 && sumY>0)
            {
				sumX = -sumX;
				ORIENT = 180 - (float)(((atan((double)(sumY)/(double)(sumX))) * (180/M_PI)));
			}
			/* Can't take invtan of angle in 4th Quad */
			else if(sumX>0 && sumY<0) {
				sumY = -sumY;
				ORIENT = 180 -(float)(((atan((float)(sumY)/(float)(sumX))) * (180/M_PI)));
			}
			/* else angle is in 1st or 3rd Quad */
			else
				ORIENT =(float)((atan((float)(sumY)/(float)(sumX))) * (180/M_PI));
				
			/***************************************************
			* Find edgeDirection by assigning ORIENT a value of
			* either 0, 45, 90 or 135 degrees, depending on which
			* value ORIENT is closest to
			****************************************************/
				
			if(ORIENT < 22.5) edgeDirection = 0;
			else if(ORIENT < 67.5) edgeDirection = 45;
			else if(ORIENT < 112.5) edgeDirection = 90;
			else if(ORIENT < 157.5) edgeDirection = 135;
			else edgeDirection = 0;
			
			/***************************************************
			* Obtain values of 2 adjacent pixels in edge
			* direction.
			****************************************************/
			if(edgeDirection == 0) {
				leftPixel = (BYTE)imgG[y*lLBytes+x-1];
				rightPixel = (BYTE)imgG[y*lLBytes+x+1];
			}
			else if(edgeDirection == 45) {
				leftPixel = (BYTE)imgG[(y+1)*lLBytes+x-1];
				rightPixel = (BYTE)imgG[(y-1)*lLBytes+x+1];
			}
			else if(edgeDirection == 90) {
				leftPixel = (BYTE)imgG[(y-1)*lLBytes+x];
				rightPixel = (BYTE)imgG[(y+1)*lLBytes+x];
			}
			else  {
				leftPixel = (BYTE)imgG[(y-1)*lLBytes+x-1];
				rightPixel = (BYTE)imgG[(y+1)*lLBytes+x+1];
			}

			/*********************************************
			* Compare current magnitude to both adjacent
			* pixel values.  And if it is less than either
			* of the 2 adjacent values - suppress it and make
			* a nonedge.
			*********************************************/
			if(SUM < leftPixel || SUM < rightPixel)
				SUM = 0;
			else {
			/**********************
			* Hysteresis
			**********************/
				if(SUM >= hThresh)
					SUM = 255; /* edge */
				else if(SUM < lThresh)
					SUM = 0;  /* nonedge */
				/* SUM is between T1 & T2 */
				else  {
					/* Determine values of neighboring pixels */
					P1 = (int)imgG[(y-1)*lLBytes+x-1];
					P2 = (int)imgG[(y-1)*lLBytes+x];
					P3 = (int)imgG[(y-1)*lLBytes+x+1];
					P4 = (int)imgG[y*lLBytes+x-1];
					P5 = (int)imgG[y*lLBytes+x+1];
					P6 = (int)imgG[(y+1)*lLBytes+x-1];
					P7 = (int)imgG[(y+1)*lLBytes+x];
					P8 = (int)imgG[(y+1)*lLBytes+x+1];
                 		
					/* Check to see if neighboring pixel values are edges */
					if(P1 > hThresh || P2 > hThresh || P3 > hThresh || P4 > hThresh 
						            || P5 > hThresh || P6 > hThresh || P7 > hThresh 
								    || P8 > hThresh)
						SUM = 255; /* make edge */
					else 
						SUM = 0; /* make nonedge */
				}
			}
			*(lpBits + x + y*lLBytes) = (unsigned char)(SUM);
	   }
    }
   	::GlobalUnlock(hMemG);	
	::GlobalFree(hMemG);
	return true;
}

void Edge_RotateInvariantOperator(BYTE *InputImage, int ImgW, int ImgH, int Masksize)
{
	int	 Imagesize=ImgW*ImgH;
	int  i,j;
	int  Area = Masksize*Masksize; 
	int  HalfArea = Area/2;
	int	 r=Masksize/2; 
	BYTE *tmpImage = new BYTE[Imagesize];
    BYTE **RowAddress = new BYTE*[ImgH];
	//initialize 
	double *Operatorx = new double[Area];
	double *Operatory = new double[Area];
	memset(tmpImage, 0, Imagesize);
	memset(Operatorx, 0, Area*sizeof(double));
	memset(Operatory, 0, Area*sizeof(double));
	
	double *tmpOpx = Operatorx , *tmpOpy = Operatory;
	double Coef1 = PI/double(r+1), Mr;
	int    dx, dy, R02=r*r, Radius, nG;
	for(dy=-r; dy<=r; dy++)
	{
		for(dx=-r; dx<=r; dx++)
		{
			Radius = dx*dx +dy*dy;
			if(Radius<=R02)
			{
				Mr = sqrt(Radius);
				*tmpOpx = cos(2*Mr*Coef1)*sin(dx*Coef1);
				*tmpOpy = cos(2*Mr*Coef1)*sin(dy*Coef1);
			}
			tmpOpx++;
			tmpOpy++;
		}
	}
	
	RowAddress[0] = InputImage;
	for(i=1 ; i<ImgH ;i++)
	{
		RowAddress[i] = RowAddress[i-1] + ImgW;
	}
	
	
	double gx, gy;
	BYTE *Image = tmpImage + r*ImgW;
	for(j=r; j<ImgH-r; j++)
	{		
		for(i=r; i<ImgW-r; i++)
		{
			gx = gy = 0;
			tmpOpx = Operatorx;
			tmpOpy = Operatory;
			for(dy=-r; dy<=r; dy++)
			{
				for(dx=-r; dx<=r; dx++)
				{
					gx += (*tmpOpx++)* RowAddress[j+dy][i+dx];
					gy += (*tmpOpy++)* RowAddress[j+dy][i+dx];
				}
			}
			nG  = int(sqrt(gx*gx+gy*gy) + 0.5);
			Image[i] = (nG>255) ? 255: nG;		
		}
		Image += ImgW;
	}
	memcpy(InputImage,tmpImage,Imagesize);
	
	if(tmpImage!=NULL)  delete []tmpImage;
	if(RowAddress!=NULL)  delete []RowAddress;
	if(Operatorx!=NULL)  delete []Operatorx;
	if(Operatory!=NULL)  delete []Operatory;
}

void Gradient_Sharper(unsigned char *image,int width,int height,int thresh=30)
{
	int lineByte = (width + 3) / 4 * 4;
	
	unsigned char*	pSrc;       // ָ��Դͼ���ָ��
	unsigned char*	pDst; 	
	unsigned char*	pSrc1;
	unsigned char*	pSrc2;	

	unsigned char*  m_temp=new unsigned char[lineByte*height]; 
	memset(m_temp,0,lineByte*height*sizeof(unsigned char));

	LONG	i,j;				// ѭ������
	int	bTemp;

	for(i = 0; i < height; i++)		// ÿ��
	{		
		for(j = 0; j < width; j++)		// ÿ��
		{
			//ָ����DIB��i�е�j�е����ص�ָ��
			pDst = m_temp + lineByte * (height -1 - i) + j;
			// �����ݶ�����
			// ָ��DIB��i�У���j�����ص�ָ��
			pSrc  = image+ lineByte * (height - 1 - i) + j;			
			// ָ��DIB��i+1�У���j�����ص�ָ��
			pSrc1 =  image + lineByte * (height - 2 - i) + j;			
			// ָ��DIB��i�У���j+1�����ص�ָ��
			pSrc2 =  image + lineByte * (height - 1 - i) + j + 1;
			
			bTemp = abs((*pSrc)-(*pSrc1)) + abs((*pSrc)-(*pSrc2));
			
			// �ж��Ƿ�С����ֵ
			if ((bTemp+120) < 255)
			{
				// �ж��Ƿ������ֵ������С��������Ҷ�ֵ���䡣
				if (bTemp >= thresh)
				{
					*pSrc = bTemp + 120;
				}
			}
			else
			*pSrc = 255;
			//�����µ�DIB����ֵ
			*pDst = *pSrc;
		}
	}
	memcpy(image,m_temp,lineByte*height*sizeof(unsigned char));
	delete[]m_temp;
}

void  Roberts_Sharper(BYTE *InputImage,int ImgW,int ImgH)
{
	int	 Imagesize = ImgW*ImgH;
	if (InputImage==NULL || Imagesize==0) return ;
	
	int i,j;
	BYTE *Image ;
	BYTE *tmpImage = new BYTE[Imagesize];
    BYTE **RowAddress = new BYTE*[ImgH];
	
	//initialize 
	memset(tmpImage, 0, Imagesize);
	RowAddress[0] = InputImage;
	for(i=1 ; i<ImgH ;i++)
	{
		RowAddress[i] = RowAddress[i-1] + ImgW;
	}
	
	int dx,dy;
	Image = tmpImage + ImgW;
	for(j=0; j<ImgH-1; j++)
	{		
		for(i=0; i<ImgW-1; i++)
		{
			dx = abs(RowAddress[ j ][ i ] - RowAddress[j+1][i+1]);
			dy = abs(RowAddress[j+1][ i ] - RowAddress[ j ][i+1]);
			Image[i] = (dx>dy) ? dx: dy;				 
		}
		Image += ImgW;
	}
	
	memcpy(InputImage, tmpImage, Imagesize);
	if(  tmpImage!=NULL)  delete []tmpImage;
	if(RowAddress!=NULL)  delete RowAddress;
}

void Laplacian_sharper(BYTE *InputImage,int ImgW,int ImgH)
{
	int	 Imagesize = ImgW*ImgH;
	if (InputImage==NULL || Imagesize==0) return ;
	
	int i,j;
	BYTE *Image ;
	BYTE *tmpImage = new BYTE[Imagesize];
    BYTE **RowAddress = new BYTE*[ImgH];
	
	//initialize 
	memset(tmpImage, 0, Imagesize);
	RowAddress[0] = InputImage;
	for(i=1 ; i<ImgH ;i++)
	{
		RowAddress[i] = RowAddress[i-1] + ImgW;
	}
	
	Image = tmpImage + ImgW;
	int nG;
	for(j=1; j<ImgH-1; j++)
	{		
		for(i=1; i<ImgW-1; i++)
		{
			nG  = RowAddress[j-1][ i ] + RowAddress[j+1][ i ];
			nG += RowAddress[ j ][i-1] + RowAddress[ j ][i+1];
			nG  = abs(5*RowAddress[j][i] - nG);
			Image[i] = (nG>255) ? 255 : nG;				 
		}
		Image += ImgW;
	}
	
	memcpy(InputImage,tmpImage,ImgW*ImgH);
	if(  tmpImage!=NULL)  delete []tmpImage;
	if(RowAddress!=NULL)  delete []RowAddress;
}


void int_to_uchar(int *r,BYTE *in,int size)
{
	int i,
		max_r=r[0],
		min_r=r[0];
	double fTemp;
	for (i=1; i<size; i++)
	{
		if ( r[i] > max_r )
			max_r=r[i];
		else if ( r[i] < min_r )
			min_r=r[i];
	}
	
	max_r-=min_r;
	if(max_r!=0)
	{
		for (i=0; i<size; i++)
		{
			fTemp=(r[i]-min_r)*255.0/max_r+0.5; 
			in[i] = (BYTE)(fTemp);
		}
	}
}

void Image_SUSAN_Principle(BYTE* InputImage , int ImgW , int ImgH ,int bt)
{
	int	 Imagesize = ImgW*ImgH;
	if (InputImage==NULL || Imagesize==0) return ;
	
	int    i, j, n, max_no=1850;
	int    Width=ImgW, Height=ImgH;
    int    *r  =new int[Imagesize];
	BYTE   *mid = new BYTE[Imagesize];
	BYTE   *in=InputImage;
	BYTE   *p,*cp,*bp,bbp[516];
    double temp; 
     
    bbp[0]=bbp[1]=bbp[515]=bbp[514]=0;
    bp=bbp+258;
   
    for(int k=-256;k<257;k++)
	{
        temp=((float)k)/((float)bt);
        temp=temp*temp;
		temp=temp*temp*temp;
        temp=100.0*exp(-temp)+0.5;
        bbp[k+258] = (BYTE)(temp);
	}
    memset (mid,100,Width*Height);
    memset (r,0,Width*Height * sizeof(int));

    for(i=3;i<Height-3;i++)
       for(j=3;j<Width-3;j++)
	   {
		   n=100;
           p=in + (i-3)*Width + j - 1;
           cp=bp + in[i*Width+j];

           n+=*(cp-*p++);
           n+=*(cp-*p++);
           n+=*(cp-*p);
           p+=Width-3; 

           n+=*(cp-*p++);
           n+=*(cp-*p++);
           n+=*(cp-*p++);
           n+=*(cp-*p++);
           n+=*(cp-*p);
           p+=Width-5;

           n+=*(cp-*p++);
           n+=*(cp-*p++);
           n+=*(cp-*p++);
           n+=*(cp-*p++);
           n+=*(cp-*p++);
           n+=*(cp-*p++);
           n+=*(cp-*p);
           p+=Width-6;

           n+=*(cp-*p++);
           n+=*(cp-*p++);
           n+=*(cp-*p);
           p+=2;
           n+=*(cp-*p++);
           n+=*(cp-*p++);
           n+=*(cp-*p);
           p+=Width-6;

           n+=*(cp-*p++);
           n+=*(cp-*p++);
           n+=*(cp-*p++);
           n+=*(cp-*p++);
		   n+=*(cp-*p++);
	       n+=*(cp-*p++);
		   n+=*(cp-*p);
		   p+=Width-5;

		   n+=*(cp-*p++);
		   n+=*(cp-*p++);
           n+=*(cp-*p++);
           n+=*(cp-*p++);
           n+=*(cp-*p);
           p+=Width-3;

           n+=*(cp-*p++);
           n+=*(cp-*p++);
           n+=*(cp-*p);

           if (n<=max_no)
               r[i*Width+j] = max_no - n;
	   }
	int_to_uchar(r,InputImage,Imagesize);
	if(mid!=NULL) delete []mid;
	if(r!=NULL) delete []r;
}
/************************************************************************/
/*                   susan�˲���                                        */
/************************************************************************/
void Image_SUSAN_Filter(BYTE* InputImage , int ImgW , int ImgH, int bt, double dt) 
{
	int	 Imagesize = ImgW*ImgH;
	if (InputImage==NULL || Imagesize==0) return ;

    int      i,j,x,y,area,brightness,tmp,centre;
	int      mask_size = ((int)(1.5 * dt)) + 1;
	BYTE    *dpt,  *cp,  *bp,  *ip,  bbp[516];
    
	double   temp,total;

	BYTE *tmpImage = new BYTE[Imagesize];
    BYTE **RAtmpImage = new BYTE*[ImgH];
	BYTE **RAimg = new BYTE*[ImgH];
	BYTE *dp = new BYTE[(2*mask_size+1)*(2*mask_size+1)];
	//initialize 
	memcpy(tmpImage, InputImage, Imagesize);
	
	RAtmpImage[0] = InputImage;
	RAimg[0]      = tmpImage;
	for(i=1 ; i<ImgH ;i++)
	{
		RAtmpImage[i] = RAtmpImage[i-1] + ImgW;
		RAimg[i]      =      RAimg[i-1] + ImgW;
	}

	//�������ұ�
    bbp[0]=bbp[1]=bbp[515]=bbp[514]=0;
    bp=bbp+258;    
   
    for(int k=-256;k<257;k++)
	{
        temp=((float)k)/((float)bt);
        temp=temp*temp;
        temp=100.0*exp(-temp)+0.5;
        bbp[k+258]= (BYTE)(temp);
	}

	dpt = dp;
	temp   = -(dt*dt);
    for(i=-mask_size; i<=mask_size; i++)
	{
        for(j=-mask_size; j<=mask_size; j++)
		{
            x = (int) (100.0 * exp( ((float)((i*i)+(j*j))) / temp )+0.5);
            *dpt++ = (BYTE)x;
		}
	}  //�������ұ����
	//��ʼ�˲�
    for(j=mask_size;  j<ImgH-mask_size;  j++)
	{
		for(i=mask_size;  i<ImgW-mask_size;  i++) 
		{
			area   = 0;
			total  = 0;
			dpt    = dp;
            centre = RAimg[j][i];
		    cp     = bp + centre;
            for(y=-mask_size; y<=mask_size; y++)
			{
			    ip = RAimg[j+y] + i - mask_size;
                for(x=-mask_size; x<=mask_size; x++)
				{
				    brightness = *ip++;
                    tmp = (*dpt++) * (*(cp-brightness));
                    area += tmp;
                    total += tmp * brightness;
				}    			  
			}
            tmp = area-10000;
            if(tmp>30000)
			{
			   RAtmpImage[j][i] = (BYTE)( (total-(centre*10000))/tmp+0.5 );
			}
		    else
			{
			    int p[8],k,l,tmp;
                p[0]=RAimg[j-1][i-1];  p[1]=RAimg[j-1][ i ];  p[2]=RAimg[j-1][i+1];
				p[3]=RAimg[ j ][i-1];                         p[4]=RAimg[ j ][i+1]; 
				p[5]=RAimg[j+1][i-1];  p[6]=RAimg[j+1][i-1];  p[7]=RAimg[j+1][i+1];

                for(k=0; k<7; k++)
                    for(l=0; l<(7-k); l++)
						if (p[l]>p[l+1])
						{
                             tmp=p[l]; p[l]=p[l+1]; p[l+1]=tmp;
						}
                RAtmpImage[j][i] = (p[3]+p[4])/2;
			}
		}

	}
	//�˲�����
	if(dp!=NULL)  delete []dp;
	if(RAimg!=NULL)  delete []RAimg;
	if(tmpImage!=NULL)  delete []tmpImage;
	if(RAtmpImage!=NULL)  delete []RAtmpImage;
}

#define  M_E 2.7182818284590

// VGaussianKernel Generate a 1D Gaussian kernel.
BOOL VGaussianKernel (int ncoeffs, float *coeffs, float s)
{
    float sum, x, b = (float)( -1.0 / (2.0 * s * s) );
    float *p1, *p2, *pend = coeffs + ncoeffs;

    // Setup depends on whether the number of coefficients asked for is odd or even:
    if (ncoeffs & 1) // ncoeffs is a odd number
	{
		p1    = & coeffs[ncoeffs / 2];
		*p1-- = 1.0;
		p2    = p1 + 2;
		x     = 1.0;
		sum   = 0.5;
    } 
	else
	{
		p1  = & coeffs[ncoeffs / 2];
		p2  = p1 + 1;
		x   = 0.5;
		sum = 0.0;
    }

    // Fill the vector with coefficients, then make them sum to 1.0:
    while (p2 < pend) 
	{
		sum += *p1-- = *p2++ = (float)exp (x * x * b);
		x += 1.0;
    }
    sum += sum;
    for (p1 = coeffs; p1 < pend; p1++)
		*p1 /= sum;

    return TRUE;
}
// VGaussianD1Kernel Generate a kernel for the first derivative of a 1D Gaussian.
BOOL VGaussianD1Kernel (int ncoeffs, float *coeffs, float s)
{
    float x, a = (float)( 1.0 / (sqrt (2.0 * PI) * s * s * s) );
    float b = (float)( -1.0 / (2.0 * s * s) );
    float *p1, *p2, *pend = coeffs + ncoeffs;

    // Setup depends on whether the number of coefficients asked for is odd or even:
    if (ncoeffs & 1) 
	{
		p1 =  coeffs + ncoeffs/ 2;
		*p1-- = 0.0;
		p2 = p1 + 2;
		x = 1.0;
    }
	else 
	{
		p1 = & coeffs[ncoeffs / 2];
		p2 = p1 + 1;
		x = 0.5;
    }
    
    while (p2 < pend) // Fill the vector with coefficients:
	{
		/*float tx = x - 0.5f;  // this code is get derivation from integration of small area
		*p2 = 0;
		for(int i=0; i <11; i++,tx += 0.1f)
		{
			*p2 += a * tx * (float)exp (tx * tx * b);
		}
		*p2 /= 10; //*/
		*p2 = a * x * (float)exp (x * x * b);
		*p1-- = -*p2++;
		x += 1.0;
    }
    return TRUE;
}
//  VGaussianSplineKernel
//
//  Generate a 1D Gaussian kernel with its tails splined to zero.
//  The kernel is Gaussian over the domain [-2s,2s].
//  Over [-4s,-2s) and over (2s,4s], it splines to zero.
//  Beyond -4s and 4s it is zero.
//  Thus ncoeffs should be at least 8s + 1 to get a smoothly splined kernel.
// 
//  Based on Mark Nitzberg, David Mumford, Takahiro Shiota, `Filtering,
//  Segmentation and Depth', Lecture Notes in Computer Science 662,
//  Springer-Verlag, 1993.
BOOL VGaussianSplineKernel (int ncoeffs, float *coeffs, float s)
{
    float sum, x, t, b = (float)( -1.0 / (2.0 * s * s) );
    float k2 = (float)( 2.0 * s ) , k4 = (float) ( 4.0 * s );
    float *p1, *p2, *pend = coeffs + ncoeffs;

    // Setup depends on whether the number of coefficients asked for is odd or even:
    if (ncoeffs & 1) 
	{
		p1    = & coeffs[ncoeffs / 2];
		*p1-- = 1.0;
		p2    = p1 + 2;
		x     = 1.0;
		sum   = 0.5;
    }
	else
	{
		p1  = & coeffs[ncoeffs / 2];
		p2  = p1 + 1;
		x   = 0.5;
		sum = 0.0;
    }

    // Fill the vector with coefficients, then make them sum to 1.0:
    while (p2 < pend) 
	{
		if (x <= k2)
			t = (float)exp (x * x * b);
		else if (x <= k4) 
		{
			t  = (float)( 4.0 - x / s );
			t *= (float)( 1.0 / (16.0 * M_E * M_E) * t * t * t );
		} 
		else 
			t = 0.0;
		sum += *p1-- = *p2++ = t;
		x   += 1.0;
    }
    sum += sum;
    for (p1 = coeffs; p1 < pend; p1++)
		*p1 /= sum;

    return TRUE;
}

/************************************************************************/
/*                                                                      */
/************************************************************************/
FloatImage::FloatImage()
{
	lp_data   = NULL;
	lp_AddRow = NULL;
	Width     = Height = 0;
}
FloatImage::~FloatImage()
{
	DeleteData();
}
FloatImage::FloatImage(unsigned int w, unsigned int h)
{
	Width  = w;
	Height = h;
	ImageSize = w * h;
	lp_data   = NULL;
	lp_AddRow = NULL;
	Construct(w, h);
}
void FloatImage::Construct(unsigned int w, unsigned int h)
{
	DeleteData();
	if( ImageSize > 0 ) lp_data   = new  float[ImageSize];
	memset(lp_data, 0, ImageSize* sizeof(float) );
	if( Height    > 0 ) lp_AddRow = new float*[Height];
	lp_AddRow[0] = lp_data;
	for( UINT i= 1; i < Height; i++)
		lp_AddRow[i] = lp_AddRow[i- 1] + Width;
	
}
void FloatImage::DeleteData()
{
	if(lp_data   != NULL) delete []lp_data;
	if(lp_AddRow != NULL) delete []lp_AddRow;
}

/************************************************************************/
/*                                                                      */
/************************************************************************/
CImage::CImage()
{
	m_ImageSize=m_ImageHeight=m_ImageWidth=0;
	m_lpDibArray=NULL;
	RowAddress=NULL;
}
CImage::~CImage()
{
	DeleteData1();
}
CImage::CImage(unsigned int w, unsigned int h)
{
	m_ImageWidth = w;
	m_ImageHeight = h;
	m_ImageSize = w * h;
	m_lpDibArray   = NULL;
	RowAddress = NULL;
	Construct1(w, h);
}
void CImage::Construct1(unsigned int w, unsigned int h)
{
	DeleteData1();
	if( m_ImageSize > 0 ) m_lpDibArray   = new  unsigned char[m_ImageSize];
	memset(m_lpDibArray, 0, m_ImageSize* sizeof(unsigned char) );
	if( m_ImageHeight   > 0 ) RowAddress = new unsigned char*[m_ImageHeight];
	RowAddress[0] = m_lpDibArray;
	for( UINT i= 1; i < m_ImageHeight; i++)
		RowAddress[i] = RowAddress[i- 1] +m_ImageWidth;
}
void CImage::DeleteData1()
{
	if(m_lpDibArray!= NULL) delete []m_lpDibArray;
	if(RowAddress != NULL) delete []RowAddress;
}

//***********************************************************
//������� ͼ���Ե��ȡ
//�������ƣ�                       
//           Make_Rfilter_Mask
//������;��
//           ����R filter[ R(x) = exp(-|x|/��)/(2��) ] 
//			 ��ԭʼ���ݶ�ģ�� 
//˵    ��:  ���� Muhittin Gokmen ������" �˦�-Space Representation 
//			 of Images and Geeralized Edge Detector", PAMI Vol.19.NO.6 JUNE 1997
//ԭʼ���ߣ� ½  ��  ΰ
//ԭʼ���ڣ� 19/12/1999
//***********************************************************
void Make_Rfilter_Mask(int masksize, float lambda,  
					   float *lpDmask, float *lpmask, float maxresponse)
{
    int	  i, hM= masksize/ 2;    
    float gN0 = 0.5f* maxresponse/ lambda;			 // 
    float gD0 = 0.5f* maxresponse/(lambda * lambda); //

    for(lpDmask[hM] = 0.0f, i = -hM; i < 0; i++)
    {
		lpDmask[hM+ i] = -float( gD0* exp(i/ lambda) );
		lpDmask[hM- i] = -lpDmask[hM+ i];
    }
	if( lpmask != NULL )
	{		
		for(lpDmask[hM] = gN0, i = -hM; i < 0; i++)
			lpmask[hM+ i] = lpmask[hM- i] = float( gN0 * exp(i/ lambda) );
    }
}
//***********************************************************
//������� ͼ���Ե��ȡ
//�������ƣ�                       
//           Make_Second_Rfilter_Mask
//������;��
//           ����2 order R filter[ R(x) = exp(-|x|/sqrt(2)/pow(��,0.25))
//                    /(2*pow(��,0.25))*cos(-|x|/sqrt(2)/pow(��,0.25)- ��/4)]
//			 ��ԭʼ���ݶ�ģ�� 
//˵    ��:  ���� Muhittin Gokmen ������" �˦�-Space Representation 
//			 of Images and Geeralized Edge Detector", PAMI Vol.19.NO.6 JUNE 1997
//ԭʼ���ߣ� ½  ��  ΰ
//ԭʼ���ڣ� 19/12/1999
//***********************************************************
void Make_Second_Rfilter_Mask(int masksize, float lambda, float tau, 
							  float *lpDmask, float *lpmask, float maxresponse)
{
	int   i, hM= masksize/ 2;
    float Tg0, Tg1;
    
    float lambdaq    = float(sqrt(sqrt(lambda)));
    float lambdasqrt = float(sqrt(lambda));
    float cons2      = float(sqrt(2.0)* lambdaq);

    float gN0 = 0.5f* maxresponse/ lambda;
    float gN1 = 0.5f* maxresponse/ lambdaq;
    float gD0 = 0.5f* maxresponse/ (lambda * lambda);
    float gD1 = 0.5f* maxresponse/ lambdasqrt;

	for(lpDmask[hM] = 0.0f, i = -hM; i < 0; i++)
    {
		Tg0 = float( gD0 * exp(i/ lambda) );
		Tg1 = float( gD1 * exp(i/ cons2) * cos( i/ cons2) );
		lpDmask[hM- i] = (1.0f - tau) * Tg0 + tau * Tg1;
		lpDmask[hM+ i] = -lpDmask[hM- i];
    }
	if( lpmask != NULL )
	{		
		for(lpDmask[hM] = gN0, i = -hM; i < 0; i++)
		{
			Tg0 = float( gN0 * exp(i/ lambda) );
			Tg1 = float( gN1 * exp(i/ cons2) * cos(-i/ cons2 - PI/4.0) );
			lpmask[hM+ i] = lpmask[hM- i] = (1.0f - tau) * Tg0 + tau * Tg1;
		}
    }
}
//***********************************************************
//������� ͼ���Ե��ȡ
//�������ƣ�                       
//           Make_GED_filter_mask
//������;��
//           ����General Edge Detector filter ��ԭʼ���ݶ�ģ�� 
//˵    ��:  ���� Muhittin Gokmen ������" �˦�-Space Representation 
//			 of Images and Geeralized Edge Detector", PAMI Vol.19.NO.6 JUNE 1997
//ԭʼ���ߣ� ½  ��  ΰ
//ԭʼ���ڣ� 20/12/1999  ������Żع���!!!
//***********************************************************
void Make_GED_filter_mask(int masksize, float lambda, float tau, 
						  float *lpDmask, float *lpmask, float maxresponse)
{
    int   i, hM= masksize/ 2;
    float gN0, gD0, d1, d2, dSqrt, x, quarta;
    float a = lambda * tau;
    float b = lambda * ( 1.0f - tau);
    float discrim = b * b - 4.0f * a;
    float abs_discrim = (discrim > 0.0) ? discrim : -discrim;
    
    if (abs_discrim < 0.01) // case III
	{	
		d1 = (float)sqrt(b / (2.0 * a));
		d2 = 1.0f/ d1;

		gN0 = 0.5f* maxresponse/ b;
		gD0 = 0.5f* d1 * maxresponse / b;

		for(lpDmask[hM] = 0.0f, i = -hM; i < 0; i++)
		{
			lpDmask[hM+ i] = float( gD0 * i * exp(i* d1) );
			lpDmask[hM- i] = -lpDmask[hM+ i];
		}
		if( lpmask != NULL )
		{		
			for(lpDmask[hM] = gN0, i = -hM; i < 0; i++)
				lpmask[hM+ i] = lpmask[hM- i] = float( gN0 * exp(i* d1) * ( d2 - i) );
		}
	}
    else if (fabs(a) < 0.000001) // case V -  membrane 
	{
		d1  = (float)sqrt(b);
		gN0 = maxresponse / d1;
		gD0 = maxresponse / b;

		for(lpDmask[hM] = 0.0f, i = -hM; i < 0; i++)
		{
			lpDmask[hM+ i] = -float( gD0 * exp(i/ d1) );
			lpDmask[hM- i] = -lpDmask[hM+ i];
		}
		if( lpmask != NULL )
		{		
			for(lpDmask[hM] = gN0, i = -hM; i < 0; i++)
				lpmask[hM+ i] = lpmask[hM- i] = float( gN0 * exp(i/ d1) );
		}
	}
    else if (fabs(b) < 0.000001) // case IV -  Plate  
	{
		dSqrt = (float)sqrt(a);
		d2	  = (float)sqrt(dSqrt);
		d1    = 1.0f / (float)( sqrt(2.0) * d2);
		gN0   = maxresponse / (4.0f * dSqrt * d1);
		gD0   = maxresponse / (2.0f * dSqrt);

		for(lpDmask[hM] = 0.0f, x= -d1, i = 1; i <= hM; i++, x -= d1)
		{
			lpDmask[hM- i] = float( gD0 * exp(x)* sin(x) );
			lpDmask[hM+ i] = -lpDmask[hM- i];
		}
		if( lpmask != NULL )
		{		
			for(lpDmask[hM] = gN0, x= -d1, i = 1; i <= hM; i++, x -= d1)
				lpmask[hM+ i] = lpmask[hM- i] = float( gN0 * exp(x)* (cos(x) - sin(x)) );
		}
	} 
    else if (discrim > 0.0 ) // case I
	{
		d1 = (float)sqrt((b + sqrt(discrim))/(2.0 * a));
		d2 = (float)sqrt((b - sqrt(discrim))/(2.0 * a));

		gN0 = 1.0f * maxresponse/(2.0f* a* (d2* d2 - d1* d1));
		gD0 = gN0;
    
		for(lpDmask[hM] = 0.0f, i = -hM; i < 0; i++)
		{
			lpDmask[hM+ i] = -float( gD0 * (exp(i* d2) - exp(i* d1)) );
			lpDmask[hM- i] = -lpDmask[hM+ i];
		}
		if( lpmask != NULL )
		{		
			for(lpDmask[hM] = gN0, i = -hM; i < 0; i++)
				lpmask[hM+ i] = lpmask[hM- i] = float( gN0 * ( exp(i* d1)/ d1 - exp(i* d2)/ d2 ) );
		}
	} 
    else if (discrim < 0.0) // case II
	{
		dSqrt  = (float)sqrt(a);
        quarta = (float)sqrt(dSqrt);
		d1	   = (float)( cos(0.5* atan(sqrt(4.0*a/(b*b) - 1)))) / quarta ;
		d2     = (float)( sin(0.5* atan(sqrt(4.0*a/(b*b) - 1)))) / quarta ;
		gN0    = maxresponse/ (4.0f * dSqrt);
		gD0    = maxresponse/ (4.0f * a * d1 * d2);
    
		for(lpDmask[hM] = 0.0f, i = -hM; i < 0; i++)
		{
			lpDmask[hM+ i] = float( gD0 * exp(i* d2) * sin(i* d1) );
			lpDmask[hM- i] = -lpDmask[hM+ i];
		}
		if( lpmask != NULL )
		{		
			for(lpDmask[hM] = gN0, i = -hM; i < 0; i++)
				lpmask[hM+ i] = lpmask[hM- i] = float( gN0 * exp(i* d2)*( cos(i* d1)/ d2 + sin(i* d1)/ d1 ) );
		}
	} 
}

//***********************************************************
//������� ͼ���Ե��ȡ
//�������ƣ�                       
//           Image_Edge_Canny_EstimateNoise
//������;��
//           Canny edge detector using function
//˵    ��:
//           ���õ������������صĴ���, ��Ч�������.
//ԭʼ���ߣ� ����  ½��ΰ
//ԭʼ���ڣ� �� ֪ �� 
//�޸����ߣ� ½ �� ΰ
//�޸����ڣ� 19/12/1999 
//***********************************************************
float Image_Edge_Canny_EstimateNoise(FloatImage &mag, float amax, float q)
{
	const int Length = 1000;
    int t, histogram[Length];
	int histogram_size = Length * sizeof(int);
    int maxValue = Length - 1;
    if (q < 0.0 || q > 1.0) 
	{
		CString msg;
		msg.Format("VCanny: Noise estimation fraction (%f) is outside [0,1]", q);
		AfxMessageBox(msg);
	    q = 0.6f;
    }

    memset(histogram, 0, histogram_size);

    // Build histogram: 
    for (unsigned int i = 0; i < mag.ImageSize; i++ )
		histogram[(int) (mag.lp_data[i] / amax * maxValue )]++;

    // Find the index "i" such that percentile % of the pixels have indices 
    // (in the histogram) less than "i" (roughly speaking): 
	t = int( q * mag.ImageSize + 0.5);
    for (i = 0; i < Length; i++ ) 
	{
		t -= histogram[i++];
		if( t <= 0 ) break;
	}

    // Determine low threshold: 
    return float(i * amax / Length);
} 

//***********************************************************
//������� ͼ���Ե��ȡ
//�������ƣ�                       
//           Image_Edge_Canny_NonMaximalSuppress
//������;��
//           Canny edge detector using function
//˵    ��:
//           ���õ������������صĴ���, ��Ч�������.
//			 For each pixel, if the gradient magnitude is a local
//			 (3x3 neighbourhood) maxima in the gradient orientation, then
//			 output the gradient magnitude (appropriately scaled) and
//			 optionally the gradient orientation, otherwise output zero.
//ԭʼ���ߣ� ����  ½��ΰ
//ԭʼ���ڣ� �� ֪ �� 
//�޸����ߣ� ½ �� ΰ
//�޸����ڣ� 19/12/1999 
//***********************************************************
void Image_Edge_Canny_NonMaximalSuppress 
				  (CImage     &dest, 
				   FloatImage &xgrad, 
				   FloatImage &ygrad,
				   FloatImage &mag,
				   float	   amax)
{
    int x, y;
    float g;               // center pixel gradient magnitude value
    float g1, g2;          // interpolated gradient magnitude values
    float ga, gb, gc, gd;  // nearest pixels' gradient magnitude values 
    float ux, uy;          // components of edge normal vector 
    float dx, dy;
    float nf = 255 / amax; // output pixel normalizing factor     

	int x1 = 1,              y1 = 1;
	int x2 = mag.Width - 1,  y2 = mag.Height - 1;
    // For each pixel (except those on the image boundary), check if 
    // the gradient magnitude is a local (3x3 neighbourhood) maxima in 
    // the gradient orientation:
    for (y = y1; y < y2; y++) 
	{
		for (x = x1; x < x2; x++)
		{
			dx = xgrad.lp_AddRow[y][x];
			dy = ygrad.lp_AddRow[y][x];
			
			if (dx == 0.0f && dy == 0.0f) continue;

			if (fabs (dy) > fabs (dx)) 
			{
				ux = float(fabs (dx) / fabs (dy) );
				uy = 1.0f;
					
				gb = mag.lp_AddRow[y+1][x];
				gd = mag.lp_AddRow[y-1][x];
					
				if( dx*dy < 0 ) 
				{
					ga = mag.lp_AddRow[y+1][x-1];
					gc = mag.lp_AddRow[y-1][x+1];
				}
				else
				{
					ga = mag.lp_AddRow[y+1][x+1];
					gc = mag.lp_AddRow[y-1][x-1];
				}
			}
			else 
			{			
				ux = float( fabs (dy) / fabs (dx) );
				uy = 1.0f;

				gb = mag.lp_AddRow[y][x+1];
				gd = mag.lp_AddRow[y][x-1];

				if( dx*dy < 0 ) 
				{
					ga = mag.lp_AddRow[y-1][x+1];
					gc = mag.lp_AddRow[y+1][x-1];
				}
				else
				{
					ga = mag.lp_AddRow[y+1][x+1];
					gc = mag.lp_AddRow[y-1][x-1];
				}
			}
			g1 = (ux * ga) + (uy - ux) * gb;
			g2 = (ux * gc) + (uy - ux) * gd;
			g  = mag.lp_AddRow[y][x];
			
			if (g > g1 && g >= g2)// Store magnitude in output image while mapping [0,amax] to [1,255]:
			{				
				dest.RowAddress[y][x] = BYTE (mag.lp_AddRow[y][x] * nf + 0.5);
			}
		}
    }	
}
//***********************************************************
//������� ͼ���Ե��ȡ
//�������ƣ�                       
//           Image_Edge_Canny
//������;��
//           ��Ե��ȡ
//˵    ��:
//           ���õ������������صĴ���, ��Ч�������.
//ԭʼ���ߣ� ����  ½��ΰ
//ԭʼ���ڣ� �� ֪ �� 
//�޸����ߣ� ½ �� ΰ
//�޸����ڣ� 19/12/1999 
//***********************************************************
void    Image_Edge_Canny(CImage		&InImage,
						 float		sigma,
						 float		&noise,						 
						 float		lambda,
						 float		tau,
						 int		MaskType,
						 FloatImage *ori)
{
    
    //Image_Edge_Canny(TempImage, 2, noise, 4, 0, 3);
	int i, j, k;
	int w = InImage.m_ImageWidth;
	int h = InImage.m_ImageHeight;

    FloatImage xgrad(w, h), ygrad(w, h), mag(w, h);
    float amax = 0.0f;
    
    // Find the gradient (after Gaussian smoothing) in x and y directions:
    int hMask_Size = int( 3.0 * sigma + 0.5);
	int Mask_Size  = 2 * hMask_Size + 1;
	int x1 = hMask_Size,         y1 = hMask_Size;
	int x2 = w - hMask_Size- 1,  y2 = h - hMask_Size- 1;


    //1���˲�
	float *Mask_x = new float[Mask_Size], *lp_mx;
	float *Mask_y = new float[Mask_Size], *lp_my;
	switch(MaskType)
	{
	case 0:
		VGaussianD1Kernel( Mask_Size, Mask_x, sigma);
		break;
	case 1:
		Make_Rfilter_Mask(Mask_Size, lambda, Mask_x);
		break;
	case 2:
		Make_Second_Rfilter_Mask(Mask_Size, lambda, tau, Mask_x);
		break;
	case 3:
		Make_GED_filter_mask(Mask_Size, lambda, tau, Mask_x);
		break;
	}	


    //2��
	for(i = 0; i< Mask_Size; i++) 
        Mask_y[i] = Mask_x[i];
 	//for(i = 0; i< Mask_Size; i++) TRACE("Mask_x[%2d] = %15.10f  Mask_y[%2d] = %15.10f\n", i, Mask_x[i], i, Mask_y[i]);
   
	register float dx, dy;    //Ϊ�˼ӿ����㣿���õ��˼Ĵ�����


//     FILE *fp1,*fp2,*fp3;
//     fp1 = fopen("GED_dx.txt","wt");
//     fp2 = fopen("GED_dy.txt","wt");
//     fp3 = fopen("GED_dx_dy.txt","wt");
// 
    //˳���dx��dy����ֵ���

	for (j = y1; j < y2; j++)
	{
		for (i = x1; i < x2; i++)
		{
			lp_mx = Mask_x;  lp_my = Mask_y; 
			dx = dy = 0.0f;
			for(k = -hMask_Size; k<= hMask_Size; k++)
			{
				dx += InImage.RowAddress[j][i+ k] * *lp_mx++;
				dy += InImage.RowAddress[j+ k][i] * *lp_my++;
			}
			xgrad.lp_AddRow[j][i] = dx;
			ygrad.lp_AddRow[j][i] = dy;
			mag.lp_AddRow[j][i] = (float)sqrt(dx* dx + dy* dy);      //�ݶȷ�ֵ��������

//             fprintf(fp1,"%f  ",dx);
//             fprintf(fp2,"%f  ",dy);
//             fprintf(fp3,"%f  ",mag.lp_AddRow[j][i]);

			if( mag.lp_AddRow[j][i] > amax)
				amax = mag.lp_AddRow[j][i];
		}
//         fprintf(fp1,"\n");
//         fprintf(fp2,"\n");
//         fprintf(fp3,"\n");
	}	

//     fclose(fp1);
//     fclose(fp2);
//     fclose(fp3);


    //������������ʲô�õģ�

	// Estimate noise 
	noise = Image_Edge_Canny_EstimateNoise(mag, amax, noise) *  255 / amax;
    // Non-maximal suppression:
	memset(InImage.m_lpDibArray, 0, InImage.m_ImageSize);
  
    //�Ǽ���ֵ����
    Image_Edge_Canny_NonMaximalSuppress(InImage, xgrad, ygrad, mag, amax);
  
/*	non_max_supp(magnitude, delta_x, delta_y, rows, cols, nms);*/
	
    

//     unsigned char *edge;
// 
//     short int *magnitude;
// 
//     unsigned char *nms;
// 
// 
//     if((nms = (unsigned char *) calloc(h*w,sizeof(unsigned char)))==NULL)
//     {
//         fprintf(stderr, "Error allocating the nms image.\n");
//         exit(1);
// 	}
// 
// 
//     if((edge = (unsigned char *) calloc(h*w,sizeof(unsigned char)))==NULL)
//     {
//         fprintf(stderr, "Error allocating the nms image.\n");
//         exit(1);
//     }
// 
//     if((magnitude = (short int *) calloc(h*w,sizeof(short int)))==NULL)
//     {
//         fprintf(stderr, "Error allocating the nms image.\n");
//         exit(1);
//     }


//     FILE *fp4;
//     fp4 = fopen("GED_non_max.txt","wt");



// 
//     //�ѽ������д��ȥ
//     for (i=0;i<h;i++)
//     {
//         for (j=0;j<w;j++)
//         {
//   //          fprintf(fp4,"%d  ",InImage.m_lpDibArray[i*w+j]);
//             magnitude[i*w+j] = (short int)InImage.m_lpDibArray[i*w+j];
//             if (magnitude>0)
//             {
//                 nms[i*w+j] = 128;
//             }
//             else
//             {
//                 nms[i*w+j] = 0;
//             }
//         }
//    //     fprintf(fp4,"\n");
//     }
  //  fclose(fp4);
// 
// 	/****************************************************************************
// 	* Use hysteresis to mark the edge pixels.
// 	****************************************************************************/
    // 
//     apply_hysteresis(magnitude, nms, h, w, 0.4, 0.95, edge);
// 
// 
//     for (i=0;i<h;i++)
//     {
//         for (j=0;j<w;j++)
//         {
//              InImage.m_lpDibArray[i*w+j] = edge[i*w+j];
//         }
//     }






    if( ori != NULL ) 
	{
		for (j = y1; j < y2; j++) 
		{
			for (i = x1; i < x2; i++)
			{
				if( InImage.RowAddress[j][i] == 0 ) continue;

				dx = xgrad.lp_AddRow[j][i];
				dy = ygrad.lp_AddRow[j][i];
				ori->lp_AddRow[j][i] = float( atan2 (dy, dx) );
			}
		}
	}
    if( Mask_x != NULL ) delete []Mask_x;
	if( Mask_y != NULL ) delete []Mask_y;


//     delete[] nms;
//     delete[] magnitude;
//     delete[] edge;
//    free(nms);
//    free(magnitude);
//    free(edge);

}


/////////////////////////////////////////////////////////////////////////      
//������ip,jp������ͼ���һά����    
//      lx��ͼ����    
//      ly��ͼ��߶�    
//      ����Ҫ��ǰ���ǰ׵㻹�Ǻڵ���������ԸĶ�    
//      ���ǰ���ǰ׵㣬������������֮������
//���ܣ�  ͼ��ɫ������ǰ�ɫ�ĵ�Ͳ���˲���     
//////////////////////////////////////////////////////////////////////////
void beforethin(unsigned char *ip, unsigned char *jp,    
                unsigned long lx, unsigned long ly)   
{   
    unsigned long i,j;   
    for(i=0; i<ly; i++)   
    {   
        for(j=0; j<lx; j++)   
        {   
            //����Ҫ��ǰ���ǰ׵㻹�Ǻڵ���������ԸĶ�    
            //���ǰ���ǰ׵㣬������������֮������    
            if(ip[i*lx+j]>0)   
                jp[i*lx+j]=1;   
            else   
                jp[i*lx+j]=0;   
        }   
    }   
}   
/////////////////////////////////////////////////////////////////////////    
//  ϸ���㷨1
//Hilditchϸ���㷨    
//���ܣ���ͼ�����ϸ��    
//������image������ͼ���һά����    
//      lx��ͼ����    
//      ly��ͼ��߶�    
//      �޷���ֵ   
//////////////////////////////////////////////////////////////////////////
 
void ThinnerHilditch(void *image, unsigned long lx, unsigned long ly)   
{   
    char *f, *g;   
    char n[10];   
    unsigned int counter;   
    short k, shori, xx, nrn;   
    unsigned long i, j;   
    long kk, kk11, kk12, kk13, kk21, kk22, kk23, kk31, kk32, kk33, size;   
    size = (long)lx * (long)ly;   
    g = (char *)malloc(size);   
   
    if(g == NULL)   
    {   
        printf("error in allocating memory!\n");   
        return;   
    }   
   
    f = (char *)image;   
    for(i=0; i<lx; i++)   
    {   
        for(j=0; j<ly; j++)   
        {   
            kk=i*ly+j;   
            if(f[kk]!=0)   
            {   
                f[kk]=1;   
                g[kk]=f[kk];   
            }   
        }   
    }   
   
    counter = 1;   
   
    do   
    {   
        printf("%4d*",counter);   
        counter++;   
        shori = 0;   
   
        for(i=0; i<lx; i++)   
        {   
            for(j=0; j<ly; j++)   
            {   
                kk = i*ly+j;   
                if(f[kk]<0)   
                    f[kk] = 0;   
                g[kk]= f[kk];   
            }   
        }   
   
        for(i=1; i<lx-1; i++)   
        {   
            for(j=1; j<ly-1; j++)   
            {   
                kk=i*ly+j;   
   
                if(f[kk]!=1)   
                    continue;   
   
                kk11 = (i-1)*ly+j-1;   
                kk12 = kk11 + 1;   
                kk13 = kk12 + 1;   
                kk21 = i*ly+j-1;   
                kk22 = kk21 + 1;   
                kk23 = kk22 + 1;   
                kk31 = (i+1)*ly+j-1;   
                kk32 = kk31 + 1;   
                kk33 = kk32 + 1;   
   
                if((g[kk12]&&g[kk21]&&g[kk23]&&g[kk32])!=0)   
                    continue;   
   
                nrn = g[kk11] + g[kk12] + g[kk13] + g[kk21] + g[kk23] +    
                    g[kk31] + g[kk32] + g[kk33];   
   
                if(nrn = 1)   
                {   
                    f[kk22] = 2;   
                    continue;   
                }   
   
                n[4] = f[kk11];   
                n[3] = f[kk12];   
                n[2] = f[kk13];   
                n[5] = f[kk21];   
                n[1] = f[kk23];   
                n[6] = f[kk31];   
                n[7] = f[kk32];   
                n[8] = f[kk33];   
                n[9] = n[1];   
                xx = 0;   
   
                for(k=1; k<8; k=k+2)   
                {   
                    if((!n[k])&&(n[k+1]||n[k+2]))   
                        xx++;   
                }   
   
                if(xx!=1)   
                {   
                    f[kk22] = 2;   
                    continue;   
                }   
   
                if(f[kk12] == -1)   
                {   
                    f[kk12] = 0;   
                    n[3] = 0;   
                    xx = 0;   
   
                    for(k=1; k<8; k=k+2)   
                    {   
                        if((!n[k])&&(n[k+1]||n[k+2]))   
                            xx++;   
                    }   
   
                    if(xx != 1)   
                    {   
                        f[kk12] = -1;   
                        continue;   
                    }   
   
                    f[kk12] = -1;   
                    n[3] = -1;   
                }   
   
                if(f[kk21]!=-1)   
                {   
                    f[kk22] = -1;   
                    shori = 1;   
                    continue;   
                }   
   
                f[kk21] = 0;   
                n[5] = 0;   
                xx = 0;   
   
                for(k=1; k<8; k=k+2)   
                {   
                    if((!n[k])&&(n[k+1]||n[k+2]))   
                    {   
                        xx++;   
                    }   
                }   
   
                if(xx == 1)   
                {   
                    f[kk21] = -1;   
                    f[kk22] = -1;   
                    shori =1;   
                }   
                else   
                    f[kk21] = -1;   
            }   
        }   
    }while(shori);   
   
    free(g);   
}   
   
/////////////////////////////////////////////////////////////////////////   
//  ϸ���㷨2 
//Pavlidisϸ���㷨    
//���ܣ���ͼ�����ϸ��    
//������image������ͼ���һά����    
//      lx��ͼ����    
//      ly��ͼ��߶�    
//      �޷���ֵ    
void ThinnerPavlidis(void *image, unsigned long lx, unsigned long ly)   
{   
    char erase, n[8];   
    char *f;   
    unsigned char bdr1,bdr2,bdr4,bdr5;   
    short c,k,b;   
    unsigned long i,j;   
    long kk,kk1,kk2,kk3;   
    f = (char*)image;   
   
    for(i=1; i<lx-1; i++)   
    {   
        for(j=1; j<ly-1; j++)   
        {   
            kk = i*ly + j;   
            if(f[kk])   
                f[kk] = 1;   
        }   
    }   
   
    for(i=0, kk1=0, kk2=ly-1; i<lx; i++, kk1+=ly, kk2+=ly)   
    {   
        f[kk1]=0;   
        f[kk2]=0;   
    }   
   
    for(j=0, kk=(lx-1)*ly; j<ly; j++,kk++)   
    {   
        f[j]=0;   
        f[kk]=0;   
    }   
   
    c=5;   
    erase =1;   
    while(erase)   
    {   
        c++;   
        for(i=1; i<lx-1; i++)   
        {   
            for(j=1; j<ly-1; j++)   
            {   
                kk=i*ly+j;   
                if(f[kk]!=1)   
                    continue;   
   
                kk1 = kk-ly -1;   
                kk2 = kk1 + 1;   
                kk3 = kk2 + 1;   
                n[3] = f[kk1];   
                n[2] = f[kk2];   
                n[1] = f[kk3];   
                kk1 = kk - 1;   
                kk3 = kk + 1;   
                n[4] = f[kk1];   
                n[0] = f[kk3];   
                kk1 = kk + ly -1;   
                kk2 = kk1 + 1;   
                kk3 = kk2 + 1;   
                n[5] = f[kk1];   
                n[6] = f[kk2];   
                n[7] = f[kk3];   
   
                bdr1 =0;   
                for(k=0; k<8; k++)   
                {   
                    if(n[k]>=1)   
                        bdr1|=0x80>>k;   
                }   
   
                if((bdr1&0252)== 0252)   
                    continue;   
                f[kk] = 2;   
                b=0;   
   
                for(k=0; k=7; k++)   
                {   
                    b+=bdr1&(0x80>>k);   
                }   
   
                if(b=1)   
                    f[kk]=3;   
   
                if((bdr1&0160)!=0&&(bdr1&07)!=0&&(bdr1&0210)==0)   
                    f[kk]=3;   
                else if((bdr1&&0301)!=0&&(bdr1&034)!=0&&(bdr1&042)==0)   
                    f[kk]=3;   
                else if((bdr1&0202)==0 && (bdr1&01)!=0)   
                    f[kk]=3;   
                else if((bdr1&0240)==0 && (bdr1&0100)!=0)   
                    f[kk]=3;   
                else if((bdr1&050)==0 && (bdr1&020)!=0)   
                    f[kk]=3;   
                else if((bdr1&012)==0 && (bdr1&04)!=0)   
                    f[kk]=3;   
            }   
        }   
   
        for(i=1; i<lx-1; i++)   
        {   
            for(j=1; j<ly-1; j++)   
            {   
                kk = i*ly + j;   
                if(!f[kk])   
                    continue;   
   
                kk1 = kk - ly -1;   
                kk2 = kk1 + 1;   
                kk3 = kk2 + 1;   
                n[3] = f[kk1];   
                n[2] = f[kk2];   
                n[1] = f[kk3];   
                kk1 = kk - 1;   
                kk2 = kk + 1;   
                n[4] = f[kk1];   
                n[0] = f[kk3];   
                kk1 = kk + ly -1;   
                kk2 = kk1 + 1;   
                kk3 = kk2 + 1;   
                n[5] = f[kk1];   
                n[6] = f[kk2];   
                n[7] = f[kk3];   
                bdr1 = bdr2 =0;   
   
                for(k=0; k=7; k++)   
                {   
                    if(n[k]>=1)   
                        bdr1|=0x80>>k;   
                    if(n[k]>=2)   
                        bdr2|=0x80>>k;   
                }   
   
                if(bdr1==bdr2)   
                {   
                    f[kk] = 4;   
                    continue;   
                }   
   
                if(f[kk]!=2)   
                    continue;   
   
                if((bdr2&0200)!=0 && (bdr1&010)==0 &&   
                    ((bdr1&0100)!=0 &&(bdr1&001)!=0 ||   
                    ((bdr1&0100)!=0 ||(bdr1 & 001)!=0) &&   
                    (bdr1&060)!=0 &&(bdr1&06)!=0))   
                {   
                    f[kk] = 4;   
                }   
   
                else if((bdr2&040)!=0 && (bdr1&02)==0 &&   
                    ((bdr1&020)!=0 && (bdr1&0100)!=0 ||   
                    ((bdr1&020)!=0 || (bdr1&0100)!=0) &&   
                    (bdr1&014)!=0 && (bdr1&0201)!=0))   
                {   
                    f[kk] = 4;   
                }   
   
                else if((bdr2&010)!=0 && (bdr1&0200)==0 &&   
                    ((bdr1&04)!=0 && (bdr1&020)!=0 ||   
                    ((bdr1&04)!=0 || (bdr1&020)!=0) &&   
                    (bdr1&03)!=0 && (bdr1&0140)!=0))   
                {   
                    f[kk] = 4;   
                }   
   
                else if((bdr2&02)!=0 && (bdr1&040)==0 &&   
                    ((bdr1&01)!=0 && (bdr1&04)!=0 ||   
                    ((bdr1&01)!=0 || (bdr1&04)!=0) &&   
                    (bdr1&0300)!=0 && (bdr1&030)!=0))   
                {   
                    f[kk] = 4;   
                }   
            }   
        }   
   
        for(i=1; i<lx-1; i++)   
        {   
            for(j=1; j<ly-1; j++)   
            {   
                kk = i*ly + j;   
                if(f[kk]!=2)   
                    continue;   
                kk1 = kk - ly -1;   
                kk2 = kk1 + 1;   
                kk3 = kk2 + 1;   
                n[3] = f[kk1];   
                n[2] = f[kk2];   
                n[1] = f[kk3];   
                kk1 = kk - 1;   
                kk2 = kk + 1;   
                n[4] = f[kk1];   
                n[0] = f[kk3];   
                kk1 = kk + ly -1;   
                kk2 = kk1 + 1;   
                kk3 = kk2 + 1;   
                n[5] = f[kk1];   
                n[6] = f[kk2];   
                n[7] = f[kk3];   
                bdr4 = bdr5 =0;   
                for(k=0; k=7; k++)   
                {   
                    if(n[k]>=4)   
                        bdr4|=0x80>>k;   
                    if(n[k]>=5)   
                        bdr5|=0x80>>k;   
                }   
                if((bdr4&010) == 0)   
                {   
                    f[kk] = 5;   
                    continue;   
                }   
                if((bdr4&040) == 0 && bdr5 ==0)   
                {   
                    f[kk] = 5;   
                    continue;   
                }   
                if(f[kk]==3||f[kk]==4)   
                    f[kk] = c;   
            }   
        }   
   
        erase = 0;   
        for(i=1; i<lx-1; i++)   
        {   
            for(j=1; j<ly-1; j++)   
            {   
                kk = i*ly +j;   
                if(f[kk]==2||f[kk] == 5)   
                {   
                    erase = 1;   
                    f[kk] = 0;   
                }   
            }   
        }   
    }   
}   
   
/////////////////////////////////////////////////////////////////////////   
/*  ϸ���㷨3 
* �������ƣ�
*     ThinnerRosenfeld
*
* ������
*   void*     image             ����ֵ��ͼ�����ǰ��ɫΪ1����ɫΪ0
*   unsigned  longlx             ��ͼ��Ŀ��
*   unsigned  longly             ��ͼ��ĸ߶�
*
* ����ֵ
*       ��
*
*�������ܣ�
*       �������ͼ�����ϸ�������ϸ�����ͼ��
***********************************************************/
void ThinnerRosenfeld(void *image, unsigned long lx, unsigned long ly)
{
    char *f, *g;
    char n[10];
    char a[5] = {0, -1, 1, 0, 0};
    char b[5] = {0, 0, 0, 1, -1};
    char nrnd, cond, n48, n26, n24, n46, n68, n82, n123, n345, n567, n781;
    short k, shori;
    unsigned long i, j;
    long ii, jj, kk, kk1, kk2, kk3, size;
    size = (long)lx * (long)ly;

    g = (char *)malloc(size);
    if(g==NULL)
    {
        printf("error in alocating mmeory!\n");
        return;
    }

    f = (char *)image;
    for(kk=0l; kk<size; kk++)
    {
        g[kk] = f[kk];
    }

    do
    {
        shori = 0;
        for(k=1; k<=4; k++)
        {
            for(i=1; i<lx-1; i++)
            {
                ii = i + a[k];

                for(j=1; j<ly-1; j++)
                {
                    kk = i*ly + j;

                    if(!f[kk])
                        continue;

                    jj = j + b[k];
                    kk1 = ii*ly + jj;

                    if(f[kk1])
                        continue;

                    kk1 = kk - ly -1;
                    kk2 = kk1 + 1;
                    kk3 = kk2 + 1;
                    n[3] = f[kk1];
                    n[2] = f[kk2];
                    n[1] = f[kk3];
                    kk1 = kk - 1;
                    kk3 = kk + 1;
                    n[4] = f[kk1];
                    n[8] = f[kk3];
                    kk1 = kk + ly - 1;
                    kk2 = kk1 + 1;
                    kk3 = kk2 + 1;
                    n[5] = f[kk1];
                    n[6] = f[kk2];
                    n[7] = f[kk3];

                    nrnd = n[1] + n[2] + n[3] + n[4]
                        +n[5] + n[6] + n[7] + n[8];
                    if(nrnd<=1)
                        continue;

                    cond = 0;
                    n48 = n[4] + n[8];
                    n26 = n[2] + n[6];
                    n24 = n[2] + n[4];
                    n46 = n[4] + n[6];
                    n68 = n[6] + n[8];
                    n82 = n[8] + n[2];
                    n123 = n[1] + n[2] + n[3];
                    n345 = n[3] + n[4] + n[5];
                    n567 = n[5] + n[6] + n[7];
                    n781 = n[7] + n[8] + n[1];

                    if(n[2]==1 && n48==0 && n567>0)
                    {
                        if(!cond)
                            continue;
                        g[kk] = 0;
                        shori = 1;
                        continue;
                    }

                    if(n[6]==1 && n48==0 && n123>0)
                    {
                        if(!cond)
                            continue;
                        g[kk] = 0;
                        shori = 1;
                        continue;
                    }

                    if(n[8]==1 && n26==0 && n345>0)
                    {
                        if(!cond)
                            continue;
                        g[kk] = 0;
                        shori = 1;
                        continue;
                    }

                    if(n[4]==1 && n26==0 && n781>0)
                    {
                        if(!cond)
                            continue;
                        g[kk] = 0;
                        shori = 1;
                        continue;
                    }

                    if(n[5]==1 && n46==0)
                    {
                        if(!cond)
                            continue;
                        g[kk] = 0;
                        shori = 1;
                        continue;
                    }

                    if(n[7]==1 && n68==0)
                    {
                        if(!cond)
                            continue;
                        g[kk] = 0;
                        shori = 1;
                        continue;
                    }

                    if(n[1]==1 && n82==0)
                    {
                        if(!cond)
                            continue;
                        g[kk] = 0;
                        shori = 1;
                        continue;
                    }

                    if(n[3]==1 && n24==0)
                    {
                        if(!cond)
                            continue;
                        g[kk] = 0;
                        shori = 1;
                        continue;
                    }

                    cond = 1;
                    if(!cond)
                        continue;
                    g[kk] = 0;
                    shori = 1;
                }
            }

            for(i=0; i<lx; i++)
            {
                for(j=0; j<ly; j++)
                {
                    kk = i*ly + j;
                    f[kk] = g[kk];
                }
            }
        }
    }while(shori);

    free(g);
}              
  
   
/////////////////////////////////////////////////////////////////////////    
//  ϸ���㷨4
//�����������ϸ��ϸ���㷨    
//���ܣ���ͼ�����ϸ��    
//������lpDIBBits������ͼ���һά����    
//      lWidth��ͼ��߶�    
//      lHeight��ͼ����    
//      �޷���ֵ    
BOOL WINAPI ThiningDIBSkeleton (LPSTR lpDIBBits, LONG lWidth, LONG lHeight)   
{      
    //ѭ������    
    long i;   
    long j;   
    long lLength;   
       
    unsigned char deletemark[256] = {   
        0,0,0,0,0,0,0,1,    0,0,1,1,0,0,1,1,   
        0,0,0,0,0,0,0,0,    0,0,1,1,1,0,1,1,   
        0,0,0,0,0,0,0,0,    1,0,0,0,1,0,1,1,   
        0,0,0,0,0,0,0,0,    1,0,1,1,1,0,1,1,   
        0,0,0,0,0,0,0,0,    0,0,0,0,0,0,0,0,   
        0,0,0,0,0,0,0,0,    0,0,0,0,0,0,0,0,   
        0,0,0,0,0,0,0,0,    1,0,0,0,1,0,1,1,   
        1,0,0,0,0,0,0,0,    1,0,1,1,1,0,1,1,   
        0,0,1,1,0,0,1,1,    0,0,0,1,0,0,1,1,   
        0,0,0,0,0,0,0,0,    0,0,0,1,0,0,1,1,   
        1,1,0,1,0,0,0,1,    0,0,0,0,0,0,0,0,   
        1,1,0,1,0,0,0,1,    1,1,0,0,1,0,0,0,   
        0,1,1,1,0,0,1,1,    0,0,0,1,0,0,1,1,   
        0,0,0,0,0,0,0,0,    0,0,0,0,0,1,1,1,   
        1,1,1,1,0,0,1,1,    1,1,0,0,1,1,0,0,   
        1,1,1,1,0,0,1,1,    1,1,0,0,1,1,0,0   
    };//������    
   
    unsigned char p0, p1, p2, p3, p4, p5, p6, p7;   
    unsigned char *pmid, *pmidtemp;   
    unsigned char sum;   
    int changed;   
    bool bStart = true;   
    lLength = lWidth * lHeight;   
    unsigned char *pTemp = (unsigned char *)malloc(sizeof(unsigned char) * lWidth * lHeight);   
       
    //    P0 P1 P2    
    //    P7    P3    
    //    P6 P5 P4    
   
    while(bStart)   
    {   
        bStart = false;   
        changed = 0;   
   
        //�������Ե��(����)    
        pmid = (unsigned char *)lpDIBBits + lWidth + 1;   
        memset(pTemp, (BYTE) 0, lLength);   
        pmidtemp = (unsigned char *)pTemp + lWidth + 1;   
        for(i = 1; i =lHeight -1; i++)   
        {   
            for(j = 1; j=  lWidth - 1; j++)   
            {   
                if( *pmid == 0)   
                {   
                    pmid++;   
                    pmidtemp++;   
                    continue;   
                }   
   
                p3 = *(pmid + 1);   
                p2 = *(pmid + 1 - lWidth);   
                p1 = *(pmid - lWidth);   
                p0 = *(pmid - lWidth -1);   
                p7 = *(pmid - 1);   
                p6 = *(pmid + lWidth - 1);   
                p5 = *(pmid + lWidth);   
                p4 = *(pmid + lWidth + 1);   
                   
                sum = p0 & p1 & p2 & p3 & p4 & p5 & p6 & p7;   
                if(sum == 0)   
                {   
                    *pmidtemp = 1;   
                }   
   
                pmid++;   
                pmidtemp++;   
            }   
            pmid++;   
            pmid++;   
            pmidtemp++;   
            pmidtemp++;   
        }   
           
        //���ڿ�ʼ����ɾ��    
        pmid = (unsigned char *)lpDIBBits + lWidth + 1;   
        pmidtemp = (unsigned char *)pTemp + lWidth + 1;   
   
        for(i = 1; i = lHeight -1; i++)   
        {   
            for(j = 1; j = lWidth - 1; j++)   
            {   
                if( *pmidtemp == 0)   
                {   
                    pmid++;   
                    pmidtemp++;   
                    continue;   
                }   
   
                p3 = *(pmid + 1);   
                p2 = *(pmid + 1 - lWidth);   
                p1 = *(pmid - lWidth);   
                p0 = *(pmid - lWidth -1);   
                p7 = *(pmid - 1);   
                p6 = *(pmid + lWidth - 1);   
                p5 = *(pmid + lWidth);   
                p4 = *(pmid + lWidth + 1);   
                   
                p1 *= 2;   
                p2 *= 4;   
                p3 *= 8;   
                p4 *= 16;   
                p5 *= 32;   
                p6 *= 64;   
                p7 *= 128;   
   
                sum = p0 | p1 | p2 | p3 | p4 | p5 | p6 | p7;   
                if(deletemark[sum] == 1)   
                {   
                    *pmid = 0;   
                    bStart = true;   
                }   
   
                pmid++;   
                pmidtemp++;   
            }   
   
            pmid++;   
            pmid++;   
            pmidtemp++;   
            pmidtemp++;   
        }   
    }   
   
    return true;   
}  


void Image_binary(unsigned char*image,int width,int height,int thresh)
{
	int linebyte=(width+3)/4*4;
	for(int y=0;y<height;y++)
		for (int x=0;x<width;x++)
		{
			int k=linebyte*y+x;
			int temp = image[k]>thresh ? 255:0;
			image[k]=temp;
		}
}

void Edge_max_block(unsigned char*image,int width,int height,int thresh)
{
	int linebyte=(width+3)/4*4;
	int skip=3;
	unsigned char*temp=new unsigned char[linebyte*height];
	memset(temp,0,linebyte*height*sizeof(unsigned char));
	
	for(int y=0;y<height;y=y+skip)
	{
		for (int x=0;x<width;x=x+skip)
		{
			int k=linebyte*y+x;

			int temp = image[k]>thresh ? 255:0;
		
			image[k]=temp;
		}
	}
}