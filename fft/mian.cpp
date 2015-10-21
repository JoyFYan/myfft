#include <stdio.h>
#include "math.h"
#include <stdlib.h>
#include<iostream> 
#include "opencv2/core/core.hpp" 
#include "opencv2/highgui/highgui.hpp" 
#include "FT_PubFunction.h"
using namespace cv;
#define W 256
#define H 256
void main()
{
	int i, j;
	//Mat img = imread("C:\\cameraman.tif");// 创建一个名为 "游戏原画"窗口
	IplImage* img = cvLoadImage("C:\\cameraman.tif", -1);
	int m = img->height;
	int n = img->width;
	namedWindow("1");// 在窗口中显示游戏原画 
	cvShowImage("1", img);
	waitKey(60);
	//uchar *p;
	
	/*p = new uchar *[m];
	t = new double *[m];
	for (i = 0; i<m; i++)
	{
		p[i] = new uchar[n];
		t[i] = new double[n];
	}*/
	uchar *p = new uchar[m*n];
	//double *t = new double[m*n];
	uchar *t = new uchar[m*n];
	uchar *t0 = new uchar[m*n];
	uchar *t1 = new uchar[m*n];
	uchar *modle = new uchar[m*n];
	uchar *ptr;

	/*for (i = 0; i<m; i++)
	{
		ptr = (uchar*)img->imageData + i*img->widthStep;
		for (j = 0; j<n; j++)
		{
			p[i][j] = (uchar)*(ptr + j);
		}
	}*/

	int k = 0;
	for (i = 0;i<m;i++)
	{
		ptr = (uchar*)img->imageData + i*img->widthStep;
		for (j = 0;j<n;j++)
		{
			unsigned int k1 = m*i + j;
			p[k1]=(uchar)*(ptr + j);    //实部    源图像可能不符合宽度为4
			//t[k1] = 0;                 //虚部
		}
	}

	FFT_all(p, t, t0, n,m);
	//for (i = 0;i < m;i++)//构建滤波模板，点乘
	//{
	//	for (j = 0;j < n;j++)
	//	{
	//		if ((sqrt((i + 1 - m / 2) * (1 + i - m / 2) + (j + 1 - n / 2) * (j + 1 - n / 2)))>50)modle[i * n + j] = 1;
	//		//D = sqrt((m - M) ^ 2 + (n - N) ^ 2);
	//		//H(m, n) = exp((-D ^ 2) / (2 * (50) ^ 2));
	//		else modle[i * n + j] = 0;
	//		//printf("%f,  ", Modle[i * n + j]);

	//	}
	//	printf("\n");
	//}
	//FFT(p, t, n - 1, m - 1);
	//fftshift_2uc(t, t1,n,m);
	
	IplImage *copy = cvCreateImage(cvGetSize(img), 8, 1);



	for (i = 0; i<m; i++)
	{
		for (j = 0; j<n; j++)
		{
			unsigned int k1 = m*i + j;
			
			cvSetReal2D(copy, i, j, t0[k1]);
		}
	}





	// 显示图像副本
	cvNamedWindow("copy", 1);
	cvShowImage("copy", copy);
	cvWaitKey();


	//
}