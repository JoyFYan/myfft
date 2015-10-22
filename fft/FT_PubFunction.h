#ifndef C_FT_PubFunction_H
#define C_FT_PubFunction_H
#include "stdafx.h"
#include "fftw3.h"   //添加fftw头文件
#include<string.h>
#include "math.h"
#include "stdio.h"
#include "stdlib.h"

void factorize(int n, int *nFact, int fact[]);
void transTableSetup(int sofar[], int actual[], int remain[],int *nFact,int *nPoints);
void permute(int nPoint, int nFact,int fact[], int remain[],double xRe[], double xIm[],double yRe[], double yIm[]);
void initTrig(int radix);
void fft_4(double aRe[], double aIm[]);
void fft_5(double aRe[], double aIm[]);
void fft_8();
void fft_10();
void fft_odd(int radix);
void twiddleTransf(int sofarRadix, int radix, int remainRadix,double yRe[], double yIm[]);
void fft(int n, double xRe[], double xIm[],double yRe[], double yIm[]);

void ifft(int n,double yRe[],double yIm[],double xRe[],double xIm[]);

void Fourier(unsigned char *TM,long lHeight,long lWidth);

void fft2d(double*TM,long lHeight,long lWidth,double*yRe,double*yIm);

void fft2d(unsigned char*TM,long lHeight,long lWidth,double*yRe,double*yIm,int DataType);
void normYi(double*yRe,double*yIm,unsigned char*Yi,int num,int DataType);
void ifft2d(unsigned char*TM,long lHeight,long lWidth,double*yRe,double*yIm,int DataType); //傅里叶变换返回unsigned char
void ifft2d_new(double*TM,long lHeight,long lWidth,double*yRe,double*yIm,int DataType); //傅里叶变换返回浮点值

void FdProc(double*meanfd,double*fdmark,int len,long height,long width,double minifreq);
void IGPI_DePeriodicNoise(unsigned char*IMG,unsigned char*IMG_OUT,long height,long width,double minifreq,int DataType);
double filter(int length,int i,double r1,double r2,double cutoff);
double filter2d(int height,int width,int i,int j,double r1,double r2,double cutoff);
void IGPI_Homomorphic(unsigned char*TM,long lHeight,long lWidth,double rl,double rh,double cutoff,int DataType);

//两幅图像在傅里叶频域里进行匹配，找出平移量
void FFT_Match(unsigned char *src1,unsigned char *src2,long lHeight,long lWidth,int &xout, int &yout);

//傅里叶正反变换的结果及显示
void FFT_all(unsigned char *image,unsigned char *fft_image,unsigned char *ifft_image,long lWidth,long lHeight);

//正常的傅里叶正变换
void FFT(unsigned char *image,double *fft_image,long lWidth,long lHeight);
void FFTUC(unsigned char *image, unsigned char *fft_image, long lWidth, long lHeight);
//傅里叶频谱搬移,对任意宽度，其中处理了长宽不为偶数的情况。
void fftshift(unsigned char *image,unsigned char *fftshift_image,long lWidth,long lHeight);

//傅里叶频谱搬移,对任意宽度，其中处理了长宽不为偶数的情况。
void fftshift_2(double *image,double *fftshift_image,long lWidth,long lHeight);

void fftshift_2uc(unsigned char *image, unsigned char *fftshift_image, long lWidth, long lHeight);
//计算长方环傅里叶功率谱的百分比
void Fourier_pu_distrbution(unsigned char *image,long lHeight,long lWidth,int N,double *pu_percent);

//整幅图像的积分图
void Image_integral(double *image,double *integral,int lwidth,int lheight);

//计算边缘的互功率谱，找出平移量,两个二值图像之间的移动
void Edge_Match(unsigned char *src1,int width1,int height1,unsigned char *src2,int width2,int height2,int &xout, int &yout);

#endif
