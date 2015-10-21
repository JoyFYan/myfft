#include "FT_PubFunction.h"



#define  maxPrimeFactor        65537//37
#define  maxPrimeFactorDiv2    (maxPrimeFactor+1)/2
#define  maxFactorCount        20

double  c3_1 = -1.5000000000000E+00;  /*  c3_1 = cos(2*pi/3)-1;          */
double  c3_2 =  8.6602540378444E-01;  /*  c3_2 = sin(2*pi/3);            */

double  u5   =  1.2566370614359E+00;  /*  u5   = 2*pi/5;                 */
double  c5_1 = -1.2500000000000E+00;  /*  c5_1 = (cos(u5)+cos(2*u5))/2-1;*/
double  c5_2 =  5.5901699437495E-01;  /*  c5_2 = (cos(u5)-cos(2*u5))/2;  */
double  c5_3 = -9.5105651629515E-01;  /*  c5_3 = -sin(u5);               */
double  c5_4 = -1.5388417685876E+00;  /*  c5_4 = -(sin(u5)+sin(2*u5));   */
double  c5_5 =  3.6327126400268E-01;  /*  c5_5 = (sin(u5)-sin(2*u5));    */
double  c8   =  7.0710678118655E-01;  /*  c8 = 1/sqrt(2);    */

double   pi;
int      groupOffset,dataOffset,blockOffset,adr;
int      groupNo,dataNo,blockNo,twNo;
double   omega, tw_re,tw_im;
double   twiddleRe[maxPrimeFactor], twiddleIm[maxPrimeFactor],
trigRe[maxPrimeFactor], trigIm[maxPrimeFactor],
zRe[maxPrimeFactor], zIm[maxPrimeFactor];
double   vRe[maxPrimeFactorDiv2], vIm[maxPrimeFactorDiv2];
double   wRe[maxPrimeFactorDiv2], wIm[maxPrimeFactorDiv2];

#define  PI  3.141592653

void factorize(int n, int *nFact, int fact[])
{
    int i,j,k;
    int nRadix;
    int radices[7];
    int factors[maxFactorCount];

    nRadix    =  6;  
    radices[1]=  2;
    radices[2]=  3;
    radices[3]=  4;
    radices[4]=  5;
    radices[5]=  8;
    radices[6]= 10;

    if (n==1)
    {
        j=1;
        factors[1]=1;
    }
    else j=0;
    i=nRadix;
    while ((n>1) && (i>0))
    {
      if ((n % radices[i]) == 0)
      {
        n=n / radices[i];
        j=j+1;
        factors[j]=radices[i];
      }
      else  i=i-1;
    }
    if (factors[j] == 2)   /*substitute factors 2*8 with 4*4 */
    {   
      i = j-1;
      while ((i>0) && (factors[i] != 8)) i--;
      if (i>0)
      {
        factors[j] = 4;
        factors[i] = 4;
      }
    }
    if (n>1)
    {
        for (k=2; k<sqrt(n)+1; k++)
            while ((n % k) == 0)
            {
                n=n / k;
                j=j+1;
                factors[j]=k;
            }
        if (n>1)
        {
            j=j+1;
            factors[j]=n;
        }
    }               
    for (i=1; i<=j; i++)         
    {
      fact[i] = factors[j-i+1];
    }
    *nFact=j;
}   /* factorize */

/****************************************************************************
  After N is factored the parameters that control the stages are generated.
  For each stage we have:
    sofar   : the product of the radices so far.
    actual  : the radix handled in this stage.
    remain  : the product of the remaining radices.
 ****************************************************************************/

void transTableSetup(int sofar[], int actual[], int remain[],
                     int *nFact,
                     int *nPoints)
{
    int i;

    factorize(*nPoints, nFact, actual);
    if (actual[1] > maxPrimeFactor)
    {
        printf("\nPrime factor of FFT length too large : %6d",actual[1]);
        printf("\nPlease modify the value of maxPrimeFactor in mixfft.c");
        exit(1);
    }
    remain[0]=*nPoints;
    sofar[1]=1;
    remain[1]=*nPoints / actual[1];
    for (i=2; i<=*nFact; i++)
    {
        sofar[i]=sofar[i-1]*actual[i-1];
        remain[i]=remain[i-1] / actual[i];
    }
}   /* transTableSetup */

/****************************************************************************
  The sequence y is the permuted input sequence x so that the following
  transformations can be performed in-place, and the final result is the
  normal order.
 ****************************************************************************/

void permute(int nPoint, int nFact,
             int fact[], int remain[],
             double xRe[], double xIm[],
             double yRe[], double yIm[])

{
    int i,j,k;
    int count[maxFactorCount]; 

    for (i=1; i<=nFact; i++) count[i]=0;
    k=0;
    for (i=0; i<=nPoint-2; i++)
    {
        yRe[i] = xRe[k];
        yIm[i] = xIm[k];
        j=1;
        k=k+remain[j];
        count[1] = count[1]+1;
        while (count[j] >= fact[j])
        {
            count[j]=0;
            k=k-remain[j-1]+remain[j+1];
            j=j+1;
            count[j]=count[j]+1;
        }
    }
    yRe[nPoint-1]=xRe[nPoint-1];
    yIm[nPoint-1]=xIm[nPoint-1];
}   /* permute */


/****************************************************************************
  Twiddle factor multiplications and transformations are performed on a
  group of data. The number of multiplications with 1 are reduced by skipping
  the twiddle multiplication of the first stage and of the first group of the
  following stages.
 ***************************************************************************/

void initTrig(int radix)
{
    int i;
    double w,xre,xim;

    w=2*pi/radix;
    trigRe[0]=1; trigIm[0]=0;
    xre=cos(w); 
    xim=-sin(w);
    trigRe[1]=xre; trigIm[1]=xim;
    for (i=2; i<radix; i++)
    {
        trigRe[i]=xre*trigRe[i-1] - xim*trigIm[i-1];
        trigIm[i]=xim*trigRe[i-1] + xre*trigIm[i-1];
    }
}   /* initTrig */

void fft_4(double aRe[], double aIm[])
{
    double  t1_re,t1_im, t2_re,t2_im;
    double  m2_re,m2_im, m3_re,m3_im;

    t1_re=aRe[0] + aRe[2]; t1_im=aIm[0] + aIm[2];
    t2_re=aRe[1] + aRe[3]; t2_im=aIm[1] + aIm[3];

    m2_re=aRe[0] - aRe[2]; m2_im=aIm[0] - aIm[2];
    m3_re=aIm[1] - aIm[3]; m3_im=aRe[3] - aRe[1];

    aRe[0]=t1_re + t2_re; aIm[0]=t1_im + t2_im;
    aRe[2]=t1_re - t2_re; aIm[2]=t1_im - t2_im;
    aRe[1]=m2_re + m3_re; aIm[1]=m2_im + m3_im;
    aRe[3]=m2_re - m3_re; aIm[3]=m2_im - m3_im;
}   /* fft_4 */


void fft_5(double aRe[], double aIm[])
{    
    double  t1_re,t1_im, t2_re,t2_im, t3_re,t3_im;
    double  t4_re,t4_im, t5_re,t5_im;
    double  m2_re,m2_im, m3_re,m3_im, m4_re,m4_im;
    double  m1_re,m1_im, m5_re,m5_im;
    double  s1_re,s1_im, s2_re,s2_im, s3_re,s3_im;
    double  s4_re,s4_im, s5_re,s5_im;

    t1_re=aRe[1] + aRe[4]; t1_im=aIm[1] + aIm[4];
    t2_re=aRe[2] + aRe[3]; t2_im=aIm[2] + aIm[3];
    t3_re=aRe[1] - aRe[4]; t3_im=aIm[1] - aIm[4];
    t4_re=aRe[3] - aRe[2]; t4_im=aIm[3] - aIm[2];
    t5_re=t1_re + t2_re; t5_im=t1_im + t2_im;
    aRe[0]=aRe[0] + t5_re; aIm[0]=aIm[0] + t5_im;
    m1_re=c5_1*t5_re; m1_im=c5_1*t5_im;
    m2_re=c5_2*(t1_re - t2_re); m2_im=c5_2*(t1_im - t2_im);

    m3_re=-c5_3*(t3_im + t4_im); m3_im=c5_3*(t3_re + t4_re);
    m4_re=-c5_4*t4_im; m4_im=c5_4*t4_re;
    m5_re=-c5_5*t3_im; m5_im=c5_5*t3_re;

    s3_re=m3_re - m4_re; s3_im=m3_im - m4_im;
    s5_re=m3_re + m5_re; s5_im=m3_im + m5_im;
    s1_re=aRe[0] + m1_re; s1_im=aIm[0] + m1_im;
    s2_re=s1_re + m2_re; s2_im=s1_im + m2_im;
    s4_re=s1_re - m2_re; s4_im=s1_im - m2_im;

    aRe[1]=s2_re + s3_re; aIm[1]=s2_im + s3_im;
    aRe[2]=s4_re + s5_re; aIm[2]=s4_im + s5_im;
    aRe[3]=s4_re - s5_re; aIm[3]=s4_im - s5_im;
    aRe[4]=s2_re - s3_re; aIm[4]=s2_im - s3_im;
}   /* fft_5 */

void fft_8()
{
    double  aRe[4], aIm[4], bRe[4], bIm[4], gem;

    aRe[0] = zRe[0];    bRe[0] = zRe[1];
    aRe[1] = zRe[2];    bRe[1] = zRe[3];
    aRe[2] = zRe[4];    bRe[2] = zRe[5];
    aRe[3] = zRe[6];    bRe[3] = zRe[7];

    aIm[0] = zIm[0];    bIm[0] = zIm[1];
    aIm[1] = zIm[2];    bIm[1] = zIm[3];
    aIm[2] = zIm[4];    bIm[2] = zIm[5];
    aIm[3] = zIm[6];    bIm[3] = zIm[7];

    fft_4(aRe, aIm); fft_4(bRe, bIm);

    gem    = c8*(bRe[1] + bIm[1]);
    bIm[1] = c8*(bIm[1] - bRe[1]);
    bRe[1] = gem;
    gem    = bIm[2];
    bIm[2] =-bRe[2];
    bRe[2] = gem;
    gem    = c8*(bIm[3] - bRe[3]);
    bIm[3] =-c8*(bRe[3] + bIm[3]);
    bRe[3] = gem;
    
    zRe[0] = aRe[0] + bRe[0]; zRe[4] = aRe[0] - bRe[0];
    zRe[1] = aRe[1] + bRe[1]; zRe[5] = aRe[1] - bRe[1];
    zRe[2] = aRe[2] + bRe[2]; zRe[6] = aRe[2] - bRe[2];
    zRe[3] = aRe[3] + bRe[3]; zRe[7] = aRe[3] - bRe[3];

    zIm[0] = aIm[0] + bIm[0]; zIm[4] = aIm[0] - bIm[0];
    zIm[1] = aIm[1] + bIm[1]; zIm[5] = aIm[1] - bIm[1];
    zIm[2] = aIm[2] + bIm[2]; zIm[6] = aIm[2] - bIm[2];
    zIm[3] = aIm[3] + bIm[3]; zIm[7] = aIm[3] - bIm[3];
}   /* fft_8 */

void fft_10()
{
    double  aRe[5], aIm[5], bRe[5], bIm[5];

    aRe[0] = zRe[0];    bRe[0] = zRe[5];
    aRe[1] = zRe[2];    bRe[1] = zRe[7];
    aRe[2] = zRe[4];    bRe[2] = zRe[9];
    aRe[3] = zRe[6];    bRe[3] = zRe[1];
    aRe[4] = zRe[8];    bRe[4] = zRe[3];

    aIm[0] = zIm[0];    bIm[0] = zIm[5];
    aIm[1] = zIm[2];    bIm[1] = zIm[7];
    aIm[2] = zIm[4];    bIm[2] = zIm[9];
    aIm[3] = zIm[6];    bIm[3] = zIm[1];
    aIm[4] = zIm[8];    bIm[4] = zIm[3];

    fft_5(aRe, aIm); fft_5(bRe, bIm);

    zRe[0] = aRe[0] + bRe[0]; zRe[5] = aRe[0] - bRe[0];
    zRe[6] = aRe[1] + bRe[1]; zRe[1] = aRe[1] - bRe[1];
    zRe[2] = aRe[2] + bRe[2]; zRe[7] = aRe[2] - bRe[2];
    zRe[8] = aRe[3] + bRe[3]; zRe[3] = aRe[3] - bRe[3];
    zRe[4] = aRe[4] + bRe[4]; zRe[9] = aRe[4] - bRe[4];

    zIm[0] = aIm[0] + bIm[0]; zIm[5] = aIm[0] - bIm[0];
    zIm[6] = aIm[1] + bIm[1]; zIm[1] = aIm[1] - bIm[1];
    zIm[2] = aIm[2] + bIm[2]; zIm[7] = aIm[2] - bIm[2];
    zIm[8] = aIm[3] + bIm[3]; zIm[3] = aIm[3] - bIm[3];
    zIm[4] = aIm[4] + bIm[4]; zIm[9] = aIm[4] - bIm[4];
}   /* fft_10 */

void fft_odd(int radix)
{
    double  rere, reim, imre, imim;
    int     i,j,k,n,max;

    n = radix;
    max = (n + 1)/2;
    for (j=1; j < max; j++)
    {
      vRe[j] = zRe[j] + zRe[n-j];
      vIm[j] = zIm[j] - zIm[n-j];
      wRe[j] = zRe[j] - zRe[n-j];
      wIm[j] = zIm[j] + zIm[n-j];
    }

    for (j=1; j < max; j++)
    {
        zRe[j]=zRe[0]; 
        zIm[j]=zIm[0];
        zRe[n-j]=zRe[0]; 
        zIm[n-j]=zIm[0];
        k=j;
        for (i=1; i < max; i++)
        {
            rere = trigRe[k] * vRe[i];
            imim = trigIm[k] * vIm[i];
            reim = trigRe[k] * wIm[i];
            imre = trigIm[k] * wRe[i];
            
            zRe[n-j] += rere + imim;
            zIm[n-j] += reim - imre;
            zRe[j]   += rere - imim;
            zIm[j]   += reim + imre;

            k = k + j;
            if (k >= n)  k = k - n;
        }
    }
    for (j=1; j < max; j++)
    {
        zRe[0]=zRe[0] + vRe[j]; 
        zIm[0]=zIm[0] + wIm[j];
    }
}   /* fft_odd */


void twiddleTransf(int sofarRadix, int radix, int remainRadix,
                    double yRe[], double yIm[])

{   /* twiddleTransf */ 
    double  cosw, sinw, gem;
    double  t1_re,t1_im, t2_re,t2_im, t3_re,t3_im;
    double  t4_re,t4_im, t5_re,t5_im;
    double  m2_re,m2_im, m3_re,m3_im, m4_re,m4_im;
    double  m1_re,m1_im, m5_re,m5_im;
    double  s1_re,s1_im, s2_re,s2_im, s3_re,s3_im;
    double  s4_re,s4_im, s5_re,s5_im;


    initTrig(radix);
    omega = 2*pi/(double)(sofarRadix*radix);
    cosw =  cos(omega);
    sinw = -sin(omega);
    tw_re = 1.0;
    tw_im = 0;
    dataOffset=0;
    groupOffset=dataOffset;
    adr=groupOffset;
    for (dataNo=0; dataNo<sofarRadix; dataNo++)
    {
        if (sofarRadix>1)
        {
            twiddleRe[0] = 1.0; 
            twiddleIm[0] = 0.0;
            twiddleRe[1] = tw_re;
            twiddleIm[1] = tw_im;
            for (twNo=2; twNo<radix; twNo++)
            {
                twiddleRe[twNo]=tw_re*twiddleRe[twNo-1]
                               - tw_im*twiddleIm[twNo-1];
                twiddleIm[twNo]=tw_im*twiddleRe[twNo-1]
                               + tw_re*twiddleIm[twNo-1];
            }
            gem   = cosw*tw_re - sinw*tw_im;
            tw_im = sinw*tw_re + cosw*tw_im;
            tw_re = gem;                      
        }
        for (groupNo=0; groupNo<remainRadix; groupNo++)
        {
            if ((sofarRadix>1) && (dataNo > 0))
            {
                zRe[0]=yRe[adr];
                zIm[0]=yIm[adr];
                blockNo=1;
                do {
                    adr = adr + sofarRadix;
                    zRe[blockNo]=  twiddleRe[blockNo] * yRe[adr]
                                 - twiddleIm[blockNo] * yIm[adr];
                    zIm[blockNo]=  twiddleRe[blockNo] * yIm[adr]
                                 + twiddleIm[blockNo] * yRe[adr]; 
                    
                    blockNo++;
                } while (blockNo < radix);
            }
            else
                for (blockNo=0; blockNo<radix; blockNo++)
                {
                   zRe[blockNo]=yRe[adr];
                   zIm[blockNo]=yIm[adr];
                   adr=adr+sofarRadix;
                }
            switch(radix) {
              case  2  : gem=zRe[0] + zRe[1];
                         zRe[1]=zRe[0] -  zRe[1]; zRe[0]=gem;
                         gem=zIm[0] + zIm[1];
                         zIm[1]=zIm[0] - zIm[1]; zIm[0]=gem;
                         break;
              case  3  : t1_re=zRe[1] + zRe[2]; t1_im=zIm[1] + zIm[2];
                         zRe[0]=zRe[0] + t1_re; zIm[0]=zIm[0] + t1_im;
                         m1_re=c3_1*t1_re; m1_im=c3_1*t1_im;
                         m2_re=c3_2*(zIm[1] - zIm[2]); 
                         m2_im=c3_2*(zRe[2] -  zRe[1]);
                         s1_re=zRe[0] + m1_re; s1_im=zIm[0] + m1_im;
                         zRe[1]=s1_re + m2_re; zIm[1]=s1_im + m2_im;
                         zRe[2]=s1_re - m2_re; zIm[2]=s1_im - m2_im;
                         break;
              case  4  : t1_re=zRe[0] + zRe[2]; t1_im=zIm[0] + zIm[2];
                         t2_re=zRe[1] + zRe[3]; t2_im=zIm[1] + zIm[3];

                         m2_re=zRe[0] - zRe[2]; m2_im=zIm[0] - zIm[2];
                         m3_re=zIm[1] - zIm[3]; m3_im=zRe[3] - zRe[1];

                         zRe[0]=t1_re + t2_re; zIm[0]=t1_im + t2_im;
                         zRe[2]=t1_re - t2_re; zIm[2]=t1_im - t2_im;
                         zRe[1]=m2_re + m3_re; zIm[1]=m2_im + m3_im;
                         zRe[3]=m2_re - m3_re; zIm[3]=m2_im - m3_im;
                         break;
              case  5  : t1_re=zRe[1] + zRe[4]; t1_im=zIm[1] + zIm[4];
                         t2_re=zRe[2] + zRe[3]; t2_im=zIm[2] + zIm[3];
                         t3_re=zRe[1] - zRe[4]; t3_im=zIm[1] - zIm[4];
                         t4_re=zRe[3] - zRe[2]; t4_im=zIm[3] - zIm[2];
                         t5_re=t1_re + t2_re; t5_im=t1_im + t2_im;
                         zRe[0]=zRe[0] + t5_re; zIm[0]=zIm[0] + t5_im;
                         m1_re=c5_1*t5_re; m1_im=c5_1*t5_im;
                         m2_re=c5_2*(t1_re - t2_re); 
                         m2_im=c5_2*(t1_im - t2_im);

                         m3_re=-c5_3*(t3_im + t4_im); 
                         m3_im=c5_3*(t3_re + t4_re);
                         m4_re=-c5_4*t4_im; m4_im=c5_4*t4_re;
                         m5_re=-c5_5*t3_im; m5_im=c5_5*t3_re;

                         s3_re=m3_re - m4_re; s3_im=m3_im - m4_im;
                         s5_re=m3_re + m5_re; s5_im=m3_im + m5_im;
                         s1_re=zRe[0] + m1_re; s1_im=zIm[0] + m1_im;
                         s2_re=s1_re + m2_re; s2_im=s1_im + m2_im;
                         s4_re=s1_re - m2_re; s4_im=s1_im - m2_im;

                         zRe[1]=s2_re + s3_re; zIm[1]=s2_im + s3_im;
                         zRe[2]=s4_re + s5_re; zIm[2]=s4_im + s5_im;
                         zRe[3]=s4_re - s5_re; zIm[3]=s4_im - s5_im;
                         zRe[4]=s2_re - s3_re; zIm[4]=s2_im - s3_im;
                         break;
              case  8  : fft_8(); break;
              case 10  : fft_10(); break;
              default  : fft_odd(radix); break;
            }
            adr=groupOffset;
            for (blockNo=0; blockNo<radix; blockNo++)
            {
                yRe[adr]=zRe[blockNo]; yIm[adr]=zIm[blockNo];
                adr=adr+sofarRadix;
            }
            groupOffset=groupOffset+sofarRadix*radix;
            adr=groupOffset;
        }
        dataOffset=dataOffset+1;
        groupOffset=dataOffset;
        adr=groupOffset;
    }
}   /* twiddleTransf */

void fft(int n, double xRe[], double xIm[],
                double yRe[], double yIm[])
{
    int   sofarRadix[maxFactorCount], 
          actualRadix[maxFactorCount], 
          remainRadix[maxFactorCount];
    int   nFactor;
    int   count;

    pi = 4*atan(1);    

    transTableSetup(sofarRadix, actualRadix, remainRadix, &nFactor, &n);
    permute(n, nFactor, actualRadix, remainRadix, xRe, xIm, yRe, yIm);

    for (count=1; count<=nFactor; count++)
      twiddleTransf(sofarRadix[count], actualRadix[count], remainRadix[count], 
                    yRe, yIm);
}   /* fft */

void ifft(int n,double yRe[],double yIm[],double xRe[],double xIm[])
{
	// 付立叶变换点数
	int	count=n;
	// 循环变量
	int		i;
	// 分配运算所需存储器
    double* yRetmp,*yImtmp;
	yRetmp=(double*)malloc(n*sizeof(double));
	yImtmp=(double*)malloc(n*sizeof(double));
	// 将频域点写入X
	memcpy(yRetmp, yRe, sizeof(double)*n);
	memcpy(yImtmp, yIm, sizeof(double)*n);
	// 求共轭
	for(i = 0; i <count; i++)
	{
		yImtmp[i] = -yIm[i];
	}
	// 调用快速付立叶变换
	fft(n,xRe,xIm,yRetmp,yImtmp);
	// 求时域点的共轭
	for(i = 0; i < count; i++)
	{
		xRe[i]/=count;xIm[i]/=count;
	
	}
	// 释放内存
	free(yRetmp);free(yImtmp);
}

void Fourier(unsigned char*TM,long lHeight,long lWidth)
{
	// 中间变量
	double	dTemp=0.0,tmp=0.0;
	// 循环变量
	long	i=0,j=0,h=lHeight,w=lWidth;
	double *xRe,*xIm,*yRe,*yIm;
	xRe=(double*)malloc(h*w*sizeof(double));
	xIm=(double*)malloc(h*w*sizeof(double));
	yRe=(double*)malloc(h*w*sizeof(double));
	yIm=(double*)malloc(h*w*sizeof(double));
	for(i=0;i<h;i++){
		for(j=0;j<w;j++){
			// 指向DIB第i行，第j个象素的指针
			tmp = double(TM[i*lWidth+j]);
			// 给时域赋值
			xRe[j+w*i]=tmp;xIm[j+w*i]=0.0;
		}
	}
	for(i=0;i<h;i++){
		// 对y方向进行快速付立叶变换
		fft(w,&xRe[w*i],&xIm[w*i],&yRe[w*i],&yIm[w*i]);
	}
	// 保存变换结果
	for(i=0;i<h;i++){
		for(j=0;j<w;j++){
			xRe[i+h*j]=yRe[j+w*i];
			xIm[i+h*j]=yIm[j+w*i];
		}
	}
	for(i=0;i<w;i++){
		// 对x方向进行快速付立叶变换
		fft(h,&xRe[h*i],&xIm[h*i],&yRe[h*i],&yIm[h*i]);
	}
	int d=1;
	for(i=0;i<h;i++){
		for(j=0;j<w;j++){
			// 计算频谱
			if(d==1){
				dTemp=sqrt(yRe[j*h+i]*yRe[j*h+i]+ 
							 yIm[j*h+i]*yIm[j*h+i])/100;
				if (dTemp>255) dTemp=255;
				TM[lWidth*(i<h/2?i+h/2:i-h/2)+(j<w/2?j+w/2:j-w/2)]=(BYTE)(dTemp);
				//TM[lWidth*i+j]=(BYTE)(dTemp);
			}
			else{
				dTemp=sqrt(yRe[j+i*w]*yRe[j+i*w]+ 
							 yIm[j+i*w]*yIm[j+i*w])/100;
				if (dTemp>255) dTemp=255;
				TM[lWidth*i+j]=(BYTE)(dTemp);
			}
			// 判断是否超过255
			// 对于超过的，直接设置为255
			// 指向DIB第(i<h/2 ? i+h/2 : i-h/2)行，第(j<w/2 ? j+w/2 : j-w/2)个象素的指针
			// 此处不直接取i和j，是为了将变换后的原点移到中心
			// 更新源图像
			
			//TM[lWidth*(i<h/2?i+h/2:i-h/2)+(j<w/2?j+w/2:j-w/2)]=(BYTE)(dTemp);
		}
	}
	// 删除临时变量
	free(xRe);free(xIm);free(yRe);free(yIm);
}

void fft2d(double*TM,long lHeight,long lWidth,double*yRe,double*yIm)
{
	// 中间变量
	double	tmp=0.0;
	// 循环变量
	LONG	i=0,j=0,h=lHeight,w=lWidth;
	double *xRe,*xIm;
	xRe=(double*)malloc(h*w*sizeof(double));
	xIm=(double*)malloc(h*w*sizeof(double));
	for(i=0;i<h;i++){
		for(j=0;j<w;j++){
			// 指向DIB第i行，第j个象素的指针
			tmp = TM[i*lWidth+j];
			// 给时域赋值
			xRe[j+w*i]=tmp;xIm[j+w*i]=0.0;
		}
	}
	for(i=0;i<h;i++){
		// 对y方向进行快速付立叶变换
		fft(w,&xRe[w*i],&xIm[w*i],&yRe[w*i],&yIm[w*i]);
	}
	// 保存变换结果
	for(i=0;i<h;i++){
		for(j=0;j<w;j++){
			xRe[i+h*j]=yRe[j+w*i];
			xIm[i+h*j]=yIm[j+w*i];
		}
	}
	for(i=0;i<w;i++){
		// 对x方向进行快速付立叶变换
		fft(h,&xRe[h*i],&xIm[h*i],&yRe[h*i],&yIm[h*i]);
	}
	free(xRe);free(xIm);

}

void fft2d(unsigned char*TM,long lHeight,long lWidth,double*yRe,double*yIm,int DataType=8)
{
	// 中间变量
	int	tmp=0;
	// 循环变量
	LONG	i=0,j=0,h=lHeight,w=lWidth;
	double *xRe,*xIm;
	xRe=(double*)malloc(h*w*sizeof(double));
	xIm=(double*)malloc(h*w*sizeof(double));
	for(i=0;i<h;i++){
		for(j=0;j<w;j++){
			// 指向DIB第i行，第j个象素的指针
			// 给时域赋值
			if(DataType==8){
				xRe[j+w*i]=TM[i*lWidth+j];xIm[j+w*i]=0.0;
			}
			if(DataType==16){
				tmp=TM[i*w+j*2+1];
				tmp=(tmp<<8)+TM[i*w+j*2];
				xRe[i*w+j]=tmp;xIm[j+w*i]=0.;
			}
		}
	}
	for(i=0;i<h;i++){
		// 对y方向进行快速付立叶变换
		fft(w,&xRe[w*i],&xIm[w*i],&yRe[w*i],&yIm[w*i]);
	}
	// 保存变换结果
	for(i=0;i<h;i++){
		for(j=0;j<w;j++){
			xRe[i+h*j]=yRe[j+w*i];
			xIm[i+h*j]=yIm[j+w*i];
		}
	}
	for(i=0;i<w;i++){
		// 对x方向进行快速付立叶变换
		fft(h,&xRe[h*i],&xIm[h*i],&yRe[h*i],&yIm[h*i]);
	}
	free(xRe);free(xIm);
}

void normYi(double*yRe,double*yIm,unsigned char*Yi,int num,int DataType=8)
{
	int i;double tmp=0.0;double max=0;double min=1.0e90;
	double*ytmp=(double*)malloc(num*sizeof(double));
	for(i=0;i<num;i++)
	{
		tmp=sqrt(yRe[i]*yRe[i]+yIm[i]*yIm[i])/num;
		if(tmp>max) max=tmp;
		if(tmp<min) min=tmp;
	}
	double interval=max-min;
	for(i=0;i<num;i++)
	{
		tmp=sqrt(yRe[i]*yRe[i]+yIm[i]*yIm[i])/num;
		if(DataType==8)
			Yi[i]=BYTE((tmp-min)*255/interval);
		if(DataType==16)
		{
			tmp=(tmp-min)*65535/interval;
			Yi[2*i]=BYTE(int(tmp) & 255);
			Yi[2*i+1]=BYTE((int(tmp)>>8) & 255);
		}
	}
	free(ytmp);
}

void ifft2d(unsigned char*TM,long lHeight,long lWidth,double*yRe,double*yIm,int DataType=8)
{
	// 中间变量
	double	dTemp=0.0,tmp=0.0;
	// 循环变量
	LONG	i=0,j=0,h=lHeight,w=lWidth;
	double *xRe,*xIm;
	xRe=(double*)malloc(h*w*sizeof(double));
	xIm=(double*)malloc(h*w*sizeof(double));
	for(i=0;i<h;i++){
		for(j=0;j<w;j++){
			yIm[j+w*i]=-yIm[w*i+j];
		}
	}
	for(i=0;i<w;i++){
		// 对y方向进行快速付立叶变换
		fft(h,&yRe[h*i],&yIm[h*i],&xRe[h*i],&xIm[h*i]);
	}
	// 保存变换结果
	for(i=0;i<h;i++){
		for(j=0;j<w;j++){
			yRe[j+w*i]=xRe[i+h*j];
			yIm[j+w*i]=xIm[i+h*j];
		}
	}
	for(i=0;i<h;i++){
		// 对x方向进行快速付立叶变换
		fft(w,&yRe[w*i],&yIm[w*i],&xRe[w*i],&xIm[w*i]);
	}
	//对TM进行变换，变换到0-255，便于显示
	normYi(xRe,xIm,TM,h*w,DataType);
	// 删除临时变量
	free(xRe);free(xIm);
}

void ifft2d_new(double*TM,long lHeight,long lWidth,double*yRe,double*yIm,int DataType=8)
{
	// 中间变量
	double	dTemp=0.0,tmp=0.0;
	// 循环变量
	LONG	i=0,j=0,h=lHeight,w=lWidth;
	double *xRe,*xIm;
	xRe=(double*)malloc(h*w*sizeof(double));
	xIm=(double*)malloc(h*w*sizeof(double));
	for(i=0;i<h;i++)
	{
		for(j=0;j<w;j++)
		{
			yIm[j+w*i]=-yIm[w*i+j];
		}
	}
	for(i=0;i<w;i++){
		// 对y方向进行快速付立叶变换
		fft(h,&yRe[h*i],&yIm[h*i],&xRe[h*i],&xIm[h*i]);
	}
	// 保存变换结果
	for(i=0;i<h;i++){
		for(j=0;j<w;j++){
			yRe[j+w*i]=xRe[i+h*j];
			yIm[j+w*i]=xIm[i+h*j];
		}
	}
	for(i=0;i<h;i++){
		// 对x方向进行快速付立叶变换
		fft(w,&yRe[w*i],&yIm[w*i],&xRe[w*i],&xIm[w*i]);
	}


	for(i=0;i<h*w;i++)
	{
		TM[i]=sqrt(yRe[i]*yRe[i]+yIm[i]*yIm[i]);
	}

	// 删除临时变量
	free(xRe);free(xIm);
}



void FdProc(double*meanfd,double*fdmark,int len,long height,long width,double minifreq)
{
	int i,j,m,n;
	int x=3,y=3;
	int xx,yy;
	float kw=float(len-1)/float(width-1);
	float kh=float(len-1)/float(height-1);
	double*tmp=(double*)malloc(len*len*sizeof(double));
	double mean=0.0,var=0.,k=1;
	int ll=(int)(len*minifreq/(2*PI));
	for(i=0;i<len;i++){
		for(j=0;j<len;j++){
			mean=0;
			for(m=-y;m<=y;m++){
				for(n=-x;n<=x;n++){
					xx=i+m;yy=j+n;
					if(xx<0) xx=xx+x;
					if(xx>len-1) xx=2*len-xx-2;
					if(yy<0) yy=yy+y;
					if(yy>len-1) yy=2*len-yy-2;
					mean+=meanfd[yy*len+xx];
				}
			}
			mean/=(2*x+1)*(2*y+1);
			tmp[i+j*len]=meanfd[i+j*len]-mean;
			if(tmp[i+j*len]<0) tmp[i+j*len]=0;
		}
	}
	for(i=0;i<len*len;i++) meanfd[i]=tmp[i];
	int count=0;
	for(i=0;i<len;i++){
		for(j=0;j<len;j++){
			if((i<ll && j<ll)||(i>len-ll && j<ll)||
				(i<ll && j>len-ll)||(i>len-ll && j>len-ll));
			else{
				mean+=meanfd[i+j*len];
				var+=meanfd[i+j*len]*meanfd[i+j*len];count++;
			}
		}
	}
	mean/=count;var/=count;
	var-=mean*mean;var=sqrt(var);
	for(i=0;i<len;i++){
		for(j=0;j<len;j++){
			if((i<ll && j<ll)||(i>len-ll && j<ll)||
				(i<ll && j>len-ll)||(i>len-ll && j>len-ll)){
				meanfd[i+j*len]=1;
			}
			else{
				if(meanfd[i+j*len]>mean+k*var)
					meanfd[i+j*len]=0;
				else meanfd[i+j*len]=1;
			}
		}
	}
//////////////////////////////////////////////////////*/
	for(i=0;i<height;i++){
		for(j=0;j<width;j++){
			x=(int)(j*kw);
			y=(int)(i*kh);
			if(x-int(x)>0.5) xx=int(x)+1;
			else xx=int(x);
			if(y-int(y)>0.5) yy=int(y)+1;
			else yy=int(y);
			fdmark[i+j*height]=meanfd[yy*len+xx];
		}
	}
	free(tmp);
}

void IGPI_DePeriodicNoise(unsigned char*IMG,unsigned char*IMG_OUT,long height,long width,double minifreq,int DataType=8)
{
	int len=128;int cn=0;
	double*TM=(double*)malloc(len*len*sizeof(double));
	double*meanfd=(double*)malloc(len*len*sizeof(double));
	double*fdmark=(double*)malloc(height*width*sizeof(double));
	double*yRe=(double*)malloc(len*len*sizeof(double));
	double*yIm=(double*)malloc(len*len*sizeof(double));
	double*R=(double*)malloc(height*width*sizeof(double));
	double*I=(double*)malloc(height*width*sizeof(double));
	int i,j,m,n;int h=int(height/len)*len;int w=int(width/len)*len;
	for(i=0;i<len*len;i++) meanfd[i]=0;
	for(i=0;i<h;i+=len){
		for(j=0;j<w;j+=len){
			for(m=0;m<len;m++){
				for(n=0;n<len;n++)
					if(DataType==8)
						TM[m*len+n]=IMG[(i+m)*width+j+n];
					if(DataType==16){
						int t=IMG[i*w+j*2+1];
						t=(t<<8)+IMG[i*w+j*2];
						TM[i*w+j]=double(t);
					}
			}
			fft2d(TM,len,len,yRe,yIm);
			for(m=0;m<len;m++)
				for(n=0;n<len;n++)
					meanfd[m+n*len]+=yRe[n+m*len]*yRe[n+m*len]+yIm[n+m*len]*yIm[n+m*len];
			cn++;
		}
	}
	for(i=0;i<len*len;i++) meanfd[i]=log(meanfd[i])/cn;
	for(i=0;i<height*width;i++) fdmark[i]=0.;
	meanfd[0]=0;meanfd[len*len-1]=0;meanfd[len-1]=0;meanfd[(len-1)*len]=0;
	////////////////////////////////////////////////
	FdProc(meanfd,fdmark,len,height,width,minifreq);
	fft2d(IMG,height,width,R,I,DataType);
	for(i=0;i<height*width;i++){
		R[i]*=fdmark[i];I[i]*=fdmark[i];
	}

	ifft2d(IMG_OUT,height,width,R,I,DataType);
	free(yRe);free(yIm);free(TM);free(meanfd);free(fdmark);free(R);free(I);
}

double filter(int length,int i,double r1,double r2,double cutoff)
{
	int ii=length-i;
	double angle1=PI*i/length;double a=0.0;
    double angle=PI*PI*i/(length*cutoff);
	if(angle1<cutoff)
        a=(r1+r2)/2-(r1-r2)*sin(angle-PI/2)/2;
    else
        a=r2;
	return a;
}

double filter2d(int height,int width,int i,int j,double r1,double r2,double cutoff)
{
	if(i==0 && j==0) return 1.05;
	int htmp=height/2;int wtmp=width/2;
	int itmp=0,jtmp=0;double ah=0.0,aw=0.0;
	if(htmp*2==height){
		itmp=(i<=htmp-1?i:height-i);
		ah=filter(htmp,itmp,r1,r2,cutoff);
	}
	else{
		itmp=(i<=htmp+1?i:height-i);
		ah=filter(htmp+1,itmp,r1,r2,cutoff);
	}
	if(wtmp*2==width){
		jtmp=(j<=wtmp-1?j:width-j);
		aw=filter(wtmp,jtmp,r1,r2,cutoff);
	}
	else{
		jtmp=(j<=wtmp+1?j:width-j);
		aw=filter(wtmp+1,jtmp,r1,r2,cutoff);
	}
	return sqrt(aw*ah);
}

void IGPI_Homomorphic(unsigned char*TM,long lHeight,long lWidth,double rl,double rh,double cutoff,int DataType=8)
{
	// 中间变量
	double	dTemp=0.0,tmp=0.0;
	// 循环变量
	LONG	i=0,j=0,h=lHeight,w=lWidth;
	double *yRe,*yIm,*TMLOG;int t=0;
	yRe=(double*)malloc(h*w*sizeof(double));
	yIm=(double*)malloc(h*w*sizeof(double));
	TMLOG=(double*)malloc(h*w*sizeof(double));
	for(i=0;i<h;i++){
		for(j=0;j<w;j++){
			if(DataType==8)	TMLOG[i*w+j]=log(double(TM[i*w+j])+1.);
			if(DataType==16){
				t=TM[i*w+j*2+1];
				t=(t<<8)+TM[i*w+j*2];
				TMLOG[i*w+j]=log(double(t+1));
			}
		}
	}
	fft2d(TMLOG,h,w,yRe,yIm);
	for(i=0;i<h;i++){
		for(j=0;j<w;j++){
			dTemp=filter2d(h,w,i,j,rl,rh,cutoff);
			yRe[i+j*h]*=dTemp;
			yIm[i+j*h]*=dTemp;
		}
	}
	ifft2d(TM,h,w,yRe,yIm,DataType);
	free(yRe);free(yIm);free(TMLOG);

}

/************************************************************************/
/* 函数名称:   
/*         FFT_all()            
/* 函数参数:  
/*   unsigned char *image          图像矩阵
/*   unsigned char *fft_image      图像傅里叶正变换结果 
/*   unsigned char *ifft_image     图像傅里叶反变换结果
/*   long lHeight                  图像宽度
/*   long lWidth                   图像高度
/* 说明：
/*    函数采用任意基的快速傅里叶变换，进行快速傅里叶正反变换并显示，
/*    注意这里的频谱并未搬移到中心        
/************************************************************************/
void FFT_all(unsigned char *image, unsigned char *fft_image,unsigned char *ifft_image,long lWidth,long lHeight)
{
   int lineByte = (lWidth+3)/4*4;

    fftw_complex *in, *out1,*out2;
    fftw_plan p,q;
    int i=0, j=0;

    double * temp = new double[lHeight*lWidth];
    memset(temp,0,lHeight*lWidth*sizeof(double));

    in   = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) *lHeight*lWidth );    //存放原始数据
    out1 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) *lHeight*lWidth );    //存放傅里叶正变换结果
    out2 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) *lHeight*lWidth );    //存放傅里叶反变换结果
	//float *modle = new float[lHeight*lWidth];
    int k = 0;
    for(i=0;i<lHeight;i++)
    {
        for (j=0;j<lWidth;j++,k++)
        {
            int k1=lineByte*i+j;  
            if((i+j)%2!=0)in[k][0] = (double)image[k1]*(-1);     //实部    源图像可能不符合宽度为4
			else in[k][0] = (double)image[k1];
			in[k][1] = 0.0;                   //虚部
        }
    }
      
    p=fftw_plan_dft_2d( lHeight, lWidth, in  , out1, FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_execute(p);
	k = 0;
	for (i = 0;i<lHeight;i++)
	{
		for (j = 0;j<lWidth;j++, k++)
		{
			//int k1 = lineByte*i + j;
			if ((sqrt((i + 1 - lHeight / 2) * (1 + i - lHeight / 2) + (j + 1 - lWidth / 2) * (j + 1 - lWidth / 2))) > 6)
			{
				out1[k][1] = 0;
				out1[k][0] = 0;
			}
		}
	}

	//for (i = 0;i < lHeight;i++)//构建滤波模板，点乘
	//{
	//	for (j = 0;j < lWidth;j++, k++)
	//	{
	//		
	//		//if ((sqrt((i + 1 - lHeight / 2) * (1 + i - lHeight / 2) + (j + 1 - lWidth / 2) * (j + 1 - lWidth / 2))) < 10)
	//		if (i+j<10)
	//		{
	//			out1[k][1] = 0;
	//			out1[k][0] = 0;
	//		}
	//		//D = sqrt((m - M) ^ 2 + (n - N) ^ 2);
	//		//H(m, n) = exp((-D ^ 2) / (2 * (50) ^ 2));
	//		//else modle[i * lHeight + j] = 0;

	//	}
	//	//printf("\n");
	//}





    q=fftw_plan_dft_2d( lHeight, lWidth, out1, out2, FFTW_BACKWARD, FFTW_ESTIMATE);
                     
    fftw_execute( q );

    //求fft变换后的频谱，并对其取log10，求出最大值
    double _max = -1;
    double _min = 10000;
    for(i=0;i<lHeight*lWidth;i++)
    {
        temp[i] = log10(sqrt(out1[i][0]*out1[i][0] + out1[i][1]*out1[i][1])+1);//之前log内没有加1
        if (_max<temp[i])
        {
            _max = temp[i];
        }
        if (_min>temp[i])
        {
            _min = temp[i];
        }
    }
    
    double bili = _max-_min;

    k = 0;
    int iamge_size = lWidth*lHeight;
    for(i=0;i<lHeight;i++)
    {
        for (j=0;j<lWidth;j++,k++)
        {
            //频谱去log10后归一化
            int k1=lineByte*i+j;  
             fft_image[k1]  = BYTE(fabs((temp[k]-_min)*255/bili));   //这个只是为了显示好看而已
            ifft_image[k1]  = BYTE(fabs(out2[k][0]/iamge_size));//fabs之后才正确显示
        }
    }
    
    fftw_destroy_plan(p);
    fftw_destroy_plan(q);
    fftw_free(in);
    fftw_free(out1);
    fftw_free(out2); 
}


/************************************************************************/
/* 函数名称:   
/*         FFT()            
/* 函数参数:  
/*   unsigned char *image          图像矩阵image为一维数组
/*   unsigned char *fft_image      图像傅里叶正变换结果 
/*   long lHeight                  图像宽度
/*   long lWidth                   图像高度
/* 说明：
/*    函数采用任意基的快速傅里叶变换，进行快速傅里叶正变   
/************************************************************************/
void FFT(unsigned char *image, double *fft_image,long lWidth,long lHeight)//image为一维数组
{
   int lineByte = (lWidth+3)/4*4;

    fftw_complex *in, *out;
    fftw_plan p;
    int i=0, j=0;

    in   = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) *lHeight*lWidth );    //存放原始数据
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) *lHeight*lWidth );    //存放傅里叶正变换结果
      
    int k = 0;
    for(i=0;i<lHeight;i++)
    {
        for (j=0;j<lWidth;j++,k++)
        {
            int k1=lineByte*i+j;  
            in[k][0] = (double)image[k1];     //实部    源图像可能不符合宽度为4
            in[k][1] = 0.0;                   //虚部
        }
    }
      
    p=fftw_plan_dft_2d( lHeight, lWidth, in  , out, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute( p );                   

    //求fft变换后的频谱
    k = 0;
    for(i=0;i<lHeight;i++)
    {
        for (j=0;j<lWidth;j++,k++)
        {
            int k1=lineByte*i+j;  
             fft_image[k1]  = sqrt(out[k][0]*out[k][0] + out[k][1]*out[k][1]);   //这个只是为了显示好看而已
          }
		k = 0;//
    }
    
    fftw_destroy_plan(p);
    fftw_free(in);
    fftw_free(out);
}


/************************************************************************/
/* 函数名称:   
/*         fftshift()            
/* 函数参数:  
/*   unsigned char *image          图像矩阵
/*   unsigned char *fftshift_image 频谱搬移后图像
/*   long lHeight                  图像宽度
/*   long lWidth                   图像高度
/* 说明：
/*    傅里叶频谱搬移方法，对角线和反对角线进行交换，长和宽不要求为偶数      
/************************************************************************/
void fftshift(unsigned char *image,unsigned char *fftshift_image,long lWidth,long lHeight)
{
    int linebyte = (lWidth+3)/4*4;

    int x,y;
    int half_height_ceil  = (int)ceil(lHeight/2.0);      //   ceil(5/2) =3
    int half_height_floor = (int)floor(lHeight/2.0);     //   floor(5/2)=2  
    int half_width_ceil   = (int)ceil(lWidth/2.0);
    int half_width_floor  = (int)floor(lWidth/2.0);

        
    //第一象限值换到第三象限
    for(y=0;y<half_height_floor;y++)
    {
        for (x=0;x<half_width_floor;x++)
        {
            int  k = linebyte*y+x;
            int  k1 =linebyte*(y+half_height_ceil)+x+half_width_ceil;
            fftshift_image [k] = image[k1];
        }
    }
    //第四象限值换到第二象限
    for(y=half_height_floor;y<lHeight;y++)
    {
        for (x=0;x<half_width_floor;x++)
        {
            int  k = linebyte*y+x;
            int  k1 =linebyte*(y-half_height_floor)+x+half_width_ceil;
            fftshift_image [k] = image[k1];
        }
    }

    //第二象限值换到第四象限
    for(y=0;y<half_height_floor;y++)
    {
        for (x=half_width_floor;x<lWidth;x++)
        {
            int  k = linebyte*y+x;
            int  k1 =linebyte*(y+half_height_ceil)+x-half_width_floor;
            fftshift_image [k] = image[k1];
        }
    }

    //第三象限值换到第一象限
    for(y=half_height_floor;y<lHeight;y++)
    {
        for (x=half_width_floor;x<lWidth;x++)
        {
            int  k = linebyte*y+x;
            int  k1 =linebyte*(y-half_height_floor)+x-half_width_floor;
            fftshift_image [k] = image[k1];
        }
    }

//     FILE * fp =fopen("fftshfit.txt","w");
//     for(y=0;y<lHeight;y++)
//     {
//         for (x=0;x<lWidth;x++)
//         {
//             int k = linebyte *y+x;
//             fprintf(fp,"%d  ",fftshift_image[k]);
//         }
//         fprintf(fp,"\n");
//     }
//     fclose(fp);

}

/************************************************************************/
/* 函数名称:   
/*         fftshift_2()            
/* 函数参数:  
/*   double *image          图像矩阵
/*   double *fftshift_image 频谱搬移后图像
/*   long   lWidth          图像宽度
/*   long   lHeight         图像高度
/* 说明：
/*    傅里叶频谱搬移方法，对角线和反对角线进行交换,其中长和宽都为偶数      
/************************************************************************/
void fftshift_2(double *image,double *fftshift_image,long lWidth,long lHeight)
{
    int linebyte = (lWidth+3)/4*4;

    int x,y;
    int half_height_floor = lHeight/2;     //   floor(5/2)=2  
    int half_width_floor  = lWidth/2;

        
    //第一象限和第三象限交换
    for(y=0;y<half_height_floor;y++)
    {
        for (x=0;x<half_width_floor;x++)
        {
            int  k = linebyte*y+x;
            int  k1 =linebyte*(y+half_height_floor)+x+half_width_floor;
            fftshift_image [k] = image[k1];
            fftshift_image [k1] = image[k];
        }
    }
    //第四象限和第二象限交换
    for(y=half_height_floor;y<lHeight;y++)
    {
        for (x=0;x<half_width_floor;x++)
        {
            int  k = linebyte*y+x;
            int  k1 =linebyte*(y-half_height_floor)+x+half_width_floor;
            fftshift_image [k]  = image[k1];
            fftshift_image [k1] = image[k];
        }
    }

}

void fftshift_2uc(UCHAR *image, UCHAR *fftshift_image, long lWidth, long lHeight)
{
	int linebyte = (lWidth + 3) / 4 * 4;

	int x, y;
	int half_height_floor = lHeight / 2;     //   floor(5/2)=2  
	int half_width_floor = lWidth / 2;


	//第一象限和第三象限交换
	for (y = 0;y<half_height_floor;y++)
	{
		for (x = 0;x<half_width_floor;x++)
		{
			int  k = linebyte*y + x;
			int  k1 = linebyte*(y + half_height_floor) + x + half_width_floor;
			fftshift_image[k] = image[k1];
			fftshift_image[k1] = image[k];
		}
	}
	//第四象限和第二象限交换
	for (y = half_height_floor;y<lHeight;y++)
	{
		for (x = 0;x<half_width_floor;x++)
		{
			int  k = linebyte*y + x;
			int  k1 = linebyte*(y - half_height_floor) + x + half_width_floor;
			fftshift_image[k] = image[k1];
			fftshift_image[k1] = image[k];
		}
	}

}

/***********************************************************************/
/* 函数名称:   
/*         FFT_Match()              
/* 函数参数:  
/*   unsigned char*src1     图像矩阵1
/*   unsigned char*src2     图像矩阵2
/*   long lHeight           图像高度  
/*   long lWidth            图像宽度
/* 说明：
/*    该函数实现频域自相关计算两幅图像的平移变量dx,dy
/*    函数采用任意基的快速傅里叶变换，实现任意图像尺寸的快速匹配       
/************************************************************************/

void  FFT_Match(unsigned char *src1,unsigned char *src2,long lHeight,long lWidth,int &xout, int &yout)
{
    fftw_complex *in1,*in2;          //导入两幅图像
    fftw_complex *out1,*out2;        //傅里叶正变换输出两个变换结果 
    fftw_complex *out3;              //一次反变换
    fftw_plan p1,p2,q;

    int i=0,j=0;

    in1   = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) *lHeight*lWidth );    //存放原始数据
    in2   = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) *lHeight*lWidth );    //存放原始数据

    out1 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) *lHeight*lWidth );    //存放傅里叶正变换结果
    out2 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) *lHeight*lWidth );    //存放傅里叶反变换结果
    out3 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) *lHeight*lWidth );    //存放傅里叶反变换结果
   
    //两幅傅里叶图像赋值
    for(i=0;i<lHeight*lWidth;i++)
    {
        in1[i][0] = (double)src1[i];        //实部    
        in1[i][1] = 0.0;                   //虚部
        in2[i][0] = (double)src2[i];        //实部    
        in2[i][1] = 0.0;                   //虚部
    }
    
    //两幅图像进行傅里叶正变换
    p1=fftw_plan_dft_2d( lHeight, lWidth, in1  , out1, FFTW_FORWARD, FFTW_ESTIMATE);
    p2=fftw_plan_dft_2d( lHeight, lWidth, in2  , out2, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute( p1 );                   
    fftw_execute( p2 );    
    
    //取图像2快速傅里叶变换的共轭

    for(i=0;i<lHeight*lWidth;i++)
    {
        //图像2傅里叶变换结果虚部取共轭
        out2[i][1]=-out2[i][1];
        //对应相乘
        out1[i][0] = out1[i][0] * out2[i][0] -out1[i][1] * out2[i][1];
        out1[i][1] = out1[i][0] * out2[i][1] +out1[i][1] * out2[i][0];
        //取模值
        double temp_abs;
        temp_abs=sqrt(out1[i][0]*out1[i][0]+out1[i][1]*out1[i][1])+0.00001;// 加一个小数，避免除0；
        //求互功率谱
        out1[i][0] =out1[i][0]/temp_abs;
        out1[i][1] =out1[i][1]/temp_abs;
    }

    //图像1 FFT变换和图像2 FFT共轭相乘,并计算互功率谱，结果放到out3中
    q=fftw_plan_dft_2d( lHeight, lWidth, out1, out3, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(q);

    
    //找出out3实部的极大值，以及对应的位置
    double max=-1;
    int k = 0;
    for(i=0;i<lHeight;i++)
    {
        for (j=0;j<lWidth;j++,k++)
        {
        //图像2傅里叶变换结果虚部取共轭    
            if (out3[k][0]>max)
            {
                xout=i;
                yout=j;
                max=out3[k][0];
            }
        }
	}
    
    fftw_destroy_plan(p1);
    fftw_destroy_plan(p2);
    fftw_destroy_plan(q);
    fftw_free(in1);
    fftw_free(in2);
    fftw_free(out1);
    fftw_free(out2); 	
    fftw_free(out3); 
}


/***********************************************************************/
/* 函数名称:   
/*         Fourier_pu_distrbution()              
/* 函数参数:  
/*   unsigned char *image     图像矩阵
/*   long lHeight             图像高度  
/*   long lWidth              图像宽度
/*   int N                    环的个数
/*   double *pu_percent       环内的能量向量
/* 说明：
/*    长方环傅里叶能量谱分布，为了便于统计与计算，这里N取20，
/*    并且长和宽都必须是20的倍数，并且长和宽相等 ,例如40*20；     
/************************************************************************/
void Fourier_pu_distrbution(unsigned char *image,long lHeight,long lWidth,int N,double *pu_percent)
{
    //长和宽相等，并且是2*N的倍数；
    int i,j;
    //先进行傅里叶变换
    int linebyte = (lWidth+3)/4*4;
    double *image_fft = new double[linebyte*lHeight];
    memset(image_fft,0,linebyte*lHeight*sizeof(double));
    double *image_fftshift = new double[linebyte*lHeight];
    memset(image_fftshift,0,linebyte*lHeight*sizeof(double));
   
    double * Int_image = new double[lHeight*lWidth];
    memset(Int_image,0,lHeight*lWidth*sizeof(double));

    //傅里叶正变换
    FFT(image,image_fft,lWidth,lHeight);

    //傅里叶频谱搬移,fftshift_2这个要快一些
    fftshift_2(image_fft,image_fftshift,lWidth,lHeight);
   
    //对傅里叶变换频谱搬移后的结果进行积分图计算
    Image_integral(image_fftshift,Int_image,lWidth,lHeight);
     
    //先存储N个中间变量
    double *temp = new double[N];
    int stepy = lHeight/(2*N);        //lHeight必须是2*N的整数倍，且长宽相等
    int stepx = lWidth/(2*N);             
    
    temp[0] = Int_image[lWidth*lHeight-1];   //第一个值是整幅图像的和值；  

    for (i =1;i<N;i++)
    {
        temp[i] = Int_image[lWidth*(lHeight-stepy*i-1)+ lWidth-stepx*i-1] 
                 +Int_image[lWidth*(stepy*i-1)+stepx*i-1] 
                 -Int_image[lWidth*(lHeight-stepy*i-1)+ stepx*i-1] 
                 -Int_image[lWidth*(stepy*i-1)+lWidth-stepx*i-1];
    }

    pu_percent [0] = 1;      //中间最高能量值；

    for (j =N-2;j>=0;j--)
    {
        pu_percent [N-j-1] = temp[j] - temp[j+1];
        pu_percent [N-j-1] /= temp [N-1];
    }
        

    delete [] image_fft;
    image_fft = NULL;
    delete []image_fftshift;
    image_fftshift = NULL;
    delete [] Int_image;
    Int_image = NULL;
    delete [] temp;
    temp = NULL;

}

void Image_integral(double *image,double *integral,int lwidth,int lheight)
{
    int x,y;
    int k = 0;
    //先赋值
    memcpy(integral,image,lwidth*lheight*sizeof(double));

    //先计算列数据
    for (y=1;y<lheight;y++)
    {
        for (x=0;x<lwidth;x++)
        {
            int k1 = y*lwidth+x;
            int k2 = k1-lwidth;
            integral[k1] += integral[k2];
        }
    }


    //再计算列数据
    for (y=0;y<lheight;y++)
    {
        for (x=1;x<lwidth;x++)
        {
            int k1 = y*lwidth+x;
            int k2 = k1-1;
            integral[k1] += integral[k2];
        }
    }

}





/***********************************************************************/
/* 函数名称:   
/*         Edge_Match()              
/* 函数参数:  
/*   unsigned char *src1           图像1的像素矩阵
/*   int width1                    图像1的宽度
/*   int height1                   图像1的高度
/*   unsigned char *src2           图像2的像素矩阵  
/*   int width2                    图像2的宽度  
/*   int height2                   图像2的高度 
/*   int &xout                     输出平移量x
/*   int &yout                     输出平移量y 
/* 说明：
/*    该函数利用椭圆傅里叶描述子提取目标的傅里叶描述向量     
/************************************************************************/
void Edge_Match(unsigned char *src1,int width1,int height1,unsigned char *src2,int width2,int height2,int &xout, int &yout)
{
    //就只是一个简单的卷积，找出卷积位置最大值就可以确定目标的最终位置
    //out = conv2(wedge,tedge,'same');
    //o       = max(max(out));
    // output  = (1/o)*out;
    //pixel = find(output == 1)
}






