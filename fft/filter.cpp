// filter.cpp: implementation of the filter class.
//
//////////////////////////////////////////////////////////////////////

#include "../main/stdafx.h"
#include "../main/IFTrack.h"
#include "filter.h"


#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////
filter::filter()
{
	coeff = NULL; size = 0; firstIndex = 0;
}

filter::filter(int size, int firstIndex, double *coeff)
{
	init(size,firstIndex,coeff);
}

filter::filter(const filter &filter)
{
	coeff=NULL; 
	copy(filter);
}

filter::~filter()
{
	if (coeff != NULL)
		delete [] coeff;
}

//////////////////////////////////////////////////////////////////////
//���ܣ��˲����ĳ�ʼ��
//���������filterSize��С���ĳ��ȣ�data��С����ϵ����filterFirst�˲���ʼ��
//���������size ��С��firstIndex ��ʼ����coeff ��ϵ������
//û�з���ֵ
//////////////////////////////////////////////////////////////////////
void filter::init(int filterSize, int filterFirst, double *data)
{
	size = filterSize;
	firstIndex = filterFirst;
	center = -firstIndex;
	
	coeff = new double [size];
	if (data != NULL) 
	{
		for (int i = 0; i < size; i++)
			coeff[i] = data[i];
	}
	else//�ж��Ƿ���С��ϵ��
	{
		for (int i = 0; i < size; i++)
			coeff[i] = 0;
	}
}

//////////////////////////////////////////////////////////////////////
//���ܣ��˲����ĳ�ʼ��
//���������filterһ��С��ʵ��
//�����������һ��С����ʵ��
//û�з���ֵ
//////////////////////////////////////////////////////////////////////
void filter::copy(const filter &filter)
{
	if (coeff != NULL)
		delete [] coeff;//��ȫ��С��ϵ��ָ����ж�
	//��ʼ��
	init (filter.size, filter.firstIndex, filter.coeff);
}



