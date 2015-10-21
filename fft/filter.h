// filter.h: interface for the filter class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_FILTER_H__71AB2400_9A88_4466_86F7_4947A25626C7__INCLUDED_)
#define AFX_FILTER_H__71AB2400_9A88_4466_86F7_4947A25626C7__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

class filter  
{
public:
	int size, firstIndex, center;
	double *coeff;
	
	filter();
	filter (int size, int firstIndex, double*coeff);
	filter (const filter &filter);
	virtual ~filter();
	
	void init (int filterSize, int filterFirst, double *data);
	
protected:
	void copy (const filter& filter);
};

#endif // !defined(AFX_FILTER_H__71AB2400_9A88_4466_86F7_4947A25626C7__INCLUDED_)
