#pragma once
#include "MyPoint.h"

#include"all.h"
using namespace std;


class Test
{
public:
	Test();//��ȡ����
	~Test();
	
	//vector <MyPoint> PIPS;//�洢��ѡ��
	//vector <double> demandcover;
	MyPoint *PIPSpoint;
	double **cover1;
	double *DemandArea;

	MyPoint* GetPIPSByIndex(int no);
	MyPoint m_point;

};

