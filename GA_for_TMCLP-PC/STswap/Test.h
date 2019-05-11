#pragma once
#include "MyPoint.h"

#include"all.h"
using namespace std;


class Test
{
public:
	Test();//读取数据
	~Test();
	
	//vector <MyPoint> PIPS;//存储候选点
	//vector <double> demandcover;
	MyPoint *PIPSpoint;
	double **cover1;
	double *DemandArea;

	MyPoint* GetPIPSByIndex(int no);
	MyPoint m_point;

};

