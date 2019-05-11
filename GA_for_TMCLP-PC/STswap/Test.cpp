#include "stdafx.h"
#include "Test.h"
#include<fstream>


Test::Test()
{
	PIPSpoint = new MyPoint[1012];
	cover1 = new double*[1012];
	for (int i = 0; i < 1012; i++)
	{
		cover1[i] = new double[621];
	}
	DemandArea = new double[621];

	
	//读取PIPS点
	int PIPSnum = 1012;
	int demandnum = 621;
	ifstream lever1("pips.txt");
	double x, y;
	for (int i = 0; i <PIPSnum; i++){
		lever1 >> x >> y;
		PIPSpoint[i].x = x;
		PIPSpoint[i].y = y;

	}
	lever1.close();	
	//需求区域面积
	ifstream demand("demand.txt");
	for (int i = 0; i < demandnum; i++)
	{
		demand >> DemandArea[i];
	}
	demand.close();
	//部分覆盖
	//ifstream Co1("cover1GApmp.txt");
	ifstream Co1("Cover1ForTpmp.txt"); 
	for (int i = 0; i < PIPSnum; i++)
	{
		for (int j = 0; j < demandnum; j++)
		{
			Co1 >> cover1[i][j];
		}
	}
	//PIPS的时空覆盖率
	ifstream st("PIPSSTrate.txt");
	for (int i = 0; i < PIPSnum; i++)
	{
		st >> PIPSpoint[i].STrate;
		PIPSpoint[i].index = i;
	}
	


}


Test::~Test()
{
	//delete[]PIPSpoint;
	//delete[]DemandArea;
	//for (int i = 0; i < 1012; i++)
	//{
	//	delete[]cover1[i];
	//}
	//delete[]cover1;
}




MyPoint* Test:: GetPIPSByIndex(int no)
{
	if (no <PIPS){
		m_point = PIPSpoint[no];
		return &m_point;
	}
	else
	{
		return NULL;
	}	
}