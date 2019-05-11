// STswap.cpp : 定义控制台应用程序的入口点。
//

#include "stdafx.h"
#include "Test.h"
#include "GA.h"
#include "Genome.h"
#include "MyPoint.h"
#include "Test.h"
#include<vector>
#include<iostream>
#include<fstream>
#include <time.h>
using namespace std;

string getTime()
{
	time_t timep;
	time(&timep);
	char tmp[64];
	strftime(tmp, sizeof(tmp), "%Y-%m-%d %H:%M:%S", localtime(&timep));
	return tmp;
}

int _tmain(int argc, _TCHAR* argv[])
{
	for(int i=0;i<5;i++)
	{
		Test ta;
		GA ga;
		srand((int)time(0));
		cout << getTime()<<" ";
		ga.Go(ta);
		cout << getTime() << " ";
		cout << ga.m_thisGeneration[m_populationSize - 1].m_fitness<<endl;
	}
	
	

	system("pause");
	
	return 0;
}


