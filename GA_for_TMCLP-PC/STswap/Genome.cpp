#include "stdafx.h"
#include "Genome.h"
#include "math.h"
#include<algorithm>

#include <vector>
#include<iostream>
using namespace std;

Genome::Genome()
{
	
}


Genome::~Genome()
{
}


/// 基因整数向量编码
///  m_genes[]随机的一个配置方案
void Genome::CreateGenesIntVector(int length)
{
	m_length = length;
	int num[PIPS];
	int i, r = PIPS - 1;
	int n;
	int tmp;

	for (i = 0; i < PIPS; i++)//初始化这个数组，0~m_length，取值范围为候选点集合candidate sites
	{
		num[i] = i;
	}

	///循环基因长度的次数
	for (i = 0; i < m_length; i++)
	{
		n = rand() % (r - 1);    //随机产生一个0~r之间的数，r的初始值是m_length-1，即候选点集合最后一个点位序号
		m_genes[i] = num[n];        //把产生的随机数当成num的下标赋给newNum数组
		tmp = num[n];               //然后把num[n]和它最后一个数交换(num[r])，这是避免重复
		num[n] = num[r];
		num[r] = tmp;
		r--;                        //自减，下次产生的随机数就可以从0到8了，		
	}
	
}

//计算适应度
double Genome::FitnessFunction(Test TA)
{
	double current_fitness;
	double best_fitness;
	double ind_fitness = 0;
	
	bool jump=true;
	for (int i = 0; i <= demandsize; i++)
	{
		best_fitness = 0;
		for (int j = 0; j < m_length; j++)
		{
			
			current_fitness =TA.cover1[m_genes[j]][i] ;
			if (current_fitness>best_fitness)
				best_fitness = current_fitness;
			if (best_fitness >= (int)(TA.DemandArea[i ]*6))
			{
				jump = false;
				break; 
			}
		}
		
		ind_fitness += best_fitness;
	}
	return ind_fitness;
	
}

//临近交叉
Genome* Genome::CrossoverNearbySwap_Intvector(Genome genome2, Test TA,Genome ta[])
{
	Genome child1;
	child1.m_length = m_genomeSize;

	Genome child2;
	child2.m_length = m_genomeSize;

	//保留共有基因,提取特有基因
	int *same_genes=new int [m_genomeSize];
	int *ps1_genes=new int[m_genomeSize];
	int *ps2_genes=new int[m_genomeSize];
	int same = 0;
	for (int i = 0; i < m_length; i++)
	{
		for (int j = 0; j < genome2.m_length; j++)
		{
			if (m_genes[i] == genome2.m_genes[j])
			{
				same_genes[same]=m_genes[i];
				same++;
				break;
			}
		}
	}
	int different = m_genomeSize - same;
	//提取ps1特有基因
	int m = 0;
	int sameCount;
	for (int i = 0; i < m_length; i++)
	{		
		if (same!= 0)
		{
			sameCount = 0;
			for (int j = 0; j < same; j++)
			{
				if (m_genes[i] == same_genes[j])
				{
					sameCount++;
				}
			}
			if (sameCount == 0)
				ps1_genes[m++]=(m_genes[i]);
		}
		else
			ps1_genes[m++]=m_genes[i];
	}
	//提取ps2特有基因
	m = 0;
	for (int i = 0; i < genome2.m_length; i++)
	{
		if (same!= 0)
		{
			sameCount = 0;
			for (int j = 0; j < same; j++)
			{
				if (genome2.m_genes[i] == same_genes[j])
				{
					sameCount++;
				}
			}
			if (sameCount == 0)
				ps2_genes[m++]=genome2.m_genes[i];
		}
		else
			ps2_genes[m++]=genome2.m_genes[i];
	}

	//交叉操作，核心！(仅当两个父代不相同时，才进行交叉)

	int* cs2_genes=new int[m_genomeSize];

	if (same<m_length)
	{
		////初始化cs1_genes
		//for (int i = 0; i < m_genomeSize-same; i++)
		//	cs1_genes[i]=ps1_genes[i];

		//初始化TA，以调用GetPointByIndex函数
		//Point[,] area = TA.AreaInitialize();

		for (int i = 0; i < different; i++)
		{
			MyPoint p1 = TA.PIPSpoint[ps1_genes[i]];
				//*TA.GetPIPSByIndex(ps1_genes[i]);

			//计算p1与ps2_genes中各个基因对应点位的距离
			double distance[m_genomeSize];
			for (int j = 0; j < different - i; j++)
			{				
				MyPoint p2 = TA.PIPSpoint[ps2_genes[j]];
				//*TA.GetPIPSByIndex(ps2_genes[j]);
				distance[j] = sqrt(pow((p1.x - p2.x), 2) + pow((p1.y - p2.y), 2));																	
			}
			//求距离数组的最小值
			double min_distance = distance[0];
			for (int j = 1; j <different - i; j++)
			{
				if (min_distance > distance[j])
					min_distance = distance[j];
			}
			//把最小距离对应的点位与p1配对paring！同时剔除已经配对的ps2_genes点位
			for (int j = 0; j < different-i; j++)
			{
				if (min_distance == distance[j])
				{
					cs2_genes[i]=ps2_genes[j];
					ps2_genes[j] = ps2_genes[different - i - 1];
					break;
				}
			}
		}

		//将共有基因赋还
		for (int i = 0; i < same; i++)
		{
			ps1_genes[different + i] = same_genes[i];
			cs2_genes[different+i]=same_genes[i];
		}

		//已配对好的cs1_genes和cs2_genes中，每一对paring的交换概率为0.5
		//Random rand = new Random();
		for (int i = 0; i <m_genomeSize; i++)
		{
			if (rand() % 100 / (float)100 < 0.5)//交换
			{
				child1.m_genes[i] = cs2_genes[i];
				child2.m_genes[i] = ps1_genes[i];
			
			}
			else//不交换
			{ 
				child1.m_genes[i] = ps1_genes[i];
				child2.m_genes[i] = cs2_genes[i];
			}
		}
	}
	else
	{
		//如果两个父代基因完全相同，不做交叉
		for (int i = 0; i < m_length; i++)
		{
			child1.m_genes[i] = m_genes[i];
			child2.m_genes[i] = genome2.m_genes[i];
		}
	}
	
	ta[0] = child1;
	ta[1] = child2;
	delete[]ps1_genes;
	delete[]ps2_genes;
	delete[]same_genes;
	delete[]cs2_genes;

	return ta;
}

bool comp1(MyPoint x, MyPoint y)
{
	return x.STrate < y.STrate;
}

bool comp2(MyPoint x, MyPoint y)
{
	return x.STrate > y.STrate;
}

//基于时空覆盖率
Genome* Genome::CrossoverST_Intvector(Genome genome2, Test TA, Genome ta[])
{

	Genome child1;
	child1.m_length = m_genomeSize;

	Genome child2;
	child2.m_length = m_genomeSize;

	//保留共有基因,提取特有基因
	int *same_genes = new int[m_genomeSize];
	int *ps1_genes = new int[m_genomeSize];
	int *ps2_genes = new int[m_genomeSize];
	int same = 0;
	for (int i = 0; i < m_length; i++)
	{
		for (int j = 0; j < genome2.m_length; j++)
		{
			if (m_genes[i] == genome2.m_genes[j])
			{
				same_genes[same] = m_genes[i];
				same++;
				break;
			}
		}
	}
	int different = m_genomeSize - same;
	//提取ps1特有基因
	int m = 0;
	int sameCount;
	for (int i = 0; i < m_length; i++)
	{
		if (same != 0)
		{
			sameCount = 0;
			for (int j = 0; j < same; j++)
			{
				if (m_genes[i] == same_genes[j])
				{
					sameCount++;
				}
			}
			if (sameCount == 0)
				ps1_genes[m++] = (m_genes[i]);
		}
		else
			ps1_genes[m++] = m_genes[i];
	}
	//提取ps2特有基因
	m = 0;
	for (int i = 0; i < genome2.m_length; i++)
	{
		if (same != 0)
		{
			sameCount = 0;
			for (int j = 0; j < same; j++)
			{
				if (genome2.m_genes[i] == same_genes[j])
				{
					sameCount++;
				}
			}
			if (sameCount == 0)
				ps2_genes[m++] = genome2.m_genes[i];
		}
		else
			ps2_genes[m++] = genome2.m_genes[i];
	}


	//交叉操作，核心！(仅当两个父代不相同时，才进行交叉)



	if (same<m_length)
	{
		//ps1、ps2排序
		MyPoint* p1 = new MyPoint[different];
		MyPoint* p2 = new MyPoint[different];
		for (int i = 0; i < different; i++)
		{
			p1[i] = TA.PIPSpoint[ps1_genes[i]];
			p2[i] = TA.PIPSpoint[ps2_genes[i]];
		}
		sort(p1, p1 + different - 1, comp1);//小->大
		sort(p2, p2 + different - 1, comp2);//大->小
		for (int i = 0; i < different; i++)
		{
			//if (p1[i].STrate < p2[i].STrate )//交换	
			//{
			//	ps1_genes[i] = p2[i].index;
			//	ps2_genes[i] = p1[i].index;
			//}
			//else
			//	{
			//		ps1_genes[i] = p1[i].index;
			//		ps2_genes[i] = p2[i].index;
			//	}
			if (p1[i].STrate < p2[i].STrate )//交换			
			{
				if ((rand() % 100 / (float)100) < 0.8)
				{
					ps1_genes[i] = p2[i].index;
					ps2_genes[i] = p1[i].index;
				}
				else
				{
					ps1_genes[i] = p1[i].index;
					ps2_genes[i] = p2[i].index;
				}
				
			}
			else
				if ((rand() % 100 / (float)100) < 0.3)
				{
					ps1_genes[i] = p2[i].index;
					ps2_genes[i] = p1[i].index;
				}
				else
				{
					ps1_genes[i] = p1[i].index;
					ps2_genes[i] = p2[i].index;
				}
		}
		delete[]p1; 
		delete[]p2;
		//将共有基因赋还
		for (int i = 0; i < same; i++)
		{
			ps1_genes[different + i] = same_genes[i];
			ps2_genes[different + i] = same_genes[i];
		}
		for (int i = 0; i < m_length; i++)
		{
			
			child1.m_genes[i] = ps1_genes[i];
			child2.m_genes[i] = ps2_genes[i];
		}
		
	}
	else
	{
		//如果两个父代基因完全相同，不做交叉
		for (int i = 0; i < m_length; i++)
		{
			child1.m_genes[i] = m_genes[i];
			child2.m_genes[i] = genome2.m_genes[i];
		}
	}

	child1.m_fitness = child1.FitnessFunction(TA);
	child2.m_fitness = child2.FitnessFunction(TA);
	if(child1.m_fitness <child2.m_fitness )
	{
		ta[0] = child1;
		ta[1] = child2;}
	else
	{
		ta[0] = child2;
		ta[1] = child1;
	}


	delete[]ps1_genes;
	delete[]ps2_genes;
	delete[]same_genes;
	

	return ta;
}


Genome* Genome::CrossoverOnePoint_Intvector(Genome genome2, Test GA,  Genome ta[])
{
	Genome child1;
	child1.m_length = m_genomeSize;

	Genome child2;
	child2.m_length = m_genomeSize;

	//保留共有基因,提取特有基因
	int *same_genes = new int[m_genomeSize];
	int *ps1_genes = new int[m_genomeSize];
	int *ps2_genes = new int[m_genomeSize];
	int same = 0;
	for (int i = 0; i < m_length; i++)
	{
		for (int j = 0; j < genome2.m_length; j++)
		{
			if (m_genes[i] == genome2.m_genes[j])
			{
				same_genes[same] = m_genes[i];
				same++;
				break;
			}
		}
	}
	int different = m_genomeSize - same;
	//提取ps1特有基因
	int m = 0;
	int sameCount;
	for (int i = 0; i < m_length; i++)
	{
		if (same != 0)
		{
			sameCount = 0;
			for (int j = 0; j < same; j++)
			{
				if (m_genes[i] == same_genes[j])
				{
					sameCount++;
				}
			}
			if (sameCount == 0)
				ps1_genes[m++] = (m_genes[i]);
		}
		else
			ps1_genes[m++] = m_genes[i];
	}
	//提取ps2特有基因
	m = 0;
	for (int i = 0; i < genome2.m_length; i++)
	{
		if (same != 0)
		{
			sameCount = 0;
			for (int j = 0; j < same; j++)
			{
				if (genome2.m_genes[i] == same_genes[j])
				{
					sameCount++;
				}
			}
			if (sameCount == 0)
				ps2_genes[m++] = genome2.m_genes[i];
		}
		else
			ps2_genes[m++] = genome2.m_genes[i];
	}

	//交叉操作，核心！(仅当两个父代不相同时，才进行交叉)

	int cs1_genes[m_genomeSize];
	int cs2_genes[m_genomeSize];
	if (same < m_genomeSize)
	{
		//随机生成交叉点位
		int pos = rand() % different - 1;
		for (int i = 0; i < m_genomeSize; i++)
		{
			if (i < pos)
			{
				cs1_genes[i]=ps1_genes[i];
				cs2_genes[i]=ps2_genes[i];
			}
			else
			{
				cs1_genes[i] = ps2_genes[i];
				cs2_genes[i] = ps1_genes[i];
			}
		}

		//将共有基因赋还
		for (int i = 0; i < same; i++)
		{
			cs1_genes[different + i] = same_genes[i];
			cs2_genes[different + i] = same_genes[i];
		}

		for (int i = 0; i < m_genomeSize; i++)
		{
			child1.m_genes[i] = cs1_genes[i];
			child2.m_genes[i] = cs2_genes[i];
		}
	}
	else
	{
		//如果两个父代基因完全相同，不做交叉
		for (int i = 0; i < m_genomeSize; i++)
		{
			child1.m_genes[i] = m_genes[i];
			child2.m_genes[i] = genome2.m_genes[i];
		}
	}
	ta[0] = child1;
	ta[1] = child2;

	return ta;
}

Genome* Genome::CrossoverUniform_Intvector(Genome genome2, Test GA, Genome ta[])
{
	Genome child1;
	child1.m_length = m_genomeSize;

	Genome child2;
	child2.m_length = m_genomeSize;

	//保留共有基因,提取特有基因
	int *same_genes = new int[m_genomeSize];
	int *ps1_genes = new int[m_genomeSize];
	int *ps2_genes = new int[m_genomeSize];
	int same = 0;
	for (int i = 0; i < m_length; i++)
	{
		for (int j = 0; j < genome2.m_length; j++)
		{
			if (m_genes[i] == genome2.m_genes[j])
			{
				same_genes[same] = m_genes[i];
				same++;
				break;
			}
		}
	}
	int different = m_genomeSize - same;
	//提取ps1特有基因
	int m = 0;
	int sameCount;
	for (int i = 0; i < m_length; i++)
	{
		if (same != 0)
		{
			sameCount = 0;
			for (int j = 0; j < same; j++)
			{
				if (m_genes[i] == same_genes[j])
				{
					sameCount++;
				}
			}
			if (sameCount == 0)
				ps1_genes[m++] = (m_genes[i]);
		}
		else
			ps1_genes[m++] = m_genes[i];
	}
	//提取ps2特有基因
	m = 0;
	for (int i = 0; i < genome2.m_length; i++)
	{
		if (same != 0)
		{
			sameCount = 0;
			for (int j = 0; j < same; j++)
			{
				if (genome2.m_genes[i] == same_genes[j])
				{
					sameCount++;
				}
			}
			if (sameCount == 0)
				ps2_genes[m++] = genome2.m_genes[i];
		}
		else
			ps2_genes[m++] = genome2.m_genes[i];
	}

	//交叉操作，核心！(仅当两个父代不相同时，才进行交叉)

	int cs1_genes[m_genomeSize];
	int cs2_genes[m_genomeSize];
	if (same <m_genomeSize)
	{
		//生成cs1_genes
		for (int i = 0; i < different; i++)
		{
			if (rand()%3 == 1)
			{
				cs1_genes[i] = ps1_genes[i];
			}
			else
			{
				cs1_genes[i] = ps2_genes[i];
			}
		}

		//生成cs2_genes
		for (int i = 0; i < different; i++)
		{
			if (rand() % 3 == 1)
			{
				cs2_genes[i]=ps1_genes[i];
			}
			else
			{
				cs2_genes[i]=ps2_genes[i];
			}
		}

		//将共有基因赋还
		for (int i = 0; i < same; i++)
		{
			cs1_genes[different + i] = same_genes[i];
			cs2_genes[different + i] = same_genes[i];
		}

		for (int i = 0; i < m_genomeSize; i++)
		{
			child1.m_genes[i] = cs1_genes[i];
			child2.m_genes[i] = cs2_genes[i];
		}
	}
	else
	{
		//如果两个父代基因完全相同，不做交叉
		for (int i = 0; i < m_genomeSize; i++)
		{
			child1.m_genes[i] = m_genes[i];
			child2.m_genes[i] = genome2.m_genes[i];
		}
	}
	ta[0] = child1;
	ta[1] = child2;

	return ta;
}

Genome* Genome::CrossoverFusion_Intvector(Genome genome2, Test GA, Genome ta[])
{
	Genome child1;
	child1.m_length = m_genomeSize;

	Genome child2;
	child2.m_length = m_genomeSize;

	//保留共有基因,提取特有基因
	int *same_genes = new int[m_genomeSize];
	int *ps1_genes = new int[m_genomeSize];
	int *ps2_genes = new int[m_genomeSize];
	int same = 0;
	for (int i = 0; i < m_length; i++)
	{
		for (int j = 0; j < genome2.m_length; j++)
		{
			if (m_genes[i] == genome2.m_genes[j])
			{
				same_genes[same] = m_genes[i];
				same++;
				break;
			}
		}
	}
	int different = m_genomeSize - same;
	//提取ps1特有基因
	int m = 0;
	int sameCount;
	for (int i = 0; i < m_length; i++)
	{
		if (same != 0)
		{
			sameCount = 0;
			for (int j = 0; j < same; j++)
			{
				if (m_genes[i] == same_genes[j])
				{
					sameCount++;
				}
			}
			if (sameCount == 0)
				ps1_genes[m++] = (m_genes[i]);
		}
		else
			ps1_genes[m++] = m_genes[i];
	}
	//提取ps2特有基因
	m = 0;
	for (int i = 0; i < genome2.m_length; i++)
	{
		if (same != 0)
		{
			sameCount = 0;
			for (int j = 0; j < same; j++)
			{
				if (genome2.m_genes[i] == same_genes[j])
				{
					sameCount++;
				}
			}
			if (sameCount == 0)
				ps2_genes[m++] = genome2.m_genes[i];
		}
		else
			ps2_genes[m++] = genome2.m_genes[i];
	}

	//交叉操作，核心！(仅当两个父代不相同时，才进行交叉)

	int cs1_genes[m_genomeSize];
	int cs2_genes[m_genomeSize];
	double fp1 = m_fitness;
	double fp2 = genome2.m_fitness;
	double pp1 = fp2 / (fp1 + fp2);
	double pp2 = 1 - pp1;

	if (same < m_genomeSize)
	{
		//生成cs1_genes
		for (int i = 0; i <different; i++)
		{
			if (rand() < pp1)
			{
				cs1_genes[i]=ps1_genes[i];
			}
			else
			{
				cs1_genes[i] = ps2_genes[i];
			}
		}

		//生成cs2_genes
		for (int i = 0; i < different; i++)
		{
			if (rand() < pp1)
			{
				cs2_genes[i] = ps1_genes[i];
			}
			else
			{
				cs2_genes[i] = ps2_genes[i];
			}
		}
		//将共有基因赋还
		for (int i = 0; i < same; i++)
		{
			cs1_genes[different + i] = same_genes[i];
			cs2_genes[different + i] = same_genes[i];
		}

		for (int i = 0; i < m_genomeSize; i++)
		{
			child1.m_genes[i] = cs1_genes[i];
			child2.m_genes[i] = cs2_genes[i];
		}
	}
	else
	{
		//如果两个父代基因完全相同，不做交叉
		for (int i = 0; i < m_genomeSize; i++)
		{
			child1.m_genes[i] = m_genes[i];
			child2.m_genes[i] = genome2.m_genes[i];
		}
	}
	ta[0] = child1;
	ta[1] = child2;

	return ta;
}

//变异
void Genome::Mutate(Test TA)
{
	int index = rand() % (m_length - 1);//rand1.Next(0, m_length - 1);      //随机确定变异基因的序号
	int candidates[PIPS];
	

	for (int i = 0; i < PIPS; i++)
	{
		if (TA.PIPSpoint[i].STrate>TA.PIPSpoint[m_genes[index]].STrate)
		{
			candidates[i]=-1;
		}
		else candidates[i] = i;
	}

	//生成不含当前个体基因的等位基因集合
	for (int i = 0; i <m_length; i++)
	{
		
		candidates[m_genes[i]] = -1;
		//candidates.erase(candidates.begin() + m_genes[i]-1);
	}
	int genes[PIPS - m_genomeSize];
	int k = 0;
	for (int i = 0; i < PIPS; i++)
	{
		if (candidates[i] != -1)
			genes[k++]=candidates[i];
	}
	int site = rand() % (PIPS-m_genomeSize - 1);// rand2.Next(0, m_length - 1);

	m_genes[index] = genes[site];
}
