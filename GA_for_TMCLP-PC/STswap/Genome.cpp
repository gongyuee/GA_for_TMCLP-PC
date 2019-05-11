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


/// ����������������
///  m_genes[]�����һ�����÷���
void Genome::CreateGenesIntVector(int length)
{
	m_length = length;
	int num[PIPS];
	int i, r = PIPS - 1;
	int n;
	int tmp;

	for (i = 0; i < PIPS; i++)//��ʼ��������飬0~m_length��ȡֵ��ΧΪ��ѡ�㼯��candidate sites
	{
		num[i] = i;
	}

	///ѭ�����򳤶ȵĴ���
	for (i = 0; i < m_length; i++)
	{
		n = rand() % (r - 1);    //�������һ��0~r֮�������r�ĳ�ʼֵ��m_length-1������ѡ�㼯�����һ����λ���
		m_genes[i] = num[n];        //�Ѳ��������������num���±긳��newNum����
		tmp = num[n];               //Ȼ���num[n]�������һ��������(num[r])�����Ǳ����ظ�
		num[n] = num[r];
		num[r] = tmp;
		r--;                        //�Լ����´β�����������Ϳ��Դ�0��8�ˣ�		
	}
	
}

//������Ӧ��
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

//�ٽ�����
Genome* Genome::CrossoverNearbySwap_Intvector(Genome genome2, Test TA,Genome ta[])
{
	Genome child1;
	child1.m_length = m_genomeSize;

	Genome child2;
	child2.m_length = m_genomeSize;

	//�������л���,��ȡ���л���
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
	//��ȡps1���л���
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
	//��ȡps2���л���
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

	//������������ģ�(����������������ͬʱ���Ž��н���)

	int* cs2_genes=new int[m_genomeSize];

	if (same<m_length)
	{
		////��ʼ��cs1_genes
		//for (int i = 0; i < m_genomeSize-same; i++)
		//	cs1_genes[i]=ps1_genes[i];

		//��ʼ��TA���Ե���GetPointByIndex����
		//Point[,] area = TA.AreaInitialize();

		for (int i = 0; i < different; i++)
		{
			MyPoint p1 = TA.PIPSpoint[ps1_genes[i]];
				//*TA.GetPIPSByIndex(ps1_genes[i]);

			//����p1��ps2_genes�и��������Ӧ��λ�ľ���
			double distance[m_genomeSize];
			for (int j = 0; j < different - i; j++)
			{				
				MyPoint p2 = TA.PIPSpoint[ps2_genes[j]];
				//*TA.GetPIPSByIndex(ps2_genes[j]);
				distance[j] = sqrt(pow((p1.x - p2.x), 2) + pow((p1.y - p2.y), 2));																	
			}
			//������������Сֵ
			double min_distance = distance[0];
			for (int j = 1; j <different - i; j++)
			{
				if (min_distance > distance[j])
					min_distance = distance[j];
			}
			//����С�����Ӧ�ĵ�λ��p1���paring��ͬʱ�޳��Ѿ���Ե�ps2_genes��λ
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

		//�����л��򸳻�
		for (int i = 0; i < same; i++)
		{
			ps1_genes[different + i] = same_genes[i];
			cs2_genes[different+i]=same_genes[i];
		}

		//����Ժõ�cs1_genes��cs2_genes�У�ÿһ��paring�Ľ�������Ϊ0.5
		//Random rand = new Random();
		for (int i = 0; i <m_genomeSize; i++)
		{
			if (rand() % 100 / (float)100 < 0.5)//����
			{
				child1.m_genes[i] = cs2_genes[i];
				child2.m_genes[i] = ps1_genes[i];
			
			}
			else//������
			{ 
				child1.m_genes[i] = ps1_genes[i];
				child2.m_genes[i] = cs2_genes[i];
			}
		}
	}
	else
	{
		//�����������������ȫ��ͬ����������
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

//����ʱ�ո�����
Genome* Genome::CrossoverST_Intvector(Genome genome2, Test TA, Genome ta[])
{

	Genome child1;
	child1.m_length = m_genomeSize;

	Genome child2;
	child2.m_length = m_genomeSize;

	//�������л���,��ȡ���л���
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
	//��ȡps1���л���
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
	//��ȡps2���л���
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


	//������������ģ�(����������������ͬʱ���Ž��н���)



	if (same<m_length)
	{
		//ps1��ps2����
		MyPoint* p1 = new MyPoint[different];
		MyPoint* p2 = new MyPoint[different];
		for (int i = 0; i < different; i++)
		{
			p1[i] = TA.PIPSpoint[ps1_genes[i]];
			p2[i] = TA.PIPSpoint[ps2_genes[i]];
		}
		sort(p1, p1 + different - 1, comp1);//С->��
		sort(p2, p2 + different - 1, comp2);//��->С
		for (int i = 0; i < different; i++)
		{
			//if (p1[i].STrate < p2[i].STrate )//����	
			//{
			//	ps1_genes[i] = p2[i].index;
			//	ps2_genes[i] = p1[i].index;
			//}
			//else
			//	{
			//		ps1_genes[i] = p1[i].index;
			//		ps2_genes[i] = p2[i].index;
			//	}
			if (p1[i].STrate < p2[i].STrate )//����			
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
		//�����л��򸳻�
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
		//�����������������ȫ��ͬ����������
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

	//�������л���,��ȡ���л���
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
	//��ȡps1���л���
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
	//��ȡps2���л���
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

	//������������ģ�(����������������ͬʱ���Ž��н���)

	int cs1_genes[m_genomeSize];
	int cs2_genes[m_genomeSize];
	if (same < m_genomeSize)
	{
		//������ɽ����λ
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

		//�����л��򸳻�
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
		//�����������������ȫ��ͬ����������
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

	//�������л���,��ȡ���л���
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
	//��ȡps1���л���
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
	//��ȡps2���л���
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

	//������������ģ�(����������������ͬʱ���Ž��н���)

	int cs1_genes[m_genomeSize];
	int cs2_genes[m_genomeSize];
	if (same <m_genomeSize)
	{
		//����cs1_genes
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

		//����cs2_genes
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

		//�����л��򸳻�
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
		//�����������������ȫ��ͬ����������
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

	//�������л���,��ȡ���л���
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
	//��ȡps1���л���
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
	//��ȡps2���л���
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

	//������������ģ�(����������������ͬʱ���Ž��н���)

	int cs1_genes[m_genomeSize];
	int cs2_genes[m_genomeSize];
	double fp1 = m_fitness;
	double fp2 = genome2.m_fitness;
	double pp1 = fp2 / (fp1 + fp2);
	double pp2 = 1 - pp1;

	if (same < m_genomeSize)
	{
		//����cs1_genes
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

		//����cs2_genes
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
		//�����л��򸳻�
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
		//�����������������ȫ��ͬ����������
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

//����
void Genome::Mutate(Test TA)
{
	int index = rand() % (m_length - 1);//rand1.Next(0, m_length - 1);      //���ȷ�������������
	int candidates[PIPS];
	

	for (int i = 0; i < PIPS; i++)
	{
		if (TA.PIPSpoint[i].STrate>TA.PIPSpoint[m_genes[index]].STrate)
		{
			candidates[i]=-1;
		}
		else candidates[i] = i;
	}

	//���ɲ�����ǰ�������ĵ�λ���򼯺�
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
