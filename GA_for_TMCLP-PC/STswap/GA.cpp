#include "stdafx.h"
#include<algorithm>
#include "GA.h"
#include <iostream>
#include <time.h> 
#include<string>
#include<fstream>

#include<cstdlib>
using namespace std;


GA::GA()
{

}


GA::~GA()
{
}



//GA::GA (int populationSize, double crossoverRate, double mutationRate, int generationLimit)
//{
//	m_populationSize = populationSize;
//	m_crossoverRate = crossoverRate;
//	m_mutationRate = mutationRate;
//	m_generationLimit = generationLimit;
//}
//ģ�庯�����������ɸ��������ϵ���������
template<typename T>
T GA:: randT(T Lower, T Upper)
{
	return rand() / (double)RAND_MAX *(Upper - Lower) + Lower;
}
void GA::Go(Test TA)
{
	//ofstream out("result.txt");
	//m_genomeSize = genomeSize;
	srand((unsigned)time(NULL));
	CreateGenomes(TA);        //���ɳ�ʼ��Ⱥ������
	//generation_Index = 0;

	RankPopulation(TA);       //������Ӧ�ȶ���Ⱥ�������������С����

	for (int i = 0; i < m_generationLimit; i++)
	{		
		if(m_thisGeneration[m_genomeSize-1].m_fitness>353753.8)
		{
			cout<<"gengeration="<<i<<endl;
			break;
		}
		CreateNextGeneration(TA);     //������һ����Ⱥ

		

		//RankPopulation(TA);           //����Ⱥ��������
	}

	
}

/// <summary>
/// ���ɳ�ʼ������
/// </summary>
 void GA::CreateGenomes(Test TA)
{
	
	for (int i = 0; i < m_populationSize; i++)
	{
		//Genome gene = new Genome(m_genomeSize, "binary");
		Genome gene;    //������������
		gene.CreateGenesIntVector(m_genomeSize);
		gene.m_fitness = gene.FitnessFunction(TA);
		m_thisGeneration[i]=gene;
		                                   //ȷ����������Ӳ�ͬ
	}
}

 bool comp(Genome x, Genome y)
 {
	 //if (((Genome)x).m_fitness > ((Genome)y).m_fitness)
	 //{
		// return 1;       //��Ӧ��x����y������1
	 //}
	 //else if (((Genome)x).m_fitness < ((Genome)y).m_fitness)
	 //{
		// return -1;       //��Ӧ��x����y������0
	 //}
	 //else
	 //{
		// return 0;      //��Ӧ��xС��y������-1
	 //}
	 return ((Genome)x).m_fitness < ((Genome)y).m_fitness;
 }

 //��Ⱥ����
 void GA::RankPopulation(Test TA)
 {
	 sort(m_thisGeneration, m_thisGeneration + m_populationSize-1, comp);

 }

 void GA::CreateNextGeneration(Test TA)
 {

	// m_wholeGeneration = new Genome[m_populationSize];
	 for (int i = 0; i < m_populationSize; i++)
		 m_wholeGeneration[i]=m_thisGeneration[i];
	 Selection();

	 for (int i = 0; i < m_populationSize; i = i + 1)
	 {
				 

		 Genome parent1, parent2, child;
		 if(i==m_populationSize-1)
		 {
			parent1 = m_thisGeneration[i];
			parent2 = m_thisGeneration[0];
		 }
		 else
		 {
			parent1 = m_thisGeneration[i];
			parent2 = m_thisGeneration[i + 1];
		 }
		

		 //**********
		 //***����***
		 //**********
		 if ((double)(rand() % 101) / 101 < m_crossoverRate)
		 {
			
			 Genome te[2];
			// parent1.CrossoverNearbySwap_Intvector(parent2, TA,te);
			 parent1.CrossoverST_Intvector(parent2, TA, te);
			// parent1.CrossoverOnePoint_Intvector(parent2, TA, te);
			// parent1.CrossoverFusion_Intvector(parent2, TA, te);
			// parent1.CrossoverUniform_Intvector(parent2, TA, te);
			 child = te[1];
			/// child2 = te[1];
		 }
		 else
		 {
			 if(parent1.m_fitness<parent2.m_fitness)
				 child=parent2;
			 else
				 child=parent1;
			
		 }
		 
		 //**********
		 //***����***
		 //**********
		 if ((double)(rand() % 101) / 101< m_mutationRate)
		 {
			 //����Ӧ�Ƚ�С�ĺ�����б���
			 /*if (child1.m_fitness < child2.m_fitness)
			 {
				 child1.Mutate();
				 child1.m_fitness = child1.FitnessFunction(TA);
			 }
			 else
			 {
				 child2.Mutate();
				 child2.m_fitness = child2.FitnessFunction(TA);
			 }*/
			 child.Mutate(TA);
			 child.m_fitness = child.FitnessFunction(TA);
			
				 
		 }				
		 m_nextGeneration[i]=child;
		////// m_nextGeneration[i+1]=child2;
	 }
	 //��̬�滻����
	 for (int i = 0; i < m_populationSize; i++)
		 m_wholeGeneration[i + m_populationSize] = m_nextGeneration[i];
	 sort(m_wholeGeneration, m_wholeGeneration + m_populationSize*2-1, comp);
	 
	
			// DistinctPopulation();
	 for (int i = 0; i < m_populationSize; i++)
	 {
		 /*if ((double)(rand() % 101) / 101< m_mutationRate)
		 {
			 m_wholeGeneration[i + m_populationSize].Mutate();
			m_wholeGeneration[i + m_populationSize].m_fitness = m_wholeGeneration[i + m_populationSize].FitnessFunction(TA);
		 }*/
		 
		  m_thisGeneration[i]=m_wholeGeneration[i + m_populationSize];
	 }
		

 }

 //(1)�����Ⱥ����ÿ���������Ӧ��f(i=1��2������M)��MΪȺ���С��
 //(2)�����ÿ�����屻�Ŵ�����һ��Ⱥ���еĸ��ʣ�
 //(3)�����ÿ��������ۻ����ʣ�
 //(4)��[0��1]�����ڲ���һ�����ȷֲ���α�����r��
 //(5)��r<q[1]����ѡ�����1������ѡ�����k��ʹ�ã�q[k - 1]<r��q[k] ������
 //(6)�ظ�(4)��(5)��M��
 //ѡ���н���Ȩ�ĸ���
 
 void GA::Selection()
 {
	  Genome NewPopulation[m_populationSize + 1];//��ʱ�����ѡ�ĺ������
	 const double a = 0.0;
	 const double b = 1.0;
	 int i;
	 int j;
	 int mem;
	 double p;
	 double sum;

	 sum = 0.0;
	 for (mem = 0; mem < m_populationSize; mem++)
	 {
		 sum = sum + m_thisGeneration[mem].m_fitness;
	 }
	
	
	 m_thisGeneration[0].SumFitness = m_thisGeneration[0].m_fitness/sum;
	 for (mem = 1; mem < m_populationSize; mem++)
	 {
		 m_thisGeneration[mem].SumFitness = m_thisGeneration[mem - 1].SumFitness +
			 m_thisGeneration[mem].m_fitness / sum;
	 }
	 
	 for (i = 0; i < m_populationSize; i++)
	 {
		 p = rand() % 100 / (float)100;//0��1֮��������
		 if (p < m_thisGeneration[0].SumFitness)
		 {
			 NewPopulation[i] = m_thisGeneration[0];
		 }
		 else
		 {
			 for (j = 0; j < m_populationSize; j++)
			 {
				 if (m_thisGeneration[j].SumFitness <= p && p < m_thisGeneration[j + 1].SumFitness)
				 {
					 NewPopulation[i] = m_thisGeneration[j + 1];
				 }
			 }
		 }
	 }
	 //���º������ 
	 for (i = 0; i < m_populationSize; i++)
	 {
		 m_thisGeneration[i] = NewPopulation[i];
	 }
	 return;
	
		
 }

 
 