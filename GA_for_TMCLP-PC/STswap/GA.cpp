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
//模板函数，用于生成各种区间上的数据类型
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
	CreateGenomes(TA);        //生成初始种群基因组
	//generation_Index = 0;

	RankPopulation(TA);       //根据适应度对种群个体进行排序，由小到大

	for (int i = 0; i < m_generationLimit; i++)
	{		
		if(m_thisGeneration[m_genomeSize-1].m_fitness>353753.8)
		{
			cout<<"gengeration="<<i<<endl;
			break;
		}
		CreateNextGeneration(TA);     //生成下一代种群

		

		//RankPopulation(TA);           //对种群进行排序
	}

	
}

/// <summary>
/// 生成初始基因组
/// </summary>
 void GA::CreateGenomes(Test TA)
{
	
	for (int i = 0; i < m_populationSize; i++)
	{
		//Genome gene = new Genome(m_genomeSize, "binary");
		Genome gene;    //整数向量编码
		gene.CreateGenesIntVector(m_genomeSize);
		gene.m_fitness = gene.FitnessFunction(TA);
		m_thisGeneration[i]=gene;
		                                   //确保随机数种子不同
	}
}

 bool comp(Genome x, Genome y)
 {
	 //if (((Genome)x).m_fitness > ((Genome)y).m_fitness)
	 //{
		// return 1;       //适应度x大于y，返回1
	 //}
	 //else if (((Genome)x).m_fitness < ((Genome)y).m_fitness)
	 //{
		// return -1;       //适应度x等于y，返回0
	 //}
	 //else
	 //{
		// return 0;      //适应度x小于y，返回-1
	 //}
	 return ((Genome)x).m_fitness < ((Genome)y).m_fitness;
 }

 //种群排序
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
		 //***交叉***
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
		 //***变异***
		 //**********
		 if ((double)(rand() % 101) / 101< m_mutationRate)
		 {
			 //对适应度较小的后代进行变异
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
	 //稳态替换方法
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

 //(1)计算出群体中每个个体的适应度f(i=1，2，…，M)，M为群体大小；
 //(2)计算出每个个体被遗传到下一代群体中的概率；
 //(3)计算出每个个体的累积概率；
 //(4)在[0，1]区间内产生一个均匀分布的伪随机数r；
 //(5)若r<q[1]，则选择个体1，否则，选择个体k，使得：q[k - 1]<r≤q[k] 成立；
 //(6)重复(4)、(5)共M次
 //选择有交配权的父代
 
 void GA::Selection()
 {
	  Genome NewPopulation[m_populationSize + 1];//临时存放挑选的后代个体
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
		 p = rand() % 100 / (float)100;//0和1之间的随机数
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
	 //更新后代个体 
	 for (i = 0; i < m_populationSize; i++)
	 {
		 m_thisGeneration[i] = NewPopulation[i];
	 }
	 return;
	
		
 }

 
 