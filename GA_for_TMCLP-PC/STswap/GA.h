#pragma once
#include<vector>
#include<stdlib.h>
#include"Genome.h"
#include"Test.h"
#include"all.h"
using namespace std;
class GA
{
public:
	GA();
	~GA();

	 //int m_populationSize;        //种群大小
	 //double m_crossoverRate;      //交叉概率
	 //double m_mutationRate;       //变异概率
	 //int m_generationLimit;       //进化代数限制

	// static int generation_Index;
	// int m_genomeSize;         //基因组大小，即染色体编码长度，整数向量编码	
	// GA(int populationSize, double crossoverRate, double mutationRate, int generationLimit);
	 void Go(Test TA);
	 void CreateGenomes(Test TA);
	 void RankPopulation(Test TA);
	 void CreateNextGeneration(Test GA);
	// int* BinaryTourmentSelection();
	 void Selection();


	 Genome m_thisGeneration[m_populationSize];
	 Genome m_nextGeneration[m_populationSize];
	 Genome m_wholeGeneration[m_populationSize*2];
	 double m_fitnessTable[m_populationSize];
	 template<typename T>
	 T randT(T Lower, T Upper); //产生任意类型随机数函数

	//double m_totalFitness;        //记录一个种群中每个个体的适应度之和
	//void DistinctPopulation();
};

