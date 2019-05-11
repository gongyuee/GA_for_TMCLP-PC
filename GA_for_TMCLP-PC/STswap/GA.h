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

	 //int m_populationSize;        //��Ⱥ��С
	 //double m_crossoverRate;      //�������
	 //double m_mutationRate;       //�������
	 //int m_generationLimit;       //������������

	// static int generation_Index;
	// int m_genomeSize;         //�������С����Ⱦɫ����볤�ȣ�������������	
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
	 T randT(T Lower, T Upper); //���������������������

	//double m_totalFitness;        //��¼һ����Ⱥ��ÿ���������Ӧ��֮��
	//void DistinctPopulation();
};

