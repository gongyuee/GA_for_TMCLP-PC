#pragma once
#include<string>
#include"Test.h"
#include"all.h"
class Genome
{
public:
	Genome();
	~Genome();
	int m_genes[m_genomeSize];                       //基因数组
	int m_length;                       //基因组编码长度
	//static double m_mutationRate;       //变异率
	double m_fitness;//融合交叉中用到
	double SumFitness;
	double FitnessFunction(Test TA);
	void CreateGenesIntVector(int length);
	Genome* CrossoverNearbySwap_Intvector(Genome genome2, Test GA, Genome ta[]);
	Genome* CrossoverST_Intvector(Genome genome2, Test GA, Genome ta[]);
	Genome* CrossoverUniform_Intvector(Genome genome2, Test GA, Genome ta[]);
	Genome* CrossoverOnePoint_Intvector(Genome genome2, Test GA,Genome ta[]);
	Genome* CrossoverFusion_Intvector(Genome genome2, Test GA, Genome ta[]);

	void Mutate(Test GA);
};

