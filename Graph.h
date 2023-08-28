#include <list>
#include <iostream>
using namespace std;
#include "DiffusionData.h"
//Graph Class for undirected Graph
class Graph
{
    public:
    	int n;//Number of vertices of the graph
        list <int> * adjList;//Adjacency List
        list <int> * radjList;//Adjacency List
        float ** weightMat;//Adjacency Matrix
	float *it, *dt;
	//Both Adjacency List and Adjacency matrix are maintained
	Graph(int);
	~Graph();
	void addEdge(int, int, float);//Function to add an edge to the graph
	int * getInDegrees();//Function to extract degree of all the vertices (as an array)
	int * getOutDegrees();
	DiffusionData * runDiffussion(int rsc, int *RS, int tsc, int *TS); //diffusion function
	void updateWeights();
	void setFixedTH(float itv, float dtv);
	void setRandomTH();
	list<int> * getInfluentialNodes(int rsc, int *RS, int tsc, int *TS);
	int * pruning(int rsc, int *RS, int tsc, int *TS);
	int * minGreedy(int rsc, int *RS, int k);
	int * minGreedyProxy(int rsc, int *RS, int k);
	DiffusionData * runDiffussionCtr(int rsc, int *RS, int *RI);
	int * contrId(int rsc, int *RS, int k);
	int * contrIdProxy(int rsc, int *RS, int k);
	int * topkCentrality(char * file, int k, int rsc, int *RS);
	int * GA(char *cFileName, int src, int *RS, int k, int initPSize, int totPSize, int noIters);
	void generatePopulation(int k, int initPSize, int rsc, int *RS, popwrc *population[], float cumprob[], int rsFlag[]);
	void increasePopolationwithCM(int k, int initPSize, int totPSize, int rsc, int *RS, popwrc *population[], float cumprob[], int rsFlag[]);
        void Mutate(int k, int *pop, float cumprob[], int *flag);
	int * getMaxDegNodes(int k);
	int * getRandNodes(int k);
	int * getRandTS(int k);
	float * getEcentrality(char *cFileName, int rsc, int *RS);
	int* topkECentrality(char*, int, int*, int);
	int * randPrune(char * file, int rsc, int * RS, int tsc, int Iter);
	int generateRandNode(float cumprob[]);
	DiffusionData* runDiffussionSigned(int, int*, int, int*);
	int * topkCentralityProxy(char *file, int k, int rsc, int *RS);
};
