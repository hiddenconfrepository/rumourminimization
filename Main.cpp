#include<iostream>
#include<list>
#include <fstream>
#include <stdlib.h>
#include <time.h>
#include<string.h>
#include "Graph.h"

using namespace std;

static void printList(list<int> *lst)
{
    cout<<""<<lst->size()<<"\n";
    list<int>::iterator iter;
    for(iter = lst->begin(); iter != lst->end(); iter++)
    {
	int v = *iter;
	cout<<v<<", ";
    }

}

static int adjustTS(int tsc, int *TS)
{
    for(int i = 0; i < tsc; i++)
    {
	if(TS[i] == -1) 
	{
	    return i;
	}
    }
    return tsc;
}
int main(int argc, char * argv[])
{
    int n, m;   
    ifstream inFile;
    inFile.open(argv[1]);
    if(!inFile)
    {
        cout << "Not able to open input file: " << argv[1];
        return 0;
    }
    srand(time(NULL));
    if(!inFile.eof())
    {
        inFile >> n >> m;
    }
    Graph *g = new Graph(n);

    while(!inFile.eof())
    {
        int u, v;
        inFile >> u >> v;
        if(u == v)
            continue;
        if(!inFile.eof())
	{
            g->addEdge(u, v, 0.0);
	}
    }
    g->updateWeights();
    int rsc = atoi(argv[2]);
    int tsc = atoi(argv[3]);
    int tsc1 = tsc;
    float it = atof(argv[4]);
    float dt = atof(argv[5]);
    //cout<<argv[1]<<","<<rsc<<","<<tsc<<",";
    if(!strcmp(argv[4], "randitdt"))
    {
	g->setRandomTH();
	//cout<<"randit,"<<"randdt,";
    }
    else
    {
	g->setFixedTH(it, dt);
	//cout<<it<<","<<dt<<",";
    }
    int *RS;
    if(!strcmp(argv[6], "randseed"))
    {
	RS = g->getRandNodes(rsc);
	//cout<<"randseed,";
    }
    else
    {
	RS = g->getMaxDegNodes(rsc);
	//cout<<"maxdegree,";
    }
    bool runmg = false, runpmg = false, runctrid = false, runpctrid = false, runcen = false, runpcen = false, runecen = false, runga = false, runrandprune = false, runall = false, centrlty = false;

    
    
    
    if(!strcmp(argv[7], "mingreedy"))
    {
	runmg = true;
    }
    else if(!strcmp(argv[7], "proxymingreedy"))
    {
	runpmg = true;
    }
    else if(!strcmp(argv[7], "contrid"))
    {
	runctrid = true;
    }
    else if(!strcmp(argv[7], "proxycontrid"))
    {
	runpctrid = true;
    }
    else if(!strcmp(argv[7], "centrality"))
    {
	centrlty = true;
	runcen = true;
    }
    else if(!strcmp(argv[7], "proxycentrality"))
    {
	centrlty = true;
	runpcen = true;
    }
    else if(!strcmp(argv[7], "ecentrality"))
    {
	centrlty = true;
	runecen = true;
    }
    else if(!strcmp(argv[7], "ga"))
    {
	centrlty = true;
	runga = true;
    }
    else if(!strcmp(argv[7], "randpruning"))
    {
	centrlty = true;
	runrandprune = true;
    }
    else if(!strcmp(argv[7], "runall"))
    {
	centrlty = true;
	runall = true;
    }
    clock_t time_req;
    int *TS = NULL;

    if(runmg || runall)
    {
	time_req = clock();
	tsc = tsc1;
	TS = g->minGreedy(rsc, RS, tsc);
	time_req = clock() - time_req;
	tsc = adjustTS(tsc, TS);
	DiffusionData *dd = g->runDiffussion(rsc, RS, tsc, TS);
	cout<<argv[1]<<","<<argv[2]<<","<<argv[3]<<","<<argv[4]<<","<<argv[5]<<","<<argv[6]<<","<<"mingreedy"<<","<<dd->rac+rsc<<","<<dd->tac+tsc<<","<<(float)time_req/CLOCKS_PER_SEC<<",";
	delete dd;

	time_req = clock();
	TS = g->pruning(rsc, RS, tsc, TS);
	time_req = clock() - time_req;
	dd = g->runDiffussion(rsc, RS, tsc, TS);
        cout<<dd->rac+rsc<<","<<dd->tac+tsc<<","<<(float)time_req/CLOCKS_PER_SEC<<endl;
	delete dd;
	delete [] TS;
    }
    if(runpmg || runall)
    {
	time_req = clock();
	tsc = tsc1;
	TS = g->minGreedyProxy(rsc, RS, tsc);
	time_req = clock() - time_req;
	tsc = adjustTS(tsc, TS);
	DiffusionData *dd = g->runDiffussion(rsc, RS, tsc, TS);
	cout<<argv[1]<<","<<argv[2]<<","<<argv[3]<<","<<argv[4]<<","<<argv[5]<<","<<argv[6]<<","<<"proxymingreedy"<<","<<dd->rac+rsc<<","<<dd->tac+tsc<<","<<(float)time_req/CLOCKS_PER_SEC<<",";
	delete dd;

	time_req = clock();
	TS = g->pruning(rsc, RS, tsc, TS);
	time_req = clock() - time_req;
	dd = g->runDiffussion(rsc, RS, tsc, TS);
        cout<<dd->rac+rsc<<","<<dd->tac+tsc<<","<<(float)time_req/CLOCKS_PER_SEC<<endl;
	delete dd;
	delete [] TS;
    }

    if(runctrid || runall)
    {
	time_req = clock();
	tsc = tsc1;
	TS = g->contrId(rsc, RS, tsc);
	time_req = clock() - time_req;
	tsc = adjustTS(tsc, TS);
	DiffusionData *dd = g->runDiffussion(rsc, RS, tsc, TS);
	cout<<argv[1]<<","<<argv[2]<<","<<argv[3]<<","<<argv[4]<<","<<argv[5]<<","<<argv[6]<<","<<"contrid"<<","<<dd->rac+rsc<<","<<dd->tac+tsc<<","<<(float)time_req/CLOCKS_PER_SEC<<",";
	delete dd;

	time_req = clock();
	TS = g->pruning(rsc, RS, tsc, TS);
	time_req = clock() - time_req;
	dd = g->runDiffussion(rsc, RS, tsc, TS);
        cout<<dd->rac+rsc<<","<<dd->tac+tsc<<","<<(float)time_req/CLOCKS_PER_SEC<<endl;
	delete dd;
	delete [] TS;
    }
    if(runpctrid || runall)
    {
	time_req = clock();
	tsc = tsc1;
	TS = g->contrIdProxy(rsc, RS, tsc);
	time_req = clock() - time_req;
	tsc = adjustTS(tsc, TS);
	DiffusionData *dd = g->runDiffussion(rsc, RS, tsc, TS);
	cout<<argv[1]<<","<<argv[2]<<","<<argv[3]<<","<<argv[4]<<","<<argv[5]<<","<<argv[6]<<","<<"proxycontrid"<<","<<dd->rac+rsc<<","<<dd->tac+tsc<<","<<(float)time_req/CLOCKS_PER_SEC<<",";
	delete dd;

	time_req = clock();
	TS = g->pruning(rsc, RS, tsc, TS);
	time_req = clock() - time_req;
	dd = g->runDiffussion(rsc, RS, tsc, TS);
        cout<<dd->rac+rsc<<","<<dd->tac+tsc<<","<<(float)time_req/CLOCKS_PER_SEC<<endl;
	delete dd;
	delete [] TS;
    }
    if(runcen || runall)
    {
	time_req = clock();
	tsc = tsc1;
	TS = g->topkCentrality(argv[8], tsc, rsc, RS);
	time_req = clock() - time_req;
	tsc = adjustTS(tsc, TS);
	DiffusionData *dd = g->runDiffussion(rsc, RS, tsc, TS);
	cout<<argv[1]<<","<<argv[2]<<","<<argv[3]<<","<<argv[4]<<","<<argv[5]<<","<<argv[6]<<","<<"centrality"<<","<<dd->rac+rsc<<","<<dd->tac+tsc<<","<<(float)time_req/CLOCKS_PER_SEC<<",";
	delete dd;

	time_req = clock();
	TS = g->pruning(rsc, RS, tsc, TS);
	time_req = clock() - time_req;
	dd = g->runDiffussion(rsc, RS, tsc, TS);
	cout<<dd->rac+rsc<<","<<dd->tac+tsc<<","<<(float)time_req/CLOCKS_PER_SEC<<","<<argv[8]<<endl;
	delete dd;
	delete [] TS;
    }
    if(runpcen || runall)
    {
	time_req = clock();
	tsc = tsc1;
	TS = g->topkCentralityProxy(argv[8], tsc, rsc, RS);
	time_req = clock() - time_req;
	tsc = adjustTS(tsc, TS);
	DiffusionData *dd = g->runDiffussion(rsc, RS, tsc, TS);
	cout<<argv[1]<<","<<argv[2]<<","<<argv[3]<<","<<argv[4]<<","<<argv[5]<<","<<argv[6]<<","<<"proxycentrality"<<","<<dd->rac+rsc<<","<<dd->tac+tsc<<","<<(float)time_req/CLOCKS_PER_SEC<<",";
	delete dd;

	time_req = clock();
	TS = g->pruning(rsc, RS, tsc, TS);
	time_req = clock() - time_req;
	dd = g->runDiffussion(rsc, RS, tsc, TS);
	cout<<dd->rac+rsc<<","<<dd->tac+tsc<<","<<(float)time_req/CLOCKS_PER_SEC<<","<<argv[8]<<endl;
	delete dd;
	delete [] TS;
    }
    /*if(runecen || runall)
    {
	time_req = clock();
	tsc = tsc1;
	TS = g->topkECentrality(argv[9], rsc, RS, tsc);
	time_req = clock() - time_req;
	tsc = adjustTS(tsc, TS);
	DiffusionData *dd = g->runDiffussion(rsc, RS, tsc, TS);
	cout<<argv[1]<<","<<argv[2]<<","<<argv[3]<<","<<argv[4]<<","<<argv[5]<<","<<argv[6]<<","<<"ecentrality"<<","<<dd->rac+rsc<<","<<dd->tac+tsc<<","<<(float)time_req/CLOCKS_PER_SEC<<",";
	delete dd;

	time_req = clock();
	TS = g->pruning(rsc, RS, tsc, TS);
	time_req = clock() - time_req;
	dd = g->runDiffussion(rsc, RS, tsc, TS);
	cout<<dd->rac+rsc<<","<<dd->tac+tsc<<","<<(float)time_req/CLOCKS_PER_SEC<<","<<argv[9]<<endl;
	delete dd;
	delete [] TS;
    }*/
    if(runga || runall)
    {
	time_req = clock();
	tsc = tsc1;
	TS = g->GA(argv[9], rsc, RS, tsc, atoi(argv[11]), atoi(argv[12]), atoi(argv[13]));
	time_req = clock() - time_req;
	tsc = adjustTS(tsc, TS);
	DiffusionData *dd = g->runDiffussion(rsc, RS, tsc, TS);
	cout<<argv[1]<<","<<argv[2]<<","<<argv[3]<<","<<argv[4]<<","<<argv[5]<<","<<argv[6]<<","<<"ga"<<","<<dd->rac+rsc<<","<<dd->tac+tsc<<","<<(float)time_req/CLOCKS_PER_SEC<<",";
	delete dd;

	time_req = clock();
	TS = g->pruning(rsc, RS, tsc, TS);
	time_req = clock() - time_req;
	dd = g->runDiffussion(rsc, RS, tsc, TS);
	cout<<dd->rac+rsc<<","<<dd->tac+tsc<<","<<(float)time_req/CLOCKS_PER_SEC<<","<<argv[9]<<endl;
	delete dd;
	delete [] TS;
    }
    /*if(runrandprune || runall)
    {
	time_req = clock();
	tsc = tsc1;
	TS = g->randPrune(argv[9], rsc, RS, tsc, atoi(argv[10]));
	time_req = clock() - time_req;
	tsc = adjustTS(tsc, TS);
	DiffusionData *dd = g->runDiffussion(rsc, RS, tsc, TS);
	cout<<argv[1]<<","<<argv[2]<<","<<argv[3]<<","<<argv[4]<<","<<argv[5]<<","<<argv[6]<<","<<"randpruning"<<","<<dd->rac+rsc<<","<<dd->tac+tsc<<","<<(float)time_req/CLOCKS_PER_SEC<<",,,,"<<argv[9]<<endl;
	delete dd;
	delete [] TS;
    }*/
}
