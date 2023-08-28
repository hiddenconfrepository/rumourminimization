#include <list>
#include <iostream>
using namespace std;
//Class to output of the diffusion
class DiffusionData
{
    public:
    	int rac, tac; // rac/tac excluding rumour/truth seed count	
	list <int> *RA, *TA;
        DiffusionData();
	~DiffusionData();
};

class popwrc
{
    public:
    	int *pop;
	int rc;
	popwrc(int k);
	~popwrc();
};

