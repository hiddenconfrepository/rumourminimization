#include <list>
#include <iostream>
using namespace std;
#include "DiffusionData.h"

DiffusionData::DiffusionData()
{
    RA = new list<int>();
    TA = new list<int>();
    rac = 0;
    tac = 0;
}

DiffusionData::~DiffusionData()
{
    RA->clear();
    delete RA;
    TA->clear();
    delete TA;
}

popwrc::popwrc(int k)
{
    pop = new int[k];
    rc = 0;
}
popwrc::~popwrc()
{
    delete [] pop;
}
