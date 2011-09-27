#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cstdlib>
#include "protein.h"

#include "multiString.h"
#include "pairString.h"

#include "multiAlign.h"

using namespace std;

int main()
{
    const int N = 3;
//    string pdbFiles[N] = {"1BAB.pdb", "1HLB.pdb", "1HLM.pdb"};
    string pdbFiles[N] = { "1HLB.pdb", "1HLM.pdb", "1BAB.pdb"};
//    string xyzFiles[N] = {"1BAB.xyz", "1HLB.xyz", "1HLM.xyz"};

    protein pro[N];
    vector<string>pcl;

    for (int i=0; i<N; i++)
    {
        pro[i].loadFile(pdbFiles[i]);
        //        pro[i].toXYZFile(xyzFiles[i]);
        pcl.push_back(pro[i].cl);
    }

    MultiAlign ma(pro, N);

    ma.outputAlignResult("multiAlign1.txt");

    ma.run();

    ma.outputAlignResult("multiAlign2.txt");


    return 0;
}
