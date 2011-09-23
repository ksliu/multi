#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cstdlib>
#include "protein.h"

#include "multiString.h"
#include "pairString.h"

using namespace std;

int main()
{
    const int N = 3;
    string pdbFiles[N] = {"d1ufaa2.pdb", "d1z7aa1.pdb", "d2cc0a1.pdb"};
    //    string xyzFiles[N] = {"d1ufaa2.xyz", "d1z7aa1.xyz", "d2cc0a1.xyz"};
    protein pro[N];
    vector<string>pcl;

    for (int i=0; i<N; i++)
    {
        pro[i].loadFile(pdbFiles[i]);
        // p[i].toXYZFile(xyzFiles[i]);
        pcl.push_back(pro[i].cl);
    }


    HSFBgr hsf(pcl);
    SFPgr afp(pcl[0], pcl[1]);

    return 0;
}
