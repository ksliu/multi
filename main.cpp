#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cstdlib>
#include "protein.h"

#include "multiAlign.h"
#include "pairAlign.h"

using namespace std;

int main()
{
    const int N = 3;
    string pdbFiles[N] = {"d1ufaa2.pdb", "d1z7aa1.pdb", "d2cc0a1.pdb"};
    //    string xyzFiles[N] = {"d1ufaa2.xyz", "d1z7aa1.xyz", "d2cc0a1.xyz"};
    protein p[N];
    vector<string>pcl;

    for (int i=0; i<N; i++)
    {
        p[i].loadFile(pdbFiles[i]);
        // p[i].toXYZFile(xyzFiles[i]);
        pcl.push_back(p[i].cl);
    }


//    MABCandidate mabc(pcl);
    AFPCandidate afp(pcl[0], pcl[1]);

    return 0;
}
