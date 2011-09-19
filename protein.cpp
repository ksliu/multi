#include <fstream>
#include <iostream>
#include <iomanip>
#include <cstdlib>

#include "protein.h"
#include "pdb.h"
#include "gencle.h"
protein::protein()
{

}

protein::protein(const std::string &fn)
{
    loadFile(fn);
}

void protein::loadFile(const std::string & fn)
{
    id = fn;
    PDB pdb(fn);
    pdb.getAA(aa);
    pdb.getCa(ca);
    cl = bio::ca2cl(ca.ptr(), ca.len());
}
void protein::toXYZFile(const std::string &fn)
{
    using namespace std;
    ofstream fs (fn.c_str());
    if (!fs)
    {
        cerr << "cannot open file " << fn << endl;
        exit(1);
    }
    fs << fixed << setprecision(3);
    int n = ca.len();
    for (int i=0; i<n; i++)
    {
        fs << aa[i] << cl[i];
        for (int j =0; j<3; j++)
            fs << setw(8) << ca[i][j];
        fs << endl;
    }
    fs.close();
}
