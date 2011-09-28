#include <fstream>
#include <iostream>
#include <iomanip>
#include <cstdlib>

#include "protein.h"
#include "pdb.h"
#include "gencle.h"
#include "bioinfo.h"
protein::protein()
{

}

protein::protein(const std::string &fn, char c)
{
    loadFile(fn, c);
}

void protein::loadFile(const std::string & fn, char c)
{
    id = fn;
    chain = c;
    PDB pdb(fn, c);
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
void protein::genFakePDB(char chain, int startAtomNo, std::ostream &os)
{
    using std::setw;
    using std::endl;
    for (int i=0; i< ca.len(); i++)
    {
        os << "ATOM  " << setw(5) << startAtomNo + i << " " << " CA " << " " << setw(3)
           << bio::aaConvert13(aa[i]) << " " << chain << setw(4) << i
           << " " << "   ";
        for (int j = 0; j < 3; j++)
            os << setw(8) << ca[i][j];
        os << endl;
    }
    os << "TER   " << endl;
    for (int i = 1; i < ca.len() ; i++)
        os << "CONECT" << setw(5) << startAtomNo + i - 1 << setw(5) << startAtomNo +i << endl;

}
