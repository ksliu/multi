#include <string>
#include <cstdlib>
#include "protein.h"
#include "multiAlign.h"

using namespace std;

void test();
bool RUN_TEST = true;

int main(int argc, char **argv)
{
    if (RUN_TEST)
    {
        test();
        return 0;
    }

    if (argc < 6)
    {
        cout << "program -o output_file pdb_file1 pdb_file2 pdb_file3..." << endl;
        exit(1);
    }

    string outputFile = argv[2];
    int numPDB = argc - 3;

    if (numPDB < 3)
    {
        cout << "input at least 3 PDB files" << endl;
        exit(1);
    }

    protein * pro = new protein [numPDB];
    for (int i=0; i <numPDB; i++)
        pro[i].loadFile(argv[i+3]);


    MultiAlign ma (pro, numPDB);
    ma.run();
    ma.outputAlignResult(outputFile);

    delete [] pro;

    return 0;
}

void test()
{
    const int N = 3;
    //    string pdbFiles[N] = {"1BAB.pdb", "1HLB.pdb", "1HLM.pdb"};
    string pdbFiles[N] = {"1HLB.pdb", "1HLM.pdb", "1BAB.pdb"};

    protein pro[N];
    for (int i=0; i<N; i++)
    {
        pro[i].loadFile(pdbFiles[i]);
        //        pro[i].toXYZFile(xyzFiles[i]);
    }

    MultiAlign ma(pro, N);
    ma.run();
    ma.outputAlignResult("multiAlign.txt");
}
