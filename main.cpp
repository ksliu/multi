#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cstdlib>
#include "protein.h"
#include "multiAlign.h"

using namespace std;

int main()
{
    const int N = 3;
    string fn[N] = {"d1ufaa2.pdb", "d1z7aa1.pdb", "d2cc0a1.pdb"};
    protein p [N];

    for (int i=0; i<N; i++)
    {
        p[i].loadFile(fn[i]);
//        p[i].toXYZFile(fn[i]+".xyz");
    }


    foo(p, N);

    return 0;
}


//int main(int argc, char *argv[])
//{
//    if (argc != 3)
//    {
//        cerr << "Usage: multi-ali input-file output-file" << endl;
//        exit(1);
//    }

//    return 0;
//}

//void readInput(const std::string & fn)
//{
//    ifstream fs(fn.c_str());
//    if (!fs)
//    {
//        cerr << "cannot open file " << fn << endl;
//        exit(1);
//    }
//    string line;
//    int n=0;
//    while (getline(fs, line))
//    {
//        cout << line << endl;
//    }
//    fs.close();
//}
