#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cstdlib>
#include "protein.h"
#include "multiAlign.h"

using namespace std;

vector< pair<string, char> > getProtein(const string &fn);

int main(int argc, char *argv[])
{
    if (argc != 3)
    {
        cout << "Usage:\n prog input_file output_file_prefix" << endl;
        exit(1);
    }
    vector< pair<string, char> > input =  getProtein(argv[1]);
    int n = input.size();
    protein *pro = new protein[n];
    for (int i=0; i<n; i++)
    {
        pro[i].loadFile(input[i].first, input[i].second);
    }

    MultiAlign ma(pro, n);
    ma.run();

    string outputPrefix = argv[2];
    ma.outputAlignResult(outputPrefix + ".txt");
    ma.genScript(outputPrefix + ".script");

    delete [] pro;
}

vector< pair<string, char> > getProtein (const string &fn)
{
    ifstream fin(fn.c_str());

    if (fin.fail())
    {
        cout << "cannot open file : " << fn << endl;
        exit(1);
    }
    vector< pair<string, char> > ret;

    string line;
    while (getline (fin, line))
    {
        stringstream ss (line);
        string file;
        char chain;
        ss >> file >> chain;
        ret.push_back(pair<string, char>(file, chain));

        cout << file << chain << endl;

    }
    fin.close();
    return ret;
}
