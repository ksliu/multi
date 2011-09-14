#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cstdlib>
#include "hfali.h"
#include "read_conf.h"
using namespace std;

void get_int_pair(const char *s, int *a, int *b);

int main(int argc, char *argv[])
{
    if (argc < 5)
    {
        cerr << "Usage: hfali pdb-file1 pdb-file2 pdb1-seg1,pdb1-seg2 pdb2-seg1,pdb2-seg2 [-f configure-file -l 0|1 -o rasmol-script-file]"
             << endl;
        cerr << "For example: hfali pdb-file1 pdb-file2 4-9,10-20 8-13,15-25 -l 0 -o rasmol-script-file "
             << endl;
        exit(1);
    }
    const char *subject = argv[1], *query = argv[2];
    typedef int rangeT[2];
    rangeT p11, p12, p21, p22;
    get_int_pair(argv[3], p11, p12);
    get_int_pair(argv[4], p21, p22);

    int level = 0;
    const char * script_file = 0;

    // default: in the working directory that contains the binary executive file
    const char *configure_file = "hform_align.conf";

    if (argc > 5)
    {
        int i = 5;
        while (argv[i] && argv[i + 1])
        {
            if (argv[i][0] == '-')
            {
                if (argv[i][1] == 'l')
                    level = atoi(argv[i + 1]);
                else if (argv[i][1] == 'o')
                    script_file = argv[i + 1];
                else if (argv[i][1] == 'f')
                    configure_file = argv[i + 1];
            }
            i++;
        }
    }

    Para::instance(configure_file);

    hfali ob(subject, query);

    if (level == 0)
    {
        ob.solve(p11, p12, p21, p22, hfali::IM_once);
    }
    else if (level == 1)
    {
        ob.solve(p11, p12, p21, p22, hfali::IM_zoom_in);
    }
    else
    {
        cerr << "only support level 0 and level 1" << endl;
        exit(1);
    }
    ob.output();

    if (script_file)
    {
        ob.output_script(script_file);
    }

    Para::destorySelf();
    return 0;
}
void get_int_pair(const char *s, int *a, int *b)
{
    string p = s;
    string::size_type comma = p.find(',');
    if (comma == string::npos || comma + 1 == p.size())
    {
        cerr << "range from command line error" << endl;
        exit(1);
    }
    string segment[2] =
    { p.substr(0, comma), p.substr(comma + 1) };

    int *range = a;
    for (int i = 0; i < 2; i++)
    {
        string::size_type dash = segment[i].find('-');
        if (dash == string::npos || dash + 1 == segment[i].size())
        {
            cerr << "range from command line error" << endl;
            exit(1);
        }
        if (i == 1)
            range = b;
        range[0] = atoi(segment[i].substr(0, dash).c_str());
        range[1] = atoi(segment[i].substr(dash + 1).c_str());
    }
}
