#include  <iostream>
#include  <fstream>
#include  <cstdlib>
#include "read_conf.h"
using namespace std;
//"/a.home2/code/z.hform.alignment/src/hform_align.conf";
//"/home2/dbu/janic/protein/tools/m.hform.alignment/src/hform_align.conf";


Para * Para::data = 0;

Para * Para::instance(const char *s)
{
    if (!data)
        data = new Para(s);
    return data;
}
Para * Para::getPointer()
{
    return data;
}
void Para::destorySelf()
{
    delete data;
}
Para::Para(const char *s)
{
    ifstream in(s);
    if (in.fail())
    {
        cerr << "cannot open configure file: " << s << endl;
        exit(1);
    }
    string line;
    while (getline(in, line))
    {
        if (line.empty() || line[0] == '#')
            continue;
        string::size_type found = line.find('=');
        if (found == string::npos || found + 1 == line.size())
        {
            cerr << "configure file format error" << endl;
            exit(1);
        }
        string key;
        for (string::size_type i = 0; i < found; ++i)
        {
            if (line[i] != ' ')
                key.push_back(toupper(line[i]));
        }

        double val = atof(line.substr(found + 1).c_str());
        if (val < 0 || val > 10)
        {
            cerr << "configure file value domain error" << val << endl;
            exit(1);
        }
        dict[key] = val;
    }
    in.close();
}

double Para::get_value(const string &key) const
{
    string my_key;
    for (string::size_type i = 0; i != key.size(); ++i)
    {
        my_key.push_back(toupper(key[i]));
    }
    map<string, double>::const_iterator it = dict.find(my_key);
    if (it == dict.end())
    {
        cerr << key << " is not a valid key" << endl;
        exit(1);
    }
    return it->second;
}
