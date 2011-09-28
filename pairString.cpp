#include <iostream>
#include <iomanip>
#include <fstream>
#include <algorithm>
#include <vector>
#include <string>
#include <cmath>
#include <cstdlib>
#include "bioinfo.h"
#include "pairString.h"
using namespace std;

SFP::SFP(int a, int b, int c) : ia(a), ib(b), score(c)
{

}
std::ostream & operator << (std::ostream &os, const SFP &sfp)
{
    const char * FS = " ";
    os << "{" << FS;
    os << "score = " << setw(5) << sfp.score << FS;
    os << "first = " << setw(5) << sfp.ia << FS;
    os << "second = " << setw(5) << sfp.ib;
    os << FS;
    os << "}";
    return os;
}

bool deSFPcmp(const SFP & a, const SFP & b)
{
    return a.score > b.score;
}

SFPGenerator::SFPGenerator()
{

}

SFPGenerator::SFPGenerator(const std::string &a, const std::string &b)
{
    set(a, b);
}
void SFPGenerator::set(const std::string &a, const std::string &b)
{
    pa = a;
    pb = b;

    generateRawList(SW, SC);
    removeRedundance(SW / 2);
    removeRedundance(SW / 2);
    sort(sq.begin(), sq.end(), deSFPcmp);
}
const std::vector<SFP> & SFPGenerator::getCandidate() const
{
    return sq;
}

void SFPGenerator::generateRawList(int width, int lower_socre)
{
    int la = pa.size(), lb = pb.size();
    vector<int> z(min(la, lb) + 1);

    if (la < width || lb < width)
    {
        cerr << "CLE code too short!" << endl;
        exit(1);
    }

    int i, j, k, m, n, s;
    for (k = 0; k <= la - width; k++)
    {
        z[0] = 0;
        for (i = k, j = 0; i < la && j < lb; i++, j++)
            z[j + 1] = z[j] + bio::cleScore(pa[i], pb[j]);

        for (m = 0, n = width; n <= j; m++, n++)
        {
            s = z[n] - z[m];
            if (s >= lower_socre)
                sq.push_back(SFP(m + k, m, s));
        }
    }
    for (k = 1; k <= lb - width; k++)
    {
        z[0] = 0;
        for (i = k, j = 0; i < lb && j < la; i++, j++)
            z[j + 1] = z[j] + bio::cleScore(pb[i], pa[j]);

        for (m = 0, n = width; n <= j; m++, n++)
        {
            s = z[n] - z[m];
            if (s >= lower_socre)
                sq.push_back(SFP(m, m + k, s));
        }
    }
}

void SFPGenerator::removeRedundance(int d)
{
    vector<SFP>::size_type i, j, ix, ns;

    for (ns = i = 0; i != sq.size(); i = j)
    {
        int t = sq[i].score;
        ix = i;
        for (j = i + 1; j != sq.size(); ++j)
        {
            if (sq[j].ia - sq[i].ia >= d || (sq[i].ia - sq[i].ib != sq[j].ia - sq[j].ib))
                break;
            if (sq[j].score > t)
            {
                t = sq[j].score;
                ix = j;
            }
        }
        sq[ns++] = sq[ix];
    }
    sq.erase(sq.begin() + ns, sq.end());
}

void SFPGenerator::saveList(const string &fn) const
{
    ofstream fout(fn.c_str());
    if (!fout)
    {
        cout << "cannot open file " <<  fn << endl;
        exit(1);
    }

    for (vector<SFP>::const_iterator it = sq.begin(); it !=sq.end(); ++it)
        fout << *it << endl;
    fout.close();
}


