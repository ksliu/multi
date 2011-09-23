#include <iostream>
#include <iomanip>
#include <fstream>
#include <algorithm>
#include <vector>
#include <string>
#include <cmath>
#include <cstdlib>
#include "bioinfo.h"
#include "pairAlign.h"
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

AFPCandidate::AFPCandidate(const std::string &a, const std::string &b)
{
    pa = a;
    pb = b;

    rawList(SW, SC);
    saveList("p1.txt");

    removeRedundance(SW / 2);
    removeRedundance(SW / 2);
    saveList("p2.txt");

    sort(sq.begin(), sq.end(), deSFPcmp);
    saveList("p3.txt");
}

void AFPCandidate::removeRedundance(int d)
{
    vector<SFP>::size_type i, j, ix, nn;
    int t;

    for (nn = i = 0; i != sq.size(); i = j)
    {
        t = sq[i].score;
        ix = i;
        for (j = i + 1; j != sq.size(); ++j)
        {
            if (sq[j].ia - sq[i].ia >= d || (sq[i].ia - sq[i].ib != sq[j].ia - sq[j].ib))
                break; // get proper j

            if (sq[j].score > t) // find maximum in sq[*].score
            {
                t = sq[j].score;
                ix = j;
            }
        }
        sq[nn++] = sq[ix];
    }
    sq.erase(sq.begin() + nn, sq.end());
}


void AFPCandidate::rawList(int width, int lower_socre)
{
    int i, j, k, m, n, s, la = pa.size(), lb = pb.size();
    vector<int> z(min(la, lb) + 1);

    if (la < width || lb < width)
    {
        cerr << "code too short!" << endl;
        exit(1);
    }
    // stand on 'b', shift 'a'
    // i, k ~ a, j ~ b, fixed
    for (k = 0; k <= la - width; k++)
    {
        z[0] = 0;
        for (i = k, j = 0; i < la && j < lb; i++, j++)
            z[j + 1] = z[j] + bio::cleScore(pa[i], pb[j]);

        for (m = 0, n = width; n <= j; m++, n++)
        {
            s = z[n] - z[m];
            if (s > lower_socre)
                sq.push_back(SFP(m + k, m, s));
        }
    }
    // stand on 'a', shift 'b'. Set k=1 to exclude head-head
    // i, k ~ b, j ~ a, fixed
    for (k = 1; k <= lb - width; k++)
    {
        z[0] = 0;
        for (i = k, j = 0; i < lb && j < la; i++, j++)
            z[j + 1] = z[j] + bio::cleScore(pb[i], pa[j]);

        for (m = 0, n = width; n <= j; m++, n++)
        {
            s = z[n] - z[m];
            if (s > lower_socre)
                sq.push_back(SFP(m, m + k, s));
        }
    }
}

void AFPCandidate::saveList(const string &fn) const
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


//void foo()
//{
//    list<SFP> sq;
//    list<SFP>::iterator it = sq.begin();
//    while (it != sq.end())
//    {
//        int t = it->score;
//        list<SFP>::iterator r = advance(it, 1);
//        while (r != sq.end())
//        {
//            if (it->ia - r->ia >=d  || (it->ia - it->ib  !=  r->ia - r->ib))
//                break;
//            if ( r->score > t)
//            {
//                itr = r;
//            }
//            ++r;
//        }
//        swap(it, r);
//        it = r;
//    }
//}


