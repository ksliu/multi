#include <iostream>
#include <iomanip>
#include <fstream>
#include <algorithm>
#include <vector>
#include <list>
#include <cassert>
#include <string>
#include <cmath>
#include <cstdlib>
#include "hfali.h"
#include "bioinfo.h"
#include "read_conf.h"
#include "multiAlign.h"
using namespace std;

const int SFP_LEN = 12, SFP_SCORE = 200;

class HSFB
{
public:
    HSFB(int n);
    int depth;
    int score;
    std::vector<int> positions;
    std::string consensus;
    static  std::vector<int> length;
    friend  ostream & operator << (ostream &os, const HSFB &) ;
};
std::vector<int> HSFB::length;

bool desHSFBCmp (const HSFB &a, const HSFB &b)
{
    if ( a.depth > b.depth )
        return true;
    if (a.depth == b.depth)
        return a.score > b.score;
    return false;
}

HSFB::HSFB(int n): depth(1), score(0), positions(n, -1)
{

}
ostream & operator<< (ostream &os, const HSFB &h)
{
    os << "{" << endl;
    os << "depth = " << h.depth << endl;
    os << "score = " << h.score << endl;
    os << "consensus = " << h.consensus << endl;
    os << "position = ";
    for (std::vector<int>::const_iterator it = h.positions.begin();
         it != h.positions.end(); ++it)
    {
        os << *it << " ";
    }
    os << endl;
    os << "}";

    return os;
}



bool findHSP(const string & subject, int subjectPos, const string &query, int &queryPos, int & score)
{
    queryPos = -1;
    score = -1000*SFP_LEN;

    int i, j, s, queryLength = query.size();

    for (i = 0;  i <= queryLength - SFP_LEN; i++)
    {
        for (s =j=0; j< SFP_LEN; j++)
            s += bio::cleScore( subject[subjectPos + j], query[i + j]);
        if (s > score)
        {
            score = s;
            queryPos = i;
        }
    }
    return (score >= SFP_SCORE);
}
void shaveHSFBList(list<HSFB> & pre)
{
    int numProtein = HSFB::length.size();
    vector< vector<int> > marker(numProtein);
    for (int i=0; i<numProtein; i++)
        marker[i].resize(HSFB::length[i]);
    for (int i=0; i<numProtein; i++)
        for (int j=0; j < HSFB::length[i]; j++)
            marker[i][j] = 0;

    list<HSFB>::iterator it = pre.begin();
    while (it != pre.end())
    {
        int overlap = 0;
        for (int i=0; i<numProtein; i++)
        {
            int pos = it->positions[i];
            if (pos == -1) continue;
            for (int j=0; j<SFP_LEN; j++)
            {
                if (marker[i][pos+j] == 0)
                    marker[i][pos+j] = 1;
                else
                    overlap++;
            }
        }
        if ( overlap > 0.5 * it->depth * SFP_LEN )
            it = pre.erase(it);
        else
            ++it;
    }
}

void foo(protein *p, int numProtein) // delete hsfb, after call
{
    assert(numProtein > 2);
    int iMin = 0;
    for (int i=0; i<numProtein; i++)
    {
        if ( p[i].ca.len() < p[iMin].ca.len() )
            iMin = i;
    }

    protein **pp = new protein * [numProtein];
    for (int i=0; i<numProtein; i++)
        pp[i] = p+i;

    pp[0] = p + iMin;
    pp[iMin] = p;

    for (int i=0; i<numProtein; i++)
        HSFB::length.push_back(pp[i]->ca.len());


    list<HSFB> hData;
    int subjectLength = pp[0]->ca.len();
    for (int subjectPos=0; subjectPos <= subjectLength - SFP_LEN; subjectPos++)
    {
        HSFB hsfb(numProtein);
        hsfb.positions[0] = subjectPos;
        for (int iQuery=1; iQuery<numProtein; iQuery++)
        {
            int queryPos, score;
            if ( findHSP(pp[0]->cl, subjectPos, pp[iQuery]->cl, queryPos, score) )
            {
                hsfb.positions[iQuery] = queryPos;
                hsfb.score += score;
                hsfb.depth++;
            }
        }
        if (hsfb.depth > 0)
        {
            string candidates[SFP_LEN];
            for (int i=0; i<SFP_LEN; i++)
                candidates[i].push_back(pp[0]->cl[subjectPos + i]);
            for (int j=1; j<numProtein; j++)
            {
                int pos = hsfb.positions[j];
                if (pos != -1)
                {
                    for (int i=0; i<SFP_LEN; i++)
                        candidates[i].push_back(pp[j]->cl[pos + i]);
                }
            }
            for (int i=0; i<SFP_LEN; i++)
                hsfb.consensus.push_back(bio::cleConsensus(candidates[i]));
            hData.push_back(hsfb);
        }
    }

    // debug
    // cout << "hDat.size() " << hData.size() << endl;
    //    cout << hData[0] << endl;
    hData.sort(desHSFBCmp);



    ofstream ftmp1("tmp1.txt");
    ofstream ftmp2("tmp2.txt");

    for (list<HSFB>::const_iterator it = hData.begin(); it !=hData.end(); ++it)
        ftmp1 << *it << endl;

    shaveHSFBList(hData);

    for (list<HSFB>::const_iterator it = hData.begin(); it !=hData.end(); ++it)
        ftmp2 << *it << endl;

    ftmp1.close();
    ftmp2.close();

    delete []pp;
}