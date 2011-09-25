#include <string>
#include <vector>
#include <list>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include "bioinfo.h"
#include "multiString.h"
using namespace std;


HSFB::HSFB(int np): depth(1), score(0), positions(np, -1)
{

}

ostream & operator<< (ostream &os, const HSFB &hsfb)
{
    const char * FS = " ";
    os << "{" << FS;
    os << "depth = " << setw(4) <<  hsfb.depth << FS;
    os << "score = " << setw(4) << hsfb.score << FS;
    os << "consensus = " << setw(15) << hsfb.consensus << FS;
    os << "position = ";
    for (std::vector<int>::const_iterator it = hsfb.positions.begin(); it != hsfb.positions.end(); ++it)
        os << setw(6) << *it;
    os << FS;
    os << "}";
    return os;
}

bool desHSFBCmp (const HSFB &a, const HSFB &b)
{
    if (a.depth > b.depth)
        return true;
    if (a.depth == b.depth)
        return a.score > b.score;
    return false;
}

HSFBgr::HSFBgr()
{

}

HSFBgr::HSFBgr(const std::vector<std::string> &inStr)
{
    set(inStr);
}
void HSFBgr::set(const std::vector<std::string> &inStr)
{
    shortestIndex = 0;
    ref = inStr;

    if(ref.size() < 2)
    {
        cout << "Too few proteins" << endl;
        exit(1);
    }
    int np = ref.size();
    for (int i=0; i < np; i++)
    {
        if (ref[i].size() < ref[shortestIndex].size())
            shortestIndex = i;
    }

    generateRawList();
    saveList("m1.txt");

    similarBlock.sort(desHSFBCmp);
    saveList("m2.txt");

    removeRedundance();
    saveList("m3.txt");

    if (similarBlock.begin() -> depth != np)
    {
        cout << "error: the first block is not full." << endl;
        exit(1);
    }
//    cout << "shortest: "<< shortestIndex << endl;
//    cout << "center: "<<  getBlockCenter(*similarBlock.begin()) << endl;
}
void HSFBgr::generateRawList()
{
    int np = ref.size();
    int subjectLength = ref[shortestIndex].size();
    for (int subjectPos=0; subjectPos <= subjectLength - SFB_WIDTH; subjectPos++)
    {
        HSFB hsfb(np);
        hsfb.positions[shortestIndex] = subjectPos;
        for (int iQuery=0; iQuery<np; iQuery++)
        {
            if (iQuery == shortestIndex) continue;
            int queryPos, score;
            if (findHSP(ref[shortestIndex], subjectPos, ref[iQuery], queryPos, score))
            {
                hsfb.positions[iQuery] = queryPos;
                hsfb.score += score;
                hsfb.depth++;
            }
        }
        if (hsfb.depth > 1)
        {
            string seqBlock[SFB_WIDTH];
            for (int i=0; i<np; i++)
            {
                int pos = hsfb.positions[i];
                if (pos != -1)
                {
                    for (int j=0; j<SFB_WIDTH; j++)
                        seqBlock[j].push_back(ref[i][pos + j]);
                }
            }
            for (int i=0; i<SFB_WIDTH; i++)
                hsfb.consensus.push_back(bio::cleConsensus(seqBlock[i]));
            similarBlock.push_back(hsfb);
        }
    }
}

bool HSFBgr::findHSP(const string &subject, int subjectPos,
                     const string &query, int &queryPos, int & score)
{
    queryPos = -1;
    score = -1000*SFB_WIDTH;

    int i, j, s, queryLength = query.size();

    for (i = 0; i <= queryLength - SFB_WIDTH; i++)
    {
        for (s=j=0; j< SFB_WIDTH; j++)
            s += bio::cleScore(subject[subjectPos + j], query[i + j]);
        if (s > score)
        {
            score = s;
            queryPos = i;
        }
    }
    return (score >= HSFB_SCORE);
}

void HSFBgr::removeRedundance()
{
    int numProtein = ref.size();
    vector< vector<int> > marker(numProtein);
    for (int i=0; i<numProtein; i++)
        marker[i].assign(ref[i].size(), 0);

    list<HSFB>::iterator it = similarBlock.begin();
    while (it != similarBlock.end())
    {
        int overlap = 0;
        for (int i=0; i<numProtein; i++)
        {
            int pos = it->positions[i];
            if (pos == -1) continue;
            for (int j=0; j<SFB_WIDTH; j++)
            {
                if (marker[i][pos+j] == 0)
                    marker[i][pos+j] = 1;
                else
                    overlap++;
            }
        }
        if ( overlap > 0.5 * it->depth * SFB_WIDTH )
            it = similarBlock.erase(it);
        else
            ++it;
    }
}
void HSFBgr::saveList(const string &fn) const
{
    ofstream fout(fn.c_str());
    if (!fout)
    {
        cout << "cannot open file " <<  fn << endl;
        exit(1);
    }

    for (list<HSFB>::const_iterator it = similarBlock.begin(); it !=similarBlock.end(); ++it)
        fout << *it << endl;
    fout.close();
}
int HSFBgr::getBlockCenter(const HSFB &hsfb) const
{
    int np = ref.size(), highestIndex = 0, highestScore = 0;
    for (int i=0; i<np; i++)
    {
        int score = 0;
        int pos = hsfb.positions[i];
        if (pos != -1)
        {
            for (int j=0; j<SFB_WIDTH; j++)
                score += bio::cleScore(ref[i][pos + j], hsfb.consensus[j]);
            if (score > highestScore)
            {
                highestIndex = i;
                highestScore = score;
            }
        }
    }
    return highestIndex;
}
const HSFB & HSFBgr::getFirstBlock() const
{
    if (!similarBlock.empty())
        return  *similarBlock.begin();
    else
    {
        cout << "empty HSFB list" << endl;
        exit(1);
    }
}
