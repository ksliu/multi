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


HSFBgr::HSFBgr(const std::vector<std::string> &inStr)
    :shortestIndex(0), ref(inStr)
{
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
            vector<string> seqBlock;
            for (int j=0; j<np; j++)
            {
                int pos = hsfb.positions[j];
                if (pos != -1)
                    seqBlock.push_back(ref[j].substr(pos, SFB_WIDTH));
            }
            blockConsensus(seqBlock, hsfb.consensus);
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
void HSFBgr::blockConsensus(const vector<string> &block, string & consensus)
{
    if (block.size() < 2)
    {
        cerr << "too few input" << endl;
        exit(1);
    }
    int i, row = block.size();
    for (i=1; i<row; i++)
    {
        if (block[i].size() != block[0].size())
        {
            cerr << "length mismatch" << endl;
            exit(1);
        }
    }

    int j, column=block[0].size();
    string columnStr;
    columnStr.resize(row);
    consensus.resize(column);

    for (j=0; j<column; j++)
    {
        for (i=0; i<row; i++)
            columnStr[i] = block[i][j];
        consensus[j] = bio::cleConsensus(columnStr);
    }
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
