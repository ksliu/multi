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


HSFB::HSFB(int np): depth(1), score(0), position(np, -1)
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
    for (std::vector<int>::const_iterator it = hsfb.position.begin(); it != hsfb.position.end(); ++it)
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

HSFBGenerator::HSFBGenerator()
{

}

HSFBGenerator::HSFBGenerator(const std::vector<std::string> &inStr)
{
    set(inStr);
}
void HSFBGenerator::set(const std::vector<std::string> &inStr)
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
    block.sort(desHSFBCmp);
    removeRedundance();

    if (block.begin() -> depth != np)
    {
        cout << "Warning: the first block is not full." << endl;
    }
}
void HSFBGenerator::generateRawList()
{
    int np = ref.size();
    int subjectLength = ref[shortestIndex].size();
    for (int subjectPos=0; subjectPos <= subjectLength - SFB_WIDTH; subjectPos++)
    {
        HSFB hsfb(np);
        hsfb.position[shortestIndex] = subjectPos;
        for (int iQuery=0; iQuery<np; iQuery++)
        {
            if (iQuery == shortestIndex) continue;
            int queryPos, score;
            if (findHSP(ref[shortestIndex], subjectPos, ref[iQuery], queryPos, score))
            {
                hsfb.position[iQuery] = queryPos;
                hsfb.score += score;
                hsfb.depth++;
            }
        }
        if (hsfb.depth > 1)
        {
            string seqBlock[SFB_WIDTH];
            for (int i=0; i<np; i++)
            {
                int pos = hsfb.position[i];
                if (pos != -1)
                {
                    for (int j=0; j<SFB_WIDTH; j++)
                        seqBlock[j].push_back(ref[i][pos + j]);
                }
            }
            for (int i=0; i<SFB_WIDTH; i++)
                hsfb.consensus.push_back(bio::cleConsensus(seqBlock[i]));
            block.push_back(hsfb);
        }
    }
}

bool HSFBGenerator::findHSP(const string &subject, int subjectPos,
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

void HSFBGenerator::removeRedundance()
{
    int numProtein = ref.size();
    vector< vector<int> > marker(numProtein);
    for (int i=0; i<numProtein; i++)
        marker[i].assign(ref[i].size(), 0);

    std::list<HSFB>::iterator it = block.begin();
    while (it != block.end())
    {
        int overlap = 0;
        for (int i=0; i<numProtein; i++)
        {
            int pos = it->position[i];
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
            it = block.erase(it);
        else
            ++it;
    }
}
void HSFBGenerator::saveList(const string &fn) const
{
    ofstream fout(fn.c_str());
    if (!fout)
    {
        cout << "cannot open file " <<  fn << endl;
        exit(1);
    }

    for (list<HSFB>::const_iterator it = block.begin(); it !=block.end(); ++it)
        fout << *it << endl;
    fout.close();
}
int HSFBGenerator::getBlockCenter(const HSFB &hsfb) const
{
    int np = ref.size(), highestIndex = 0, highestScore = 0;
    for (int i=0; i<np; i++)
    {
        int score = 0;
        int pos = hsfb.position[i];
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
const HSFB & HSFBGenerator::getFirstBlock() const
{
    if (!block.empty())
        return  *block.begin();
    else
    {
        cout << "empty HSFB list" << endl;
        exit(1);
    }
}
