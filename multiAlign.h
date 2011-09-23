#ifndef MULTIALIGN_H
#define MULTIALIGN_H
#include <string>
#include <vector>
#include <list>
#include <iostream>

class HSFB
{
public:
    HSFB(int np);
    int depth;
    int score;
    std::string consensus;
    std::vector<int> positions;
};

std::ostream & operator << (std::ostream &, const HSFB &);
bool desHSFBCmp (const HSFB &a, const HSFB &b);

class MABCandidate
{
public:
    MABCandidate(const std::vector<std::string> &);

private:
    void rawList();
    void shaveRedundance();
    void saveList(const std::string &filename) const;
    bool findHSP(const std::string &subject, int subjectPos,
                 const std::string &query, int &queryPos, int & score);

    std::vector<std::string> ref;
    std::list<HSFB> similarBlock;
    static const int SFB_WIDTH = 12, HSFB_SCORE = 200;
};


#endif // MULTIALIGN_H
