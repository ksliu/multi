#ifndef MULTISTRING_H
#define MULTISTRING_H
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

class HSFBgr
{
public:
    HSFBgr();

    HSFBgr(const std::vector<std::string> &);
    void set(const std::vector<std::string> &);
    const HSFB & getFirstBlock() const;
    int getBlockCenter(const HSFB &) const;


    static const int SFB_WIDTH = 12, HSFB_SCORE = 200;

private:
    void generateRawList();
    void removeRedundance();
    void saveList(const std::string &filename) const;
    bool findHSP(const std::string &subject, int subjectPos,
                 const std::string &query, int &queryPos, int & score);


    int shortestIndex;
    std::vector<std::string> ref;
    std::list<HSFB> similarBlock;

};

#endif // MULTIALIGN_H
