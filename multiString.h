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
    std::vector<int> position;
};

std::ostream & operator << (std::ostream &, const HSFB &);
bool desHSFBCmp (const HSFB &a, const HSFB &b);

class HSFBGenerator
{
public:
    HSFBGenerator();

    HSFBGenerator(const std::vector<std::string> &);
    void set(const std::vector<std::string> &);
    const HSFB & getFirstBlock() const;
    int getBlockCenter(const HSFB &) const;
    void saveList(const std::string &filename) const;

    static const int SFB_WIDTH = 12, HSFB_SCORE = 200;

private:
    void generateRawList();
    void removeRedundance();

    bool findHSP(const std::string &subject, int subjectPos,
                 const std::string &query, int &queryPos, int & score);


    int shortestIndex;
    std::vector<std::string> ref;
    std::list<HSFB> block;
};

#endif // MULTIALIGN_H
