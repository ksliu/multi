#ifndef MULTIALIGN_H
#define MULTIALIGN_H
#include <string>
#include <vector>
#include <iostream>
#include "protein.h"

class HSFB
{
public:
    HSFB(int np);
    int depth;
    int score;
    std::string consensus;
    std::vector<int> positions;

    static  std::vector<protein *> ref;
    friend  std::ostream & operator << (std::ostream &, const HSFB &) ;
};


void foo(protein *p, int numProtein);
#endif // MULTIALIGN_H
