#ifndef PAIRALIGN_H
#define PAIRALIGN_H

#include <vector>
#include <string>
#include <iostream>


class SFP
{
public:
    int ia, ib, score;
    SFP(int a = 0, int b = 0, int c = 0);
};

std::ostream & operator << (std::ostream &, const SFP &);

bool deSFPcmp(const SFP & a, const SFP & b);

class AFPCandidate
{
public:
    AFPCandidate(const std::string &a, const std::string &b);

private:
    void rawList(int width, int lower_socre);
    void removeRedundance(int);
    void saveList(const std::string &filename) const;

    std::string pa, pb;
    std::vector<SFP> sq;
    static const int SW = 8, SC = 100;
};

#endif // PAIRALIGN_H
