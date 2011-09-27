#ifndef PAIRSTRING_H
#define PAIRSTRING_H

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

class SFPGenerator
{
public:
    SFPGenerator();
    SFPGenerator(const std::string &a, const std::string &b);
    void set(const std::string &a, const std::string &b);

    const std::vector<SFP> & getCandidate() const;
    static const int SW = 8, SC = 100;

private:
    void generateRawList(int width, int lower_socre);
    void removeRedundance(int);
    void saveList(const std::string &filename) const;


    std::string pa, pb;
    std::vector<SFP> sq;
};

#endif
