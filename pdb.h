#ifndef PDB_H
#define PDB_H

#include <vector>
#include <string>
#include "dar.h"

class PDB
{
public:
    PDB(const std::string &file, char chain='_');

    void getAA(std::string &) const;
    void getCa(Dar &) const;
    void getCaCb(Dar &, Dar &) const;
    void getBB(Dar bb[5]) const;

private:
    void cacheFirst(const std::string &);
    void cacheRange(const std::string &, char);

    int getCaNumber() const;
    int bbAtomType(const std::string &) const;

    std::vector<std::string> domain;
};

#endif
