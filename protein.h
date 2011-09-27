#ifndef  PROTEIN_H
#define  PROTEIN_H
#include <string>
#include <iostream>
#include "dar.h"

class protein
{
public:
    protein ();
    protein (const std::string &);
    void loadFile (const std::string &);
    void toXYZFile(const std::string &);

    void genFakePDB(char chain, int startAtomNo, std::ostream &);

    std::string id, aa, cl;
    Dar ca;

private:
    protein(const protein &);
    protein operator=(const protein &);
};

#endif  /*PROTEIN_H*/
