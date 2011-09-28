#ifndef  PROTEIN_H
#define  PROTEIN_H
#include <string>
#include <iostream>
#include "dar.h"

class protein
{
public:
    protein ();
    protein (const std::string &, char);
    void loadFile (const std::string &, char);
    void toXYZFile(const std::string &);
    void genFakePDB(char chain, int startAtomNo, std::ostream &);

    std::string aa, cl;
    Dar ca;

private:
    std::string id, chain;

    protein(const protein &);
    protein operator=(const protein &);


};

#endif  /*PROTEIN_H*/
