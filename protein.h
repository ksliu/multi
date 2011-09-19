#ifndef  PROTEIN_H
#define  PROTEIN_H
#include <string>
#include "dar.h"

class protein
{
public:
    protein ();
    protein (const std::string &);
    void loadFile (const std::string &);
    void toXYZFile(const std::string &);

    std::string id, aa, cl;
    Dar ca;

private:
    protein(const protein &);
    protein operator=(const protein &);
};

#endif  /*PROTEIN_H*/
