#ifndef  PROTEIN_H
#define  PROTEIN_H
#include <string>

class protein
{
private:
    static const int N = 3000;
    protein(const protein &);
    protein operator=(const protein &);

public:
    protein(const char *fn);

    std::string id, aa, cl;
    double ca[N][3];
    int len;
};

#endif  /*PROTEIN_H*/
