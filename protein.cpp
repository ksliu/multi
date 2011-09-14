#include "protein.h"
#include "bioinfo.h"

const int protein::N;

protein::protein(const char *fn)
{
    id = fn;
    len = bio::readPDB(fn, aa, ca, N);
    cl = bio::ca2cl(ca, len);
}
