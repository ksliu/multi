#include "protein.h"
#include "pdb.h"

const int protein::N;

protein::protein(const char *fn)
{
    id = fn;
    PDB pdb(fn);
    pdb.getAA(aa);
    pdb.getCa(ca);
    cl = bio::ca2cl(ca.ptr(), ca.len());
}
