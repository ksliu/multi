#ifndef PAIRALIGN_H
#define PAIRALIGN_H
#include <iostream>
#include <vector>
#include <string>
#include "dar.h"
#include "protein.h"
#include "pairString.h"

class PairAlign
{
public:
    PairAlign();
    PairAlign(protein *s, protein *q);
    void set(protein *s, protein *q);
    void update(std::vector<int> &d);

private:
    std::vector<int> corres;
    SFPgr similarFragments;
    Dar *pSubjectP, *pQuery; // shadow;

    void updateCoordinatesByCorres(); // math: matrix rotation
    void updateCorresByCoordinates(); // CLE string: simliar pair fragment list

    int getAlignedLength() const;
    void refreshCorres(std::vector<int> & exclusionMarker, double topFraction, double absErr);
    bool absDev(int k1, int k2, double err) const;
    //    double rotation[3][3], translation[3];
};



#endif // PAIRALIGN_H
