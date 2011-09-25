#ifndef MULTIALIGN_H
#define MULTIALIGN_H
#include <iostream>
#include <vector>
#include <string>
#include "dar.h"
#include "protein.h"
#include "multiString.h"
#include "pairAlign.h"

class MultiAlign
{
public:
    MultiAlign(protein *p, int n);

    void moveUseHighestBlock();

    void allMoveToAverageSubject();

    void findMissingMotif();

    void updateAverageSubject();

    void outputAlignResult() const;


    ~MultiAlign();
private:
    protein *pro;
    int numProtein, subjectIndex;

    Dar averageSubject;
    HSFBgr hsfbgr;
    std::vector < std::vector<int> > multiCorres;

    PairAlign *mpa;

    void ddtransform(const Dar &subject, Dar &query, std::vector<int> &corres);

};
#endif // MULTIALIGN_H
