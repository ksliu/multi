#ifndef MULTIALIGN_H
#define MULTIALIGN_H
#include <iostream>
#include <vector>
#include <string>
#include "dar.h"
#include "protein.h"
#include "multiString.h"
#include "pairString.h"

class MultiAlign
{
public:
    MultiAlign(protein *, int);

    void run();
    void outputAlignResult(const std::string & fn) const;
    void genScript(const std::string &fn) const;

private:
    std::vector<protein *> pro;
    int np, subjectIndex;

    std::vector < std::vector<int> > multiCorres;
    std::vector <  SFPGenerator >  multiSFP;

    Dar averageSubject, centers;

    void moveToCenter();
    void updateAverageSubject();

    void transToSubject(int index);
    void transToAverage(int index);

    void fillCorresByCLE(int index);
    void tuneCorresByCLE(int index);

    void moveToAverageSubject();
    void refreshSubjectCorres();

    void findMissingMotif();

    void getRMSD(int minDepth, int &alignedSize, double &rmsd) const;


    static void ddtransform(const Dar &subject, Dar &query, const std::vector<int> &corres);
    static bool absDev(const double p[], const double q[], double err);
    static bool disDev(const double p[], const double q[], double err);
    static bool colinear(const std::vector<int> &corres, int s, int q);
};

#endif // MULTIALIGN_H
