#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include "geometry.h"
#include "protein.h"
#include "pairString.h"
#include "multiString.h"

using namespace std;

class PairAlign
{
public:
    PairAlign(protein *s, protein *q);
    void update(const std::vector<int> &d);

private:
    std::vector<int> corres;
    SFPgr similarFragments;
    Dar subject, query, shadow;

    void updateCoordinatesByCorres(); // math: matrix rotation
    void updateCorresByCoordinates(); // CLE string: simliar pair fragment list

    int getAlignedLength() const;
    void refreshCorres(std::vector<int> & exclusionMarker, double topFraction, double absErr);
    bool absDev(int k1, int k2, double err) const;
    //    double rotation[3][3], translation[3];
};

PairAlign::PairAlign (protein *s, protein *q)
    :  corres(s->ca.len(), -1), similarFragments(s->cl, q->cl)
{
    subject.reset(s->ca.len());
    query.reset(q->ca.len());
    shadow.reset(q->ca.len());

    int i, j;
    for (i=0; i<s->ca.len(); i++)
        for (j=0; j<3; j++)
            subject[i][j] = s->ca[i][j];

    for (i=0; i<q->ca.len(); i++)
        for (j=0; j<3; j++)
            query[i][j] = q->ca[i][j];
}

void PairAlign::update(const vector<int> & d)
{
    corres = d;
    updateCoordinatesByCorres();
    updateCorresByCoordinates();
    //    return corres;
}

void PairAlign::updateCoordinatesByCorres()
{
    int nAligned = getAlignedLength();
    Dar subjectAliged(nAligned), queryAligned(nAligned);
    double subjectAligedCenter[3]={0}, queryAlignedCenter[3]={0};

    int i, j, k, n;
    for (n=i=0; i<subject.len(); i++)
    {
        if (corres[i] != -1)
        {
            for (j=0; j<3; j++)
            {
                subjectAliged[n][j] = subject[i][j];
                queryAligned[n][j] = query[i][j];

                subjectAligedCenter[j] += subject[i][j];
                queryAlignedCenter[j] += query[i][j];
            }
            n++;
        }
    }

    for (j=0; j<3; j++)
    {
        subjectAligedCenter[j] /= nAligned;
        queryAlignedCenter[j] /= nAligned;
    }
    for (i=0; i<nAligned; i++)
    {
        for ( j=0; j<3; j++)
        {
            subjectAliged[i][j] -= subjectAligedCenter[j];
            queryAligned[i][j] -= queryAlignedCenter[j];
        }
    }

    double rotation[3][3];
    geom::kabsch(subjectAliged.ptr(), queryAligned.ptr(), nAligned, rotation);

    //shift
    for ( j = 0; j < 3; j++)
    {
        for (i=0; i<subject.len(); i++)
            subject[i][j] -= subjectAligedCenter[j];
        for (i=0; i<query.len(); i++)
            query[i][j] -= queryAlignedCenter[j];
    }
    // rotate
    for (k = 0; k < query.len(); k++)
        for (i = 0; i < 3; i++)
            for (shadow[k][i] = 0., j = 0; j < 3; j++)
                shadow[k][i] += rotation[i][j] * query[k][j];
}
int PairAlign::getAlignedLength() const
{
    int nAligned = 0;
    for (vector<int>::size_type ix=0; ix!=corres.size(); ++ix)
    {
        if (corres[ix] != -1)
            nAligned++;
    }
    return nAligned;
}
void PairAlign::updateCorresByCoordinates()
{
    vector<int> scon(similarFragments.sq.size(), 0);
    refreshCorres(scon, 0.8, 7.5);
}

void PairAlign::refreshCorres(vector<int> &exclusionMarker, double topFraction, double absErr)
{
    corres.assign(subject.len(), -1);
    int part = static_cast<int> (topFraction * similarFragments.sq.size());
    int i, k, k1, k2, wrongPoints;
    for (i = 0; i < part; i++)
    {
        if (exclusionMarker[i] == 1) continue;
        k1 = similarFragments.sq[i].ia;
        k2 = similarFragments.sq[i].ib;
        for (wrongPoints = k = 0; k < SFPgr::SW; k++)
        {
            if (absDev(k1 + k, k2 + k, absErr))
            {
                if (corres[k1 + k] == -1)
                    corres[k1 + k] = k2 + k;
            }
            else
            {
                wrongPoints++;
            }
            if (wrongPoints > SFPgr::SW / 2)
            {
                exclusionMarker[i] = 1;
                break;
            }
        }
    }
}
bool PairAlign::absDev(int i, int j, double err) const
{
    return (fabs(subject[i][0] - shadow[j][0]) < err
            && fabs(subject[i][1]- shadow[j][1]) < err
            && fabs(subject[i][2] - shadow[j][2]) < err );
}
class MultiAlign
{
public:
    MultiAlign(protein *p, int numProtein);

    void moveToAverageSubject(protein *);

    void findMissingMotif();

    void outputAlignResult() const;
private:
//    HSFBgr similarBlock;

    Dar averageSubject;
    vector < vector<int> > mab;
    void alignPair(protein *subject, protein *query, vector<int> & coores);
};


MultiAlign::MultiAlign(protein *p, int numProtein)
{
    vector<string> pcl;
    for (int i=0; i<numProtein; i++)
        pcl.push_back(p[i].cl);

//    HSFBgr hs(pcl);

}

//void foo()
//{
//    protein *p;
//    int numProtein;

//    vector<string> pcl;
//    for (int i=0; i<numProtein; i++)
//    {
//        pcl.push_back(p[i].cl);
//    }

//    HSFBgr h(pcl);

//    // recenter the Highest score block
//    h.reCenterByHighestScore();

//    int center = h.center();

//    for (int loop = 0; loop < 3; loop++)
//    {
//        for (int i=0; i<numProtein; i++)
//        {
//            if ( i != center )
//            {
//                alignPair(p[i], p[center]);
//            }
//        }
//        alignToAverage(p, numProtein);
//    }
//    // find missing motif

////    cellRegister(p, numProtein);

//    // output
//    // 1. aligned block, coorespondence, rmsd, transformation matrix
//    // 2. coordinates file after aligned
//}
