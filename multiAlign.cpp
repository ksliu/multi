#include <iostream>
#include <vector>
#include <string>
#include "dar.h"
#include "protein.h"
#include "multiString.h"
#include "multiAlign.h"
#include "pairAlign.h"
using namespace std;



MultiAlign::MultiAlign(protein *p, int n) :pro(p), numProtein(n)
{
    vector<string> pcl;
    for (int i=0; i<numProtein; i++)
        pcl.push_back(pro[i].cl);
    hsfbgr.set(pcl);

    multiCorres.resize(numProtein);
    mpa = new PairAlign [numProtein];
}

MultiAlign::~MultiAlign()
{
    delete []mpa;
}


void MultiAlign::moveUseHighestBlock()
{
    const HSFB & hsfb = hsfbgr.getFirstBlock();
    subjectIndex = hsfbgr.getBlockCenter(hsfb);
    cout << "center: "<<  subjectIndex << endl;

    for (int i=0; i<numProtein; i++)
    {
        multiCorres[i].assign(pro[subjectIndex].ca.len(), -1);

        for (int j=0; j<HSFBgr::SFB_WIDTH; j++)
            multiCorres[i][hsfb.positions[subjectIndex] + j] = hsfb.positions[i] + j;

        if (i != subjectIndex)
        {
            mpa[i].set(pro + subjectIndex, pro + i);
            mpa[i].update(multiCorres[i]);
        }
    }
}

void MultiAlign::updateAverageSubject()
{
    int n = pro[subjectIndex].ca.len();
    averageSubject.reset(n);
    int i, j, k;
    for (i=0; i<n; i++)
        for (j=0; j<3; j++)
            averageSubject[i][j] = pro[subjectIndex][i][j];

    std::vector<int> abDepth(n, 0);
    for (k=0; k<numProtein; k++)
    {
        if (k != subjectIndex) continue;
        for (i=0; i<n; i++)
        {
            int p = multiCorres[k][i];
            if (p != -1)
            {
                for (j=0; j<3; j++)
                    averageSubject[i][j] += pro[k][p][j];
                abDepth[i]++;
            }
        }
    }
    for (i=0; i<n; i++)
    {
        if (abDepth[i] > 0)
        {
            for (j=0; j<3; j++)
                averageSubject[i][j] /= abDepth[i];
        }
    }
}

void MultiAlign::allMoveToAverageSubject()
{
    int i, j;

    for (i=0; i<numProtein; i++)
    {
        ddtransform(averageSubject, pro[i].ca, multiCorres[i]);
    }
}

void MultiAlign::ddtransform(const Dar &subject, Dar &query, vector<int> &corres)
{
    int sublen = corres.size();
    if ( subject.len() != sublen )
    {
        cout << "length mismatch " << endl;
        exit(1);
    }
    int numAligned = 0;
    for (int i=0; i < sublen; i++)
    {
        if (corres[i] != -1)
            numAligned++;
    }

    Dar subjectAliged(numAligned), queryAligned(numAligned);
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
        subjectAligedCenter[j] /= numAligned;
        queryAlignedCenter[j] /= numAligned;
    }
    for (i=0; i<numAligned; i++)
    {
        for ( j=0; j<3; j++)
        {
            subjectAliged[i][j] -= subjectAligedCenter[j];
            queryAligned[i][j] -= queryAlignedCenter[j];
        }
    }

    double rotation[3][3];
    geom::kabsch(subjectAliged.ptr(), queryAligned.ptr(), numAligned, rotation);

    //  A - A_c ~ R*( B - B_c )   =>   A ~  R*B + A_c - R*B_c
    Dar shadow(query.len());
    for (k = 0; k < query.len(); k++)
    {
        for (i = 0; i < 3; i++)
        {
            shadow[k][i] = 0;
            for (j = 0; j < 3; j++)
                shadow[k][i] += rotation[i][j] * query[k][j];
        }
    }
    double shift[3];
    for (i = 0; i < 3; i++)
    {
        shift[i] = 0;
        for (j = 0; j < 3; j++)
            shift[i] += rotation[i][j] * queryAlignedCenter[j];
        shift[i] = subjectAligedCenter[i] - shift[i];
    }
    for (k = 0; k < query.len(); k++)
        for (j = 0; j < 3; j++)
            query[k][j] = shadow[k][j] + shift[j];
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
