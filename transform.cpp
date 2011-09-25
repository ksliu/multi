#include "transform.h"
#include <vector>
#include <iostream>
#include <cstdlib>
#include "geometry.h"
using std::vector;

transform::transform()
{
    subject = query = subjectAliged = queryAligned = queryShadow = 0;
    subjectLength = queryLength = sqLength = 0;
}
transform::~transform()
{
    delete [] subjectAliged;
    delete [] queryAligned;
    delete [] queryShadow;
}
void transform::set(double a[][3], int la,  double b[][3], int lb)
{
    subject = a;
    subjectLength = la;

    query = b;
    lq = lb;

    queryShadow = new double [lb][3];
}
void transform::run(const vector<int> &corres)
{
    int sz = alignCorres.size();
    if ( sz != queryLength )
    {
        cout << "error " << endl;
        exit(1);
    }
    numAliged = 0;
    for (int i=0; i<sz; i++)
        if (alignCorres[i] != -1)
            numAliged++;


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

}
void transform::getPartial(int n)
{
    delete [] subjectAliged;
    delete [] queryAligned;
    subjectAliged = double [n][3];
    queryAligned = double [n][3];

}

#include "dar.h"
typedef double PA [][3];
void foo(PA subject, int subjectLength,  PA query, int queryLength, const vector<int> &corres)
{
}
void foo()
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
