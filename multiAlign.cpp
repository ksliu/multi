#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <string>
#include <cstdlib>
#include <cmath>
#include "geometry.h"
#include "dar.h"
#include "protein.h"
#include "multiString.h"
#include "multiAlign.h"
using namespace std;

MultiAlign::MultiAlign(protein *allPro, int numAllPro)
{
    vector<string> pcl;

    for (int i=0; i<numAllPro; i++)
        pcl.push_back(allPro[i].cl);

    HSFBGenerator similarBlock(pcl);

    const HSFB & hsfb = similarBlock.getFirstBlock();

    if (hsfb.depth < numAllPro)
    {
        cout << "Warning: fist block is not full"  << endl;
    }
    else if (hsfb.depth > numAllPro)
    {
        cout << "error " << endl;
        exit(1);
    }

    np = hsfb.depth;
    multiCorres.resize(np);
    multiSFP.resize(np);

    int subjectIndexInAll = similarBlock.getBlockCenter(hsfb);

    for (int i=0, ni=0; i<numAllPro; i++)
    {
        if (hsfb.position[i] != -1)
        {
            pro.push_back(allPro + i);

            multiCorres[ni].assign(allPro[subjectIndexInAll].ca.len(), -1);
            for (int j=0; j<similarBlock.SFB_WIDTH; j++)
                multiCorres[ni][hsfb.position[subjectIndexInAll] + j] = hsfb.position[i] + j;

            if (i == subjectIndexInAll)
                subjectIndex = ni;
            else
                multiSFP[ni].set(allPro[subjectIndexInAll].cl, allPro[i].cl);

            ni++;
        }
        else
        {
            cout << "throw protein: " << i << endl;
        }
    }
    averageSubject.reset(allPro[subjectIndexInAll].ca.len());
    moveToCenter();
}

void MultiAlign::moveToCenter()
{
    centers.reset(np);
    int i, j, k;
    for (k=0; k<np; k++)
    {
        for (j=0; j<3; j++) centers[k][j] = 0;

        for (i=0; i<pro[k]->ca.len(); i++)
            for (j=0; j<3; j++)
                centers[k][j] += pro[k]->ca[i][j];

        for (j=0; j<3; j++) centers[k][j] /= pro[k]->ca.len();

        for (i=0; i<pro[k]->ca.len(); i++)
            for (j=0; j<3; j++)
                pro[k]->ca[i][j] -= centers[k][j];
    }
}


void MultiAlign::run()
{
    int i;
    for (i=0; i<np; i++)
    {
        if (i != subjectIndex)
        {
            transToSubject(i);
            fillCorresByCLE(i);
        }
    }
    refreshSubjectCorres();

    for (int loop = 0; loop <0; loop++)
    {
        updateAverageSubject();
        for (i=0; i<np; i++)
        {
            transToAverage(i);
        }
        for (i=0; i<np; i++)
        {
            if(i != subjectIndex)
                fillCorresByCLE(i);
        }
        refreshSubjectCorres();
    }

    updateAverageSubject();
    for (i=0; i<np; i++)
    {
        transToAverage(i);
    }
    for (i=0; i<np; i++)
    {
        if(i != subjectIndex)
            tuneCorresByCLE(i);
    }
    refreshSubjectCorres();

    findMissingMotif();

}
void MultiAlign::transToSubject(int index)
{
    ddtransform(pro[subjectIndex]->ca, pro[index]->ca, multiCorres[index]);
}
void MultiAlign::transToAverage(int index)
{
    ddtransform(averageSubject, pro[index]->ca, multiCorres[index]);
}

void MultiAlign::fillCorresByCLE(int index)
{
    vector<int> & corres = multiCorres[index];
    const vector<SFP> & sq = multiSFP[index].getCandidate();
    corres.assign(corres.size(), -1);

    int part = sq.size();
    const  double absErr = 7.5;

    for (int i = 0; i < part; i++)
    {
        int s = sq[i].ia;
        int q = sq[i].ib;
        int wrongPoints = 0;
        for (int k = 0; k < SFPGenerator::SW; k++, s++, q++)
        {
            if (corres[s] == -1
                    && absDev(pro[subjectIndex]->ca[s], pro[index]->ca[q], absErr)
                    && colinear(corres, s, q))
            {
                corres[s] = q;
            }
            else
            {
                wrongPoints++;
                if (wrongPoints > SFPGenerator::SW / 2)
                    break;
            }
        }
    }
}


void MultiAlign::tuneCorresByCLE(int index)
{
    vector<int> & corres = multiCorres[index];
    const vector<SFP> & sq = multiSFP[index].getCandidate();
    corres.assign(corres.size(), -1);

    int part = sq.size();
    const  double disErr = 5.0 * 5.0;

    for (int i = 0; i < part; i++)
    {
        int s = sq[i].ia ;
        int q = sq[i].ib ;
        for (int k = 0; k < SFPGenerator::SW; k++, s++, q++)
        {
            if ( corres[s] == -1
                    && disDev(pro[subjectIndex]->ca[s], pro[index]->ca[q], disErr)
                    && colinear(corres, s, q) )
            {
                corres[s] = q;
            }
        }
    }

    int i, j;

    for (i = 1; i <  pro[subjectIndex]->ca.len(); i++)
    {
        if (corres[i - 1] != -1 && corres[i] == -1)
        {
            j = corres[i - 1] + 1;
            if (disDev(pro[subjectIndex]->ca[i], pro[index]->ca[j], disErr))
                corres[i] = j;
        }
    }

    for (i =  pro[subjectIndex]->ca.len() - 1; i > 0; i--)
    {
        if (corres[i - 1] == -1 && corres[i] != -1)
        {
            j = corres[i];
            if (disDev(pro[subjectIndex]->ca[i - 1], pro[index]->ca[j - 1], disErr))
                corres[i - 1] = j - 1;
        }
    }
}

void MultiAlign::refreshSubjectCorres()
{
    for (int i=0; i<pro[subjectIndex]->ca.len(); i++)
    {
        int depth = 0;
        for (int k=0; k<np; k++)
        {
            if (k != subjectIndex && multiCorres[k][i] != -1)
                depth++;
        }
        if (depth > 0)
            multiCorres[subjectIndex][i] = i;
        else
            multiCorres[subjectIndex][i] = -1;
    }
}

void MultiAlign::updateAverageSubject()
{
    if (averageSubject.len() != pro[subjectIndex]->ca.len())
    {
        cout << "length mismatch" << endl;
        exit(1);
    }

    int i, j, k;
    for (i=0; i<averageSubject.len(); i++)
    {
        for (j=0; j<3; j++)
            averageSubject[i][j] = 0;

        int depth = 0;
        for (k=0; k<np; k++)
        {
            int p = multiCorres[k][i];
            if (p != -1)
            {
                for (j=0; j<3; j++)
                    averageSubject[i][j] += pro[k]->ca[p][j];
                depth++;
            }
        }
        if (depth > 0)
        {
            for (j=0; j<3; j++)
                averageSubject[i][j] /= depth;
        }
    }
}

bool MultiAlign::absDev(const double p[], const double q[], double err)
{
    return (fabs(p[0] - q[0]) < err
            && fabs(p[1]- q[1]) < err
            && fabs(p[2] - q[2]) < err );
}

bool MultiAlign::disDev(const double p[], const double q[], double err)
{
    double s = 0;
    for (int i=0; i<3; i++)
        s += (p[i] - q[i]) * (p[i] - q[i]);
    return  s < err;
}

bool MultiAlign::colinear(const vector<int> &corres, int s, int q)
{
    int i, sublen = corres.size();
    for (i=0; i<s; i++)
    {
        int p = corres[i];
        if (p != -1)
        {
            if ( p >= q )
                return false;
        }
    }
    for (i=s+1; i<sublen; i++)
    {
        int p = corres[i];
        if (p != -1)
        {
            if ( p <= q )
                return false;
        }
    }
    return true;
}


void MultiAlign::ddtransform(const Dar &subject, Dar &query, const vector<int> &corres)
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
        k = corres[i];
        if (k!= -1)
        {
            for (j=0; j<3; j++)
            {
                subjectAliged[n][j] = subject[i][j];
                queryAligned[n][j] = query[k][j];

                subjectAligedCenter[j] += subject[i][j];
                queryAlignedCenter[j] += query[k][j];
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
    double shift[3];
    for (i = 0; i < 3; i++)
    {
        shift[i] = 0;
        for (j = 0; j < 3; j++)
            shift[i] += rotation[i][j] * queryAlignedCenter[j];
        shift[i] = subjectAligedCenter[i] - shift[i];
    }

    double rq [3];
    for (k = 0; k < query.len(); k++)
    {
        for (i = 0; i < 3; i++)
        {
            rq[i] = 0;
            for (j = 0; j < 3; j++)
                rq[i] += rotation[i][j] * query[k][j];
        }
        for (i = 0; i < 3; i++)
            query[k][i] = rq[i] + shift[i];
    }
}

void MultiAlign::outputAlignResult(const std::string & fn) const
{
    // output
    // 1. aligned block, coorespondence, rmsd, transformation matrix
    // 2. coordinates file after aligned
    ofstream fout(fn.c_str());
    if (fout.fail())
    {
        cout << "cannot open file " << fn << endl;
        exit(1);
    }
    fout.setf(ios::fixed);
    fout.precision(3);

    int alignedSize[2];
    double rmsd[2];

    getRMSD(np, alignedSize[0], rmsd[0]);
    getRMSD(int(np * 0.6), alignedSize[1], rmsd[1]);

    fout << "subject protein: " << subjectIndex << endl
         << "aligned size:" << setw(8) <<  alignedSize[0] << setw(8) << alignedSize[1] << endl
         << "rmsd:        " << setw(8) <<rmsd[0] <<  setw(8) << rmsd[1] << endl;


    for (int j=0; j < pro[subjectIndex]->ca.len(); j++)
    {
        for (int i=0; i<np; i++)
            fout << setw(5) << multiCorres[i][j];
        fout << endl;
    }

    fout.close();
}
void MultiAlign::getRMSD(int minDepth, int &alignedSize, double &rmsd) const
{
    alignedSize = 0;
    rmsd = 0;
    if ( minDepth < 0 || minDepth > np)
    {
        cout << "minDepth ~ [0, number of protein] " << endl;
        exit(1);
    }
    int nPair = 0;
    for (int i=0; i< pro[subjectIndex]->ca.len(); i++)
    {
        int depth = 0;
        for (int k=0; k < np; k++)
        {
            if (multiCorres[k][i] != -1)
                depth++;
        }
        if (depth >= minDepth)
        {
            alignedSize++;
            nPair += depth;
            for (int k=0; k < np; k++)
            {
                int p = multiCorres[k][i];
                if ( p != -1)
                {
                    for (int j=0; j<3; j++)
                    {
                        double d = pro[k]->ca[p][j] - averageSubject[i][j];
                        rmsd += d * d;
                    }
                }
            }
        }
    }
    rmsd = sqrt(rmsd) / nPair;
}

void MultiAlign::findMissingMotif()
{
    // cellRegister
    //    class cell
    //    {
    //    public:
    //        cell():depth(0) { }
    //        int pos[3], depth;
    //    };
    //    class cellcmp
    //    {
    //        bool operator() (const cell &a, const cell &b)
    //        {
    //            return a.depth > b.depth;
    //        }
    //    };

    //    double min[3] = {0}, max[3] = {0};

    //    int i, j, k;

    //    for (k=0; k<np; k++)
    //    {
    //        for (i=0; i<pro[k]->ca.len(); i++)
    //        {
    //            for (j=0; j<3; j++)
    //            {
    //                double d = pro[k]->ca[i][j];
    //                if (d < min[j])  min[j] = d;
    //                if (d > max[j])  max[j] = d;
    //            }
    //        }
    //    }
    //    const double cubeSize = 6.0;

    //    int nSide[3];

    //    for (j=0; j<3; j++)
    //        nSide[j] =  int( (max[j] -  min[j]) / cubeSize ) + 1 ;

    //    int nTotal = 1;
    //    for (j=0; j<3; j++)
    //        nTotal *= nSide[j];


    //    // index = i * N_j + j + k * (N_i * N_j).
    //    vector< cell > reg;
    //    cell c;
    //    for (k=0; k<np; k++)
    //    {
    //        for (i=0; i<pro[k]->ca.len(); i++)
    //        {
    //            for (j=0; j<3; j++)
    //            {
    //                double d = pro[k]->ca[i][j];
    //                c.pos[j] = int ((d - min[j] ) / cubeSize);
    //            }
    //            ++[index];
    //        }
    //    }

}
void MultiAlign::genScript(const std::string &fn) const
{
    ofstream fout(fn.c_str());
    if (fout.fail())
    {
        cout << "cannot open file: " << fn << endl;
        exit(1);
    }
    fout << "load inline\n"
            "color white\n"
            "select all\n"
            "wireframe .45\n";

    fout.setf(ios::fixed);
    fout.precision(3);

    vector<int> atomStart;
    int allLen = 0;
    for (int i=0; i<np; i++)
    {
        atomStart.push_back(allLen);
        allLen += pro[i]->ca.len();
    }

    for (int i=0; i < pro[subjectIndex]->ca.len(); i++)
    {
        for (int k=0; k<np; k++)
        {
            int p = multiCorres[k][i];
            if ( p!= -1)
            {
                fout << "select atomno=" << setw(4) << atomStart[k] + p << endl;
                fout << "color red\n";
            }
        }
    }
    fout << "select all\nexit\n";

    for (int k=0; k<np; k++)
    {
        pro[k]->genFakePDB(' ', atomStart[k], fout);
    }
    fout.close();
}

