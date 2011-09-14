#include <iostream>
#include <cmath>
#include <algorithm>
#include <sstream>
#include "corres.h"
#include "geometry.h"

using namespace std;
using geom::kabsch;

const int corres::N;

corres::corres() :  sublen(0), quelen(0), ali(0)
{
}

void corres::set(int a[], double sub[][3], int sl, double que[][3], int ql)
{
    int i, j;
    for (i = 0; i < sl; i++)
        for (j = 0; j < 3; j++)
            subject[i][j] = sub[i][j];
    for (i = 0; i < ql; i++)
        for (j = 0; j < 3; j++)
            query[i][j] = que[i][j];
    sublen = sl;
    quelen = ql;
    ali = a;
}
void corres::get_rotate_matrix(double rot[][3]) const
{
    int i, j;
    for (i = 0; i < 3; i++)
        for (j = 0; j < 3; j++)
            rot[i][j] = rotation[i][j];
}
void corres::full_strr(const string &rs, const string &rq, string &rs_out,
                       string &rq_out) const
{
    rs_out.clear();
    rq_out.clear();

    int i, j, k, pi = 0, pj = 0;
    for (i = 0; i < sublen; i++)
    {
        j = ali[i];
        if (j == -1)
            continue;
        // The alignment strings are extended when two match points met.
        for (k = pi; k < i; k++)
        {
            rs_out.push_back(rs[k]);
            rq_out.push_back('-');
        }
        for (k = pj; k < j; k++)
        {
            rs_out.push_back('-');
            rq_out.push_back(rq[k]);
        }
        rs_out.push_back(rs[i]);
        rq_out.push_back(rq[j]);

        pi = i + 1;
        pj = j + 1;
    }
    // The tail parts are counted when the above loop is terminated.
    for (k = pi; k < sublen; k++)
    {
        rs_out.push_back(rs[k]);
        rq_out.push_back('-');
    }
    for (k = pj; k < quelen; k++)
    {
        rs_out.push_back('-');
        rq_out.push_back(rq[k]);
    }
}
void corres::fullAlignString(const string &rs, const string &rq,
                             double cutoff1, double cutoff2, string &rs_out, string &rq_out,
                             string & mid) const
{
    rs_out.clear();
    rq_out.clear();
    mid.clear();

    double sqc[2], dist;
    sqc[0] = cutoff1 * cutoff1;
    sqc[1] = cutoff2 * cutoff2;
    if (sqc[0] > sqc[1])
        swap(sqc[0], sqc[1]);

    int i, j, k, pi = 0, pj = 0;
    for (i = 0; i < sublen; i++)
    {
        j = ali[i];
        if (j == -1)
            continue;
        // The alignment strings are extended when two match points met.
        for (k = pi; k < i; k++)
        {
            rs_out.push_back(rs[k]);
            rq_out.push_back('-');
            mid.push_back(' ');
        }
        for (k = pj; k < j; k++)
        {
            rs_out.push_back('-');
            rq_out.push_back(rq[k]);
            mid.push_back(' ');
        }
        rs_out.push_back(rs[i]);
        rq_out.push_back(rq[j]);

        dist = sq_dist(i, j);

        if (dist < sqc[0])
            mid.push_back(':');
        else if (dist < sqc[1])
            mid.push_back('.');
        else
            mid.push_back(' ');
        pi = i + 1;
        pj = j + 1;
    }
    // The tail parts are counted when the above loop is terminated.
    for (k = pi; k < sublen; k++)
    {
        rs_out.push_back(rs[k]);
        rq_out.push_back('-');
    }
    for (k = pj; k < quelen; k++)
    {
        rs_out.push_back('-');
        rq_out.push_back(rq[k]);
    }

}

void corres::concise_strr(const string &rs, const string &rq, string &rs_out,
                          string &rq_out, int threshold) const
{
    rs_out.clear();
    rq_out.clear();

    int i, j, pi = 0, pj = 0;
    string ta, tb;

    for (i = 0; i < sublen; i++)
        if ((j = ali[i]) != -1)
            break;
    pi = i;
    pj = j;

    for (i = pi + 1; i < sublen; i++)
    {
        if ((j = ali[i]) == -1)
            continue;
        if (i == pi + 1 && j == pj + 1)
        {
            if (ta.empty())
            {
                ta.push_back(rs[pi]);
                tb.push_back(rq[pj]);
            }
            ta.push_back(rs[i]);
            tb.push_back(rq[j]);
        }
        else // when mismatch
        {
            if (static_cast<int> (ta.size()) >= threshold)
            {
                if (!rs_out.empty())
                {
                    rs_out.push_back('-');
                    rq_out.push_back('-');
                }
                rs_out += ta;
                rq_out += tb;
            }
            ta.clear();
            tb.clear();
        }
        pi = i;
        pj = j;
    }
    // tail part
    if (static_cast<int> (ta.size()) >= threshold)
    {
        if (!rs_out.empty())
        {
            rs_out.push_back('-');
            rq_out.push_back('-');
        }
        rs_out += ta;
        rq_out += tb;
    }
}

void corres::minAlignPair(const std::string &rs, const std::string &rq,
                          std::string &rs_out, std::string &rq_out, std::string & summary,
                          int threshold) const
{
    rs_out.clear();
    rq_out.clear();

    stringstream ss;

    int i, j, pi = 0, pj = 0, si = 0, sj = 0;
    string ta, tb;

    for (i = 0; i < sublen; i++)
        if ((j = ali[i]) != -1)
            break;
    pi = i;
    pj = j;

    for (i = pi + 1; i < sublen; i++)
    {
        if ((j = ali[i]) == -1)
            continue;
        if (i == pi + 1 && j == pj + 1)
        {
            if (ta.empty())
            {
                ta.push_back(rs[pi]);
                tb.push_back(rq[pj]);
                si = pi;
                sj = pj;
            }
            ta.push_back(rs[i]);
            tb.push_back(rq[j]);
        }
        else // when mismatch
        {
            if (static_cast<int> (ta.size()) >= threshold)
            {
                if (!rs_out.empty())
                {
                    rs_out.push_back('-');
                    rq_out.push_back('-');
                }
                rs_out += ta;
                rq_out += tb;
                ss << si << "-" << si + ta.size() -1<< ","
                   << sj << "-" << sj + tb.size() -1<< " "
                   << ta << " " << tb << endl;
            }
            ta.clear();
            tb.clear();
        }
        pi = i;
        pj = j;
    }
    // tail part
    if (static_cast<int> (ta.size()) >= threshold)
    {
        if (!rs_out.empty())
        {
            rs_out.push_back('-');
            rq_out.push_back('-');
        }
        rs_out += ta;
        rq_out += tb;
        ss << si << "-" << si + ta.size() -1 << ","
           << sj << "-" << sj + tb.size() -1 << " "
           << ta << " " << tb << endl;
    }
    summary = ss.str();

}
int corres::extract(double px[][3], double xc[3], double py[][3], double yc[3]) const
{
    int i, j, k, s;
    for (s = i = 0; i < sublen; i++)
    {
        if ((k = ali[i]) == -1)
            continue;
        for (j = 0; j < 3; j++)
        {
            px[s][j] = subject[i][j];
            py[s][j] = query[k][j];
        }
        s++;
    }
    for (j = 0; j < 3; j++)
    {
        xc[j] = 0;
        yc[j] = 0;
        for (i = 0; i < s; i++)
        {
            xc[j] += px[i][j];
            yc[j] += py[i][j];
        }
        xc[j] /= s;
        yc[j] /= s;
    }
    for (i = 0; i < s; i++)
        for (j = 0; j < 3; j++)
        {
            px[i][j] -= xc[j];
            py[i][j] -= yc[j];
        }
    return s; // number of points aligned
}

void corres::ontop()
{
    static double px[N][3], py[N][3];
    static double xc[3], yc[3]; // centroid of the aligned points

    int i, j, k, ali_len;
    ali_len = extract(px, xc, py, yc);
    kabsch(px, py, ali_len, rotation);
    //shift
    for (i = 0; i < sublen; i++)
        for (j = 0; j < 3; j++)
            subject[i][j] -= xc[j];

    for (i = 0; i < quelen; i++)
        for (j = 0; j < 3; j++)
            query[i][j] -= yc[j];
    //rotate
    for (k = 0; k < quelen; k++)
        for (i = 0; i < 3; i++)
            for (shadow[k][i] = 0., j = 0; j < 3; j++)
                shadow[k][i] += rotation[i][j] * query[k][j];
}

bool corres::increased() const
{
    bool ret = true;
    int i, j;
    for (i = 0, j = 1; j < sublen; j++)
    {
        if (ali[j] != -1)
        {
            if (ali[j] < ali[i])
            {
                cerr << "Warning: not sequential at: " << j << " align to  "
                     << ali[j] << " and " << i << " align to  " << ali[i]
                     << endl;

                if (ret)
                    ret = false;
            }
            i = j;
        }
    }
    return ret;
}

bool corres::colinear(int k1, int k2) const
{
    int i;
    for (i = k1 + 1; i < sublen; i++)
        if (ali[i] != -1 && k2 > ali[i])
            return false;
    for (i = k1 - 1; i >= 0; i--)
        if (ali[i] != -1)
            return k2 > ali[i];
    return true;
}

double corres::rmsd() const
{
    int i, j, k;
    double rd;
    for (rd = 0., k = i = 0; i < sublen; i++)
    {
        j = ali[i];
        if (j != -1)
        {
            k++;
            rd += sq_dist(i, j);
        }
    }
    return sqrt(rd / k);
}

double corres::TMscore() const
{
    if (quelen < 15)
    {
        cerr << "target protein too short" << endl;
        return 0;
    }

    int i, j;
    double t, rd, d = 1.24 * pow(quelen - 15, 1. / 3) - 1.8;
    for (rd = 0., i = 0; i < sublen; i++)
    {
        j = ali[i];
        if (j != -1)
        {
            t = sqrt(sq_dist(i, j)) / d;
            rd += 1 / (1 + t * t);
        }
    }
    return rd / quelen;
}

bool corres::abs_dev(int k1, int k2, double err) const
{
    return (fabs(subject[k1][0] - shadow[k2][0]) < err
            && fabs(subject[k1][1]- shadow[k2][1]) < err
            && fabs(subject[k1][2] - shadow[k2][2]) < err );
}

bool corres::sq_dev(int k1, int k2, double err) const
{
    int i;
    double t, s;
    for (s = 0.0, i = 0; i < 3; i++)
    {
        t = subject[k1][i] - shadow[k2][i];
        s += t * t;
    }
    return s < err; // Square root is not extracted.
}

double corres::sq_dist(int k1, int k2) const
{
    int i;
    double t, s;
    for (s = 0.0, i = 0; i < 3; i++)
    {
        t = subject[k1][i] - shadow[k2][i];
        s += t * t;
    }
    return s;
}
int corres::aligned_size() const
{
    int i, s;
    for (s = i = 0; i < sublen; i++)
        if (ali[i] != -1)
            s++;
    return s;
}
