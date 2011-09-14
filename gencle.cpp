/*
 * =============================================================================
 * Description: generate comformational letters form angles and coordinates
 * utilities:
 *	    1) translate 3 angles(btb) to 1 CL code
 *	    Source: Wei-Mou Zheng's paper,
 *		    Protein Conformational alphabets, Appendix A2
 *	    2) translate a series of CA atom coordinates in the main chain to CL codes
 * =============================================================================
 */
#include <iostream>
#include <cstdlib>
#include <cmath>
#include "gencle.h"

namespace
{
//pi: mode weight; sig: square root of |sigma|^-1; mid: means (mu)
//var: inverse convariance matrices.
int id[17]=
{'I','J','H','K','F','E','C','D','A','B','G','L','M','N','O','P','Q'};
double pi[17] =
{0.0819,0.0727,0.1617,0.0594,0.0491,0.1158,0.0749,0.0544,0.0433,0.0387,0.0556,0.0533,0.0371,0.0315,0.0214,0.0318,0.0173};
double sig[17] =
{1881.06,1796.72,10424.63,254.11,104.55,108.95,99.85,77.74,202.76,66.29,132.85,40.02,144.28,73.52,247.16,206.13,24.82};
double mid[17][3] =
{
    {1.517,  0.833, 1.517}, {1.577,  1.045, 1.550}, {1.550,  0.878, 1.552},
    {1.479,  0.699, 1.426}, {1.094, -2.722, 0.909}, {1.021, -2.977, 0.954},
    {1.014, -1.878, 1.138}, {0.786, -2.303, 1.028}, {1.025, -2.003, 1.551},
    {1.057, -2.940, 1.344}, {1.491,  2.087, 1.046}, {1.396,  0.754, 0.840},
    {1.466,  1.643, 1.440}, {1.121,  0.142, 1.490}, {1.537, -1.891, 1.484},
    {1.235, -2.982, 1.493}, {0.863, -0.372, 1.008}
};
double var[17][3][3] =
{
    {{275.430 , -28.322 , 106.938}, {-28.322 , 84.345  ,-46.140 }, {106.938 ,-46.140  ,214.445 }},
    {{314.298 , -10.252 , 37.793 }, {-10.252 , 46.021  , -69.977}, {37.793  ,-69.977  ,332.804 }},
    {{706.634 , -93.940 , 128.933}, {-93.940 , 245.505 ,-171.758}, {128.933 , -171.758,786.116 }},
    {{73.820  , -13.740 , 15.478 }, {-13.740 ,  21.548 , -25.281}, {15.478  ,-25.281  , 75.727 }},
    {{24.130  , 1.909   , -11.179}, {1.909   , 10.921  ,  -8.787}, {-11.179 , -8.787  ,  53.040}},
    {{34.340  , 4.171   , -9.288 }, {4.171   , 15.243  , -22.456}, {-9.288  , -22.456 ,  56.837}},
    {{28.008  , 4.114   , 2.318  }, {4.114   ,  6.183  , -5.128 }, {2.318   , -5.128  , 69.368 }},
    {{56.205  , 3.764   , -10.830}, {3.764   ,  4.168  , -2.148 }, {-10.830 ,  -2.148 , 30.057 }},
    {{30.524  , 9.083   ,  6.041 }, {9.083   , 8.692   , 5.674  }, {6.041   , 5.674   ,228.584 }},
    {{26.872  , 4.632   , 9.545  }, {4.632   , 4.865   ,  -4.959}, {9.545   , -4.959  ,  54.317}},
    {{163.926 , 0.615   , 1.975  }, {0.615   , 3.777   , -3.748 }, {1.975   ,-3.748   ,32.287  }},
    {{43.653  ,  2.493  ,  -6.966}, {2.493   ,  1.425  ,  -2.864}, {-6.966  , -2.864  , 34.454 }},
    {{72.892  , 2.130   , 1.870  }, {2.130   , 4.842   ,-7.864  }, {1.870   , -7.864  , 72.914 }},
    {{25.283  , 3.210   , 9.908  }, {3.210   , 3.113   , 0.872  }, {9.908   , 0.872   , 82.959 }},
    {{170.760 , -0.715  , -4.138 }, {-0.715  , 3.726   , 3.069  }, {-4.138  , 3.069   , 98.705 }},
    {{47.969  , 8.244   , -4.912 }, {8.244   , 7.340   , -6.570 }, {-4.912  , -6.570  ,155.581 }},
    {{28.408  ,  1.513  , 3.380  }, {1.513   , 1.214   , 0.059  }, {3.380   , 0.059   , 19.549 }},
};
/*
 * Gaussian distribution of # ic
 * @input: ic: for a specific mode
 * @input: btb: bending angle, torsion angle, bending angle
 * @return: probability of btb belong to ic mode
 */
double normal(int ic, double btb[3])
{
    int i, j;
    double y[3], z = 0;

    for (i=0; i<3; i++)
        y[i] = btb[i] - mid[ic][i];

    if (y[1] > 3.1416) y[1] -= 6.2832;
    else if (y[1] < -3.1416) y[1] += 6.2832;
    for (i = 0; i < 3; i++)
        for (j = 0; j <= i; j++)
            z += (i-j ? 1. : .5) * y[i] * var[ic][i][j] * y[j];
    return (exp(-z) * sig[ic]);
}
/*
 * @brief  find the mode which has maximum probability among all modes
 * @input: btb: bending angle, torsion angle, bending angle
 * @return: one CL code corresponding to the mode
 */
char btb2cl(double btb[3])
{
    int i, ip=0;
    double q, qq=-1.0;
    for(i=0; i<17; i++)
    {
        if((q=pi[i]*normal(i,btb)) > qq)
        {
            qq=q;
            ip=i;
        }
    }
    return id[ip];
}

// translate part
void normalise(double x[3])
{
    int i;
    double s;
    for (s=0., i=0; i<3; i++)
        s += x[i] *x[i];
    s = sqrt(s);
    for (i=0; i<3; i++)
        x[i] /= s;
}
double dot(const double x[3], const double y[3])
{
    int i;
    double s;
    for (s=0., i=0; i<3; i++)
        s += x[i] *y[i];
    return s;
}
void cross(const double x[3], const double y[3], double r[3])
{
    r[0] = x[1] * y[2] - x[2] *y[1];
    r[1] = x[2] * y[0] - x[0] *y[2];
    r[2] = x[0] * y[1] - x[1] *y[0];
}
double nbend(const double x[3], const double y[3]) // assuming these vectors are normalised
{
    int i;
    double s;
    for (s=0.0, i=0; i<3; i++)
        s += x[i] *y[i];
    return acos(s);
}

char normalise_range(double x[3])
{
    int i;
    double s;

    for (s=0., i=0; i<3; i++) s += x[i] * x[i];
    s = sqrt(s);

    for (i=0; i<3; i++) x[i] /= s;

    // CA distance \in (2.5, 4.5)
    return (s < 2.5 || s > 4.5) ? '1': '0';
}
}

namespace bio
{
std::string ca2cl(const double p[][3], int len)
{
    if (len < 4)
        return std::string(len, 'R');

    int i, j;
    // n: normal vector of the first and second plane
    // a: bending, torsion, bending angles
    // d: vectors of adjacent points
    double a[3], n[2][3], d[2][3];

    // cl: return value
    // brk: rectify break points
    std::string cl("RR"), brk("0");

    for (i=0; i<2; i++)
    {
        for (j=0; j<3; j++)
            d[i][j] = p[i+1][j] - p[i][j]; // ~ 3.8A
        brk.push_back(normalise_range(d[i]));
    }

    // a[2], n[1] were set before the first run
    a[2] = nbend(d[0], d[1]);
    cross(d[0], d[1], n[1]);
    normalise(n[1]);

    for (i=3; i<len; i++)
    {
        //shift
        a[0] = a[2];
        for (j=0; j<3; j++)
        {
            n[0][j] = n[1][j];
            d[0][j] = d[1][j];
        }

        for (j=0; j<3; j++)
            d[1][j] = p[i][j] - p[i-1][j]; // ~ 3.8A
        brk.push_back(normalise_range(d[1]));

        a[2] = nbend(d[0], d[1]);
        cross(d[0], d[1], n[1]);
        normalise(n[1]);

        a[1] = nbend(n[0], n[1]);
        if (dot(n[0], d[1]) < 0.)
            a[1] *= -1;

        cl.push_back(btb2cl(a));
    }

    cl.push_back('R');

    for (i=1; i<len; i++)
    {
        if (brk[i]=='1')
        {
            cl[i-1] = 'R';
            cl[i] = 'R';
            if (i+1 < len)
                cl[i+1] = 'R';
        }
    }

    return cl;
}
}
