#include <iostream>
#include <cmath>
#include <cstdlib>
#include <cassert>
#include "geometry.h"
namespace
{
const double DMIN = 1.e-9;
}
namespace geom
{
bool fit_helix(const double p[][3], double *n, double *b)
{
    int i;
    double s, x[3], y[3], t[3];
    for (i = 0; i < 3; i++)
    {
        x[i] = p[2][i] - 2 * p[1][i] + p[0][i];
        y[i] = p[3][i] - 2 * p[2][i] + p[1][i];
        t[i] = p[3][i] - p[0][i];
    }

    cross(x, y, n);
    s = len(n);

    // These 4 points do not seem coplanar.
    if (s > 1e-6)
    {
        for (i = 0; i < 3; i++)
            n[i] /= s;
        *b = dot(t, n) / 3;
        if (*b < 0)
        {
            *b *= -1;
            for (i = 0; i < 3; i++)
                n[i] *= -1;
        }
        return true;
    }
    else
    {
        *b = 0;
        for (i = 0; i < 3; i++)
            n[i] = y[i] - x[i];
        normalize(n);
        std::cerr << "Warning: coplanar detected" << std::endl;
        return false;
    }
}

void H_angles(const double *ha, const double *mid, const double *hb, double *c)
{
    c[0] = bend(ha, mid);
    c[1] = torsion(ha, mid, hb);
    c[2] = bend(mid, hb);
}

double dot(const double *a, const double *b)
{
    double s = 0;
    int i;
    for (i = 0; i < 3; i++)
        s += a[i] * b[i];
    return s;
}

void cross(const double *a, const double *b, double *c)
{
    c[0] = a[1] * b[2] - a[2] * b[1];
    c[1] = a[2] * b[0] - a[0] * b[2];
    c[2] = a[0] * b[1] - a[1] * b[0];
}

double len(const double *a)
{
    int i;
    double s = 0;
    for (i = 0; i < 3; i++)
        s += a[i] * a[i];
    return sqrt(s);
}

double dist(const double *a, const double *b)
{
    int i;
    double s = 0;
    for (i = 0; i < 3; i++)
        s += (a[i] - b[i]) * (a[i] - b[i]);
    return sqrt(s);

}

double bend(const double *a, const double *b)
{
    return acos(dot(a, b) / len(a) / len(b));
}

// a, b already normalized
double nbend(const double *a, const double *b)
{
    return acos(dot(a, b));
}

double torsion(const double *a, const double *b, const double *c)
{
    double u[3], v[3];
    cross(a, b, u);
    cross(b, c, v);
    return (dot(u, c) < 0 ? -1.0 : 1.0) * bend(u, v);
}

void normalize(double *a)
{
    double s = len(a);
    int i;
    if (s > DMIN)
    {
        for (i = 0; i < 3; i++)
            a[i] /= s;
    }
}
bool IsZero(double x)
{
    return (x > -1. * DMIN && x < DMIN);
}
double mdet(const double m[3][3])
{
    return m[0][0] * m[1][1] * m[2][2] + m[0][1] * m[1][2] * m[2][0] + m[1][0]
            * m[2][1] * m[0][2] - m[0][0] * m[1][2] * m[2][1] - m[0][2]
            * m[1][1] * m[2][0] - m[0][1] * m[1][0] * m[2][2];
}

}

namespace
{
using geom::IsZero;
int cubic_roots(double c[4], double s[3])
{
    int i, num;
    double sub;
    double A, B, C;
    double A2, p, q;
    double p3, D;

    /* normal form: x^3 + Ax^2 + Bx + C = 0 */

    A = c[2] / c[3];
    B = c[1] / c[3];
    C = c[0] / c[3];

    /*  substitute x = y - A/3 to eliminate quadric term:
  x^3 +px + q = 0 */

    A2 = A * A;
    p = 1.0 / 3 * (-1.0 / 3 * A2 + B);
    q = 1.0 / 2 * (2.0 / 27 * A * A2 - 1.0 / 3 * A * B + C);

    /* use Cardano's formula */

    p3 = p * p * p;
    D = q * q + p3;

    if (IsZero(D))
    { /* one triple solution */
        if (IsZero(q))
        {
            s[0] = 0;
            num = 1;
        }
        else
        { /* one single and one double solution */
            double u = cbrt(-q);
            s[0] = 2 * u;
            s[1] = -u;
            num = 2;
        }
    }
    else if (D < 0)
    { /* Casus irreducibilis: three real solutions */

        double phi = 1.0 / 3 * acos(-q / sqrt(-p3));
        double t = 2 * sqrt(-p);

        s[0] = t * cos(phi);
        s[1] = -t * cos(phi + M_PI / 3);
        s[2] = -t * cos(phi - M_PI / 3);
        num = 3;
    }
    else
    { /* one real solution */

        double sqrt_D = sqrt(D);
        double u = cbrt(sqrt_D - q);
        double v = -cbrt(sqrt_D + q);

        s[0] = u + v;
        num = 1;
    }

    /* resubstitute */

    sub = 1.0 / 3 * A;

    for (i = 0; i < num; ++i)
        s[i] -= sub;

    return num;
}
/*
 * Find the eigen values and vectors for the matrix m.
 */
void eigen_values(double m[3][3], double values[3], double vectors[3][3])
{
    double a1 = m[0][0], b1 = m[0][1], c1 = m[0][2];
    double a2 = m[1][0], b2 = m[1][1], c2 = m[1][2];
    double a3 = m[2][0], b3 = m[2][1], c3 = m[2][2];
    double c[4];
    double x, y, z, norm;
    double roots[3], l;
    int nroots, iroot;
    int i;

    /*
  * Expanding the characteristic equation of a 3x3 matrix...
  * Maple tells us the following is the cubic in l
  */

    c[0] = a1 * (b2 * c3 - b3 * c2) + c1 * (a2 * b3 - a3 * b2) - b1 * (a2 * c3
                                                                       - a3 * c2);
    c[1] = (a1 * (-b2 - c3) - b2 * c3 + b3 * c2 + c1 * a3 + b1 * a2);
    c[2] = (a1 + b2 + c3);
    c[3] = -1.;

    nroots = cubic_roots(c, roots);

    /* Degenerate roots are not returned individually. */

    for (i = 0; i < nroots; i++)
    {
        iroot = i > nroots ? nroots : i;
        l = roots[iroot];
        values[i] = l;
        /*
   * Find the eigen vectors by solving pairs of the
   * three simultaneous equations.  From `Mathematical Methods
   * in Science and Engineering', Heiding, pg.19
   *
   * Sometimes we get x = y = z = 0.0, so try the other two
   * pairs of equations and hope that one of them gives a solution.
   */

        x = b1 * c2 - (b2 - l) * c1;
        y = -((a1 - l) * c2 - a2 * c1);
        z = ((a1 - l) * (b2 - l) - a2 * b1);

        if (IsZero(x) && IsZero(y) && IsZero(z))
        {
            x = b1 * (c3 - l) - b3 * c1;
            y = -((a1 - l) * (c3 - l) - a3 * c1);
            z = ((a1 - l) * b3 - a3 * b1);

            if (IsZero(x) && IsZero(y) && IsZero(z))
            {
                x = (b2 - l) * (c3 - l) - b3 * c2;
                y = -(a2 * (c3 - l) - a3 * c2);
                z = (a2 * b3 - a3 * (b2 - l));
                if (IsZero(x) && IsZero(y) && IsZero(z))
                {
                    std::cerr << "eigen: no solution for eigen vector " << i
                              << std::endl;
                    exit(1);
                }
            }
        }

        norm = sqrt(x * x + y * y + z * z);

        if (!IsZero(norm))
        {
            vectors[i][0] = x / norm;
            vectors[i][1] = y / norm;
            vectors[i][2] = z / norm;
        }
    }
}
}

namespace geom
{
//void opt_rotation(double x[][3], double y[][3], int n, double r[3][3]);
//void svd(double x[][3], double y[][3], int n, double r[3][3]);
void kabsch(double x[][3], double y[][3], int n, double r[3][3])
{
    double m[3][3], m2[3][3];
    double eig[3], a[3][3], b[3][3];
    int i, j, k;
    for (i = 0; i < 3; i++)
        for (j = 0; j < 3; j++)
            for (m[i][j] = 0.0, k = 0; k < n; k++)
                m[i][j] += x[k][i] * y[k][j];

    for (i = 0; i < 3; i++)
        for (j = 0; j < 3; j++)
            for (m2[i][j] = 0.0, k = 0; k < 3; k++)
                m2[i][j] += m[i][k] * m[j][k];
    /*
  * Get the eigenvalues and vectors.
  * Reform a[2] as cross product of a[0] and a[1] to ensure
  * right handed system.
  */

    eigen_values(m2, eig, a);
    cross(a[0], a[1], a[2]);

    /* Transform first two eigenvectors and normalise them. */
    for (i = 0; i < 2; i++)
    {
        for (j = 0; j < 3; j++)
            for (b[i][j] = 0.0, k = 0; k < 3; k++)
                b[i][j] += a[i][k] * m[k][j];
        normalize(b[i]);
    }

    /* Make right handed set. */
    cross(b[0], b[1], b[2]);

    for (i = 0; i < 3; i++)
        for (j = 0; j < 3; j++)
            for (r[i][j] = 0.0, k = 0; k < 3; k++)
                r[i][j] += a[k][i] * b[k][j];
}
}
////////////////////
//
//void opt_rotation (double x[][3], double y[][3], int n, double r[3][3])
//{
//    int i, j, k, i0, i1, i2;
//    double det, a, b, c, z[3], m[3][3], m2[3][3], eig[3], u[3][3], v[3][3];
//
//    for (i = 0; i < 3; i++)
//	for (j = 0; j < 3; j++)
//	    for (m[i][j] = 0., k = 0; k < n; k++)
//		m[i][j] += x[k][i] * y[k][j];
//    det = mdet (m);
//    //debug
//    if (det < 0)
//	std::cerr << "Warning: negative determinate" << det <<std::endl;
//    //debug
//    for (i = 0; i < 3; i++)
//	for (j = 0; j < 3; j++)
//	    for (m2[i][j] = 0., k = 0; k < 3; k++)
//		m2[i][j] += m[i][k] * m[j][k];	//M2 = M*M^T
//
//    // coefficients of eigenpolynomial of m2
//    a = -(m2[0][0] + m2[1][1] + m2[2][2]);
//    b = m2[0][0] * (m2[1][1] + m2[2][2]) - m2[0][1] * m2[1][0]
//	- m2[0][2] * m2[2][0] + m2[2][2] * m2[1][1] - m2[2][1] * m2[1][2];
//    c = -det * det;
//
//    cubic (a, b, c, eig);		// eig[3]: eigenvalues
//
//    for (i = 0; i < 3; i++)
//    {
//	for (a = eig[i], j = 0; j < 3; j++)
//	{
//	    i0 = j % 3;
//	    i1 = (j + 1) % 3;
//	    i2 = (j + 2) % 3;
//
//	    z[i2] = (m2[i0][i0] - a) * (m2[i1][i1] - a) - m2[i0][i1] * m2[i1][i0];
//	    if (fabs (z[i2]) < 1.e-3)
//		continue;
//	    else
//	    {
//		z[i0] = m2[i0][i1] * m2[i1][i2] - m2[i0][i2] * (m2[i1][i1] - a);
//		z[i1] = m2[i0][i2] * m2[i1][i0] - (m2[i0][i0] - a) * m2[i1][i2];
//		break;
//	    }
//	}
//	a = sqrt (z[0] * z[0] + z[1] * z[1] + z[2] * z[2]);
//	for (j = 0; j < 3; j++)
//	    u[j][i] = z[j] / a;	// u[3][3]: column of eigenvectors
//    }
//
//    for (i = 0; i < 3; i++)
//	eig[i] = sqrt (eig[i]);
//
//    for (i = 0; i < 3; i++)
//	for (j = 0; j < 3; j++)
//	    for (v[i][j] = 0., k = 0; k < 3; k++)
//		v[i][j] += u[k][i] * m[k][j] / eig[i];
//
//    for (i = 0; i < 3; i++)
//	for (j = 0; j < 3; j++)
//	{
//	    for (r[i][j] = 0., k = 0; k < 3; k++)
//		r[i][j] += u[i][k] * v[k][j];
//	    //r[i][j] += (det > 0. ? 1. : -1.) * u[i][2] * v[2][j];
//	}
//}
//
//void svd (double x[][3], double y[][3], int n, double r[3][3])
//{
//    int i, j, k, i0, i1, i2;
//    //double a,b,c, det, z1[3],xx,m[3][3], m2[3][3];
//    double det, a, b, c, w, z1[3], xx, m[3][3], m2[3][3], eig[3], u[3][3],
//	   v[3][3];
//    for (i = 0; i < 3; i++)
//	for (j = 0; j < 3; j++)
//	    for (m[i][j] = 0., k = 0; k < n; k++)	//%%%
//		m[i][j] += x[k][i] * y[k][j];	//xy[][]:correspondence
//    det = mdet (m);		//%%%%%
//    for (i = 0; i < 3; i++)
//	for (j = 0; j < 3; j++)
//	    for (m2[i][j] = 0., k = 0; k < 3; k++)
//		m2[i][j] += m[i][k] * m[j][k];	//M2 = M*M^T
//    //svd_u(m2); 
//    a = -(m2[0][0] + m2[1][1] + m2[2][2]);
//    b =
//	m2[0][0] * (m2[1][1] + m2[2][2]) - m2[0][1] * m2[1][0] -
//	m2[0][2] * m2[2][0] + m2[2][2] * m2[1][1] - m2[2][1] * m2[1][2];
//    c = -det * det;		//-mdet(m);//
//    cubic (a, b, c, eig);
//
//    for (i = 0; i < 3; i++)
//    {
//	for (xx = eig[i], j = 0; j < 3; j++)
//	{			//%%% these m2
//	    i0 = j % 3;
//	    i1 = (j + 1) % 3;
//	    i2 = (j + 2) % 3;
//	    w = (m2[i0][i0] - xx) * (m2[i1][i1] - xx) - m2[i0][i1] * m2[i1][i0];
//	    if (fabs (w) < 1.e-3)
//		continue;
//	    else
//	    {
//		z1[i0] = m2[i0][i1] * m2[i1][i2] - m2[i0][i2] * (m2[i1][i1] - xx);
//		z1[i1] = m2[i0][i2] * m2[i1][i0] - (m2[i0][i0] - xx) * m2[i1][i2];
//		z1[i2] = w;
//		break;
//	    }
//	}
//	//if(fabs(w)<1.e-8) {printf("fails 2!\n"); exit(0);}
//	xx = sqrt (z1[0] * z1[0] + z1[1] * z1[1] + z1[2] * z1[2]);
//	for (j = 0; j < 3; j++)
//	    u[j][i] = z1[j] / xx;
//    }
//
//    for (i = 0; i < 3; i++)
//	eig[i] = sqrt (eig[i]);
//    //printf("eig = %6.2f %6.2f %6.2f\n",eig[0],eig[1],eig[2]);//
//    for (i = 0; i < 3; i++)
//	for (j = 0; j < 3; j++)
//	    for (v[i][j] = 0., k = 0; k < 3; k++)
//		v[i][j] += u[k][i] * m[k][j] / eig[i];	// m, original 
//    for (i = 0; i < 3; i++)
//	for (j = 0; j < 3; j++)
//	{
//	    for (r[i][j] = 0., k = 0; k < 2; k++)
//		r[i][j] += u[i][k] * v[k][j];
//	    r[i][j] += (det > 0. ? 1. : -1.) * u[i][2] * v[2][j];
//	}
//}//R=USV^T,|C|=|R^TM|=|M|=|VDSV^T|=|D|*|S|, |D|>0,->s_z =sgn(|C|)-> R=USV^T.
//void cubic (double a, double b, double c, double eig[3])
//{
//    double w, q, q1, r, t, a1;
//    t = a * a;
//    if ((w = t - 3. * b) < 0.)
//    {
//	std::cerr << "cubic failed 1" << std::endl;
//	exit (1);
//    }
//    q1 = sqrt (w);
//    q = q1 * q1 * q1;
//    q1 /= 3.;
//    a1 = -a / 3.;
//    r = a * t - 4.5 * a * b + 13.5 * c;
//    if (q1 < 1.e-9 || r > q)
//    {
//	std::cerr << "cubic failed 2" << std::endl;
//	exit (1);
//    }
//    t = (fabs (r / q) > 1. ? q / r : r / q);
//    t = acos (t) / 3.;		//%%%% acos(r/q)
//    eig[0] = a1 - 2. * q1 * cos (t);
//    eig[1] = a1 - 2. * q1 * cos (t + 2.094395);
//    eig[2] = a1 - 2. * q1 * cos (t - 2.094395);
//}
////sort a[] in descending order, record the rank of a[i] as: ix[rank]=i 
//void ishell(int a[], int n, int ix[])//simple sorting
//{
//    int i,j,inc,id;
//    int v;
//    for(i=0;i<n;i++) ix[i]=i; 
//    inc=1; 
//    do { inc *= 3; inc++; } while(inc <= n);
//    do 
//    { 
//	inc /= 3;
//	for(i=inc;i<n;i++) 
//	{ 
//	    id=ix[i];
//	    v=a[id];
//	    j=i;
//	    while(a[ix[j-inc]] < v) 
//	    { //straight insertion.
//		ix[j]=ix[j-inc]; 
//		j -= inc;
//		if (j < inc) break; 
//	    }
//	    ix[j]=id; 
//	}
//    } while(inc > 1);
//}
