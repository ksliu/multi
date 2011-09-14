#ifndef GEOMETRY_H
#define GEOMETRY_H

namespace geom
{
/**
 * @brief fit a helix curve by 4 points
 * @param [in] p 4 points p[0], p[1], p[2], p[3]
 * @param [out] n direction vector of the helix axis
 * @param [out] b pitch of the helix
 * @return 4 points noncoplanar?
 */
bool fit_helix(const double p[][3], double *n, double *b);
//void fit_helix(const mat p[][3], double *n, double *b);
/**
 * @brief angles of a H-form
 * @param [in] ha a helix aixs
 * @param [in] mid vector from a to b
 * @param [in] hb the other helix aixs
 * @param [out] c bending, torsion, bending angles
 */
void H_angles(const double *ha, const double *mid, const double *hb, double *c);

// easy part
double dot(const double *a, const double *b);
void cross(const double *a, const double *b, double *c);
double len(const double *a);
double dist(const double *a, const double *b);
double bend(const double *a, const double *b);
double nbend(const double *a, const double *b);
double torsion(const double *a, const double *b, const double *c);
void normalize(double *a);
double mdet(const double m[3][3]);

//rotation part
// move y to x, by left multiply x by r
// x ~ ry
//       Kabsch (1976) superposition routines
void kabsch(double x[][3], double y[][3], int n, double r[3][3]);
}

#endif /* GEOMETRY_H */
