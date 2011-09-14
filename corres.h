/*
 * =====================================================================================
 *
 *       Filename:  corres.h
 *
 *    Description:  rotation between two point set controlled by alignment
 *
 *        Created:  05/10/2011 10:27:47 PM
 *
 *         Author:  Ke-song LIU (ksliu), ksliu10@gmail.com
 *
 * =====================================================================================
 */
#ifndef  CORRES_H
#define  CORRES_H
#include "protein.h"
#include <string>

class corres
{
public:
    corres();

    void set(int a[], double sub[][3], int sl, double que[][3], int ql);

    void full_strr(const std::string &rs, const std::string &rq,
                   std::string &rs_out, std::string &rq_out) const;

    void fullAlignString(const std::string &rs, const std::string &rq,
                         double cutoff1, double cutoff2, std::string &rs_out,
                         std::string &rq_out, std::string & mid) const;

    void concise_strr(const std::string &rs, const std::string &rq,
                      std::string &rs_out, std::string &rq_out, int threshold) const;

    void minAlignPair(const std::string &rs, const std::string &rq,
                      std::string &rs_out, std::string &rq_out, std::string & summary,
                      int threshold) const;

    void get_rotate_matrix(double rot[][3]) const;

    void ontop();
    bool increased() const;
    bool colinear(int k1, int k2) const; // greedy colinearity check: is colinear?

    double rmsd() const;
    double TMscore() const;
    int aligned_size() const;

    bool abs_dev(int k1, int k2, double err) const;
    bool sq_dev(int k1, int k2, double err) const;
    double sq_dist(int k1, int k2) const;

private:
    static const int N = 3000;
    corres(const corres &);
    corres operator=(const corres &);
    int
    extract(double px[][3], double xc[3], double py[][3], double yc[3]) const;

    double subject[N][3], query[N][3], shadow[N][3], rotation[3][3];
    int sublen, quelen, *ali;
};

#endif  /*CORRES_H*/
