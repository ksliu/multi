//#include <iostream>
//#include <iomanip>
//#include <fstream>
//#include <algorithm>
//#include <vector>
//#include <string>
//#include <cmath>
//#include <cstdlib>
//#include "hfali.h"
//#include "bioinfo.h"
//#include "read_conf.h"

//using namespace std;

//const int hfali::SW, hfali::SC, hfali::N;

//hfali::hfali(const char *subject, const char *query)
//{
//    a = new protein(subject);
//    b = new protein(query);

//    adjacent_window_score(a->cl, b->cl, SW, bio::cleScore, SC, sq);
//    shave_adjacent_window_score(sq, SW / 2);
//    shave_adjacent_window_score(sq, SW / 2);
//    sort(sq.begin(), sq.end(), deSFPcmp);

//    cor.set(ali, a->ca, a->len, b->ca, b->len);
//}
//hfali::~hfali()
//{
//    delete b;
//    delete a;
//}

//void hfali::add(int fromBegin, int fromEnd, int toBegin, int toEnd,
//		int *outAlign)
//{
//    bool r1 = (fromBegin >= 0 && fromBegin < a->len) //
//            && (fromEnd >= 0 && fromEnd < a->len);
//    bool r2 = (toBegin >= 0 && toBegin < b->len) //
//            && (toEnd >= 0 && toEnd < b->len);

//    bool d = (fromBegin < fromEnd) && (fromEnd - fromBegin == toEnd - toBegin);

//    if (r1 && r2 && d)
//    {
//        int i = fromBegin, j = toBegin;
//        for (; i <= fromEnd; i++, j++)
//            outAlign[i] = j;
//    }
//    else
//    {
//        cerr << "range domain error" << endl;
//        exit(1);
//    }
//}
//void hfali::solve(int *a1, int *a2, int *b1, int *b2, IM method)
//{
//    method_log = method;

//    for (int i = 0; i < a->len; i++)
//        init_ali[i] = core_ali[i] = ali[i] = -1;

//    add(a1[0], a1[1], b1[0], b1[1], ali);
//    add(a2[0], a2[1], b2[0], b2[1], ali);

//    add(a1[0], a1[1], b1[0], b1[1], init_ali);
//    add(a2[0], a2[1], b2[0], b2[1], init_ali);

//    scon.assign(sq.size(), 0); //(is i-th sSFP irrelevant?), initial: scon[]=0

//    Para &conf = *Para::getPointer();
//    double c1 = conf.get_value("core_1"),//
//            c2 = conf.get_value("core_2"),//
//            e1 = conf.get_value("elong_1"), //
//            e2 = conf.get_value("elong_2");

//    switch (method)
//    {
//    case IM_once:
//        tune(c1, e1);
//        tune(c2, e2);
//        break;
//    case IM_zoom_in:
//        fill(0.8, 6.);
//        fill(1.0, 4.);
//        tune(c1, e1);
//        tune(c2, e2);
//        break;
//    default:
//        cerr << "unknown method" << endl;
//        break;
//    }
//}

//void hfali::output(ostream & fs) const
//{
//    cor.increased();
//    fs << fixed << setprecision(2);
//    fs << "subject: " << a->id << "  len=" << a->len << endl << "query:   "
//       << b->id << "  len=" << b->len << endl;

//    fs << "Method: base on initial aligned H-forms, method_num: #"
//       << method_log << endl;
//    fs << "Result: aligned_points=" << cor.aligned_size() << endl
//       << "aligned_ratio=" << 2.0 * cor.aligned_size() / (a->len + b->len)
//       << endl;
//    fs << "rmsd=" << cor.rmsd() << endl;
//    fs << "TMscore=" << cor.TMscore() << endl << endl;

//    // transformation
//    double rot[3][3], trans[3];
//    move(rot, trans);

//    int i, j;
//    fs << "rotation matrix" << endl;
//    for (i = 0; i < 3; i++)
//    {
//        for (j = 0; j < 3; j++)
//        {
//            fs << setw(8) << rot[i][j];
//        }
//        fs << endl;
//    }
//    fs << "translation vector" << endl;
//    for (i = 0; i < 3; i++)
//    {
//        fs << setw(8) << trans[i];
//    }
//    fs << endl;

//    fs << endl << "sequences" << endl;
//    fs << a->aa << endl << b->aa << endl;

//    // string representation
//    string ra, rb, sra, srb, mid;
//    // initial alignment
//    fs << "\ninitial H-form" << endl;
//    bool in = false;
//    for (i = 0; i < a->len; i++)
//    {
//        j = init_ali[i];
//        if (j != -1)
//        {
//            if (!in)
//            {
//                if (!ra.empty())
//                {
//                    ra.push_back('-');
//                    rb.push_back('-');
//                }
//                in = true;
//            }
//            ra.push_back(a->aa[i]);
//            rb.push_back(b->aa[j]);
//        }
//        else
//        {
//            in = false;
//        }
//    }

//    fs << ra << endl << rb << endl;

//    Para &conf = *Para::getPointer();
//    double c1 = conf.get_value("core_2"),//
//            c2 = conf.get_value("elong_2");

//    fs << "\nfull string representation ( " << ": " << "[0, " << c1 << ")   "
//       << ". " << "[" << c1 << ", " << c2 << ") )" << endl;

//    cor.fullAlignString(a->aa, b->aa, c1, c2, ra, rb, mid);
//    fs << ra << endl << mid << endl << rb << endl;

//    string summary;
//    const int min_w = 4;
//    fs << "\nconcise string representation (min_w=" << min_w << "): " << endl;

//    cor.minAlignPair(a->aa, b->aa, ra, rb, summary, min_w);
//    fs << ra << endl << rb << endl << summary << endl;
//}

//void hfali::fill(double fract, double abs_err)
//{
//    int i, k, nwr, k1, k2, part;

//    cor.ontop();
//    for (i = 0; i < a->len; i++)
//        ali[i] = -1;

//    part = static_cast<int> (fract * sq.size());
//    for (i = 0; i < part; i++)
//    {
//        if (scon[i])
//            continue;
//        k1 = sq[i].ia;
//        k2 = sq[i].ib;
//        // number of wrong points
//        for (nwr = k = 0; k < SW; k++)
//        {
//            if (cor.abs_dev(k1 + k, k2 + k, abs_err))
//            {
//                if (ali[k1 + k] == -1)
//                    ali[k1 + k] = k2 + k; //ignore confliction
//            }
//            else
//            {
//                nwr++;
//            }
//            if (nwr > SW / 2)
//            {
//                scon[i] = 1;
//                break;
//            }
//        }
//    }
//}

//void hfali::tune(double core, double elong)
//{
//    int i, j, k, k1, k2, ila[N];
//    double sq_core_err = core * core;
//    double sq_elong_err = elong * elong;

//    cor.ontop();
//    for (i = 0; i < a->len; i++)
//        ali[i] = -1;
//    for (i = 0; i < b->len; i++)
//        ila[i] = -1;

//    if (method_log == IM_once)
//    {
//        for (i = 0; i < a->len; i++)
//        {
//            j = init_ali[i];
//            if (j != -1)
//            {
//                ali[i] = j;
//                ila[j] = i;
//            }
//        }
//    }

//    //core// sq[] is sorted in descending order.
//    for (i = 0; i < static_cast<int> (sq.size()); i++)
//    {
//        if (scon[i])
//            continue;
//        for (k = 0; k < SW; k++)
//        {
//            k1 = sq[i].ia + k;
//            k2 = sq[i].ib + k;

//            if (ali[k1] == -1 && ila[k2] == -1 // first come, first served
//                    && cor.colinear(k1, k2) && cor.sq_dev(k1, k2, sq_core_err))
//            {
//                ali[k1] = k2;
//                ila[k2] = k1;
//            }
//        }
//    }
//    for (i = 0; i < a->len; i++)
//        core_ali[i] = ali[i];

//    //elongation right hand
//    for (i = 1; i < a->len; i++)
//    {
//        if (ali[i - 1] != -1 && ali[i] == -1)
//        {
//            j = ali[i - 1] + 1;
//            if (ila[j] == -1 && cor.sq_dev(i, j, sq_elong_err))
//            {
//                ali[i] = j;
//                ila[j] = i;
//            }
//        }
//    }
//    // elongation left hand
//    for (i = a->len - 1; i > 0; i--)
//    {
//        if (ali[i - 1] == -1 && ali[i] != -1)
//        {
//            j = ali[i];
//            if (ila[j - 1] == -1 && cor.sq_dev(i - 1, j - 1, sq_elong_err))
//            {
//                ali[i - 1] = j - 1;
//                ila[j - 1] = i - 1;
//            }
//        }
//    }

//    // pick the remaining
//    for (i = 0; i < static_cast<int> (sq.size()); i++)
//    {
//        if (scon[i])
//            continue;
//        for (k = 0; k < SW; k++)
//        {
//            k1 = sq[i].ia + k;
//            k2 = sq[i].ib + k;
//            if (ali[k1] == -1 && ila[k2] == -1 && cor.colinear(k1, k2)
//                    && cor.sq_dev(k1, k2, sq_elong_err))
//            {
//                ali[k1] = k2;
//                ila[k2] = k1;
//            }
//        }
//    }

//}

//void hfali::output_script(const char *fn) const
//{
//    int i, j, k;
//    double rot[3][3], trans[3];
//    move(rot, trans);

//    const int SEP = 3000;
//    double mb[N][3];

//    for (i = 0; i < b->len; i++)
//    {
//        for (j = 0; j < 3; j++)
//        {
//            mb[i][j] = trans[j];
//            for (k = 0; k < 3; k++)
//                mb[i][j] += rot[j][k] * b->ca[i][k];
//        }
//    }

//    ofstream fs(fn);
//    if (fs.fail())
//    {
//        cerr << "Warning: write file failed" << endl;
//        exit(1);
//    }
//    fs << fixed << setprecision(3);
//    // output header
//    fs << "load inline\n"
//          "select atomno<3000\n"
//          "wireframe .45\n"
//          "select none\n"
//          "select atomno>3000\n"
//          "wireframe .20\n"
//          "color white\n";
//    for (i = 0; i < a->len; i++)
//    {
//        if (ali[i] == -1)
//            continue;
//        fs << "select atomno=" << setw(4) << i << endl;
//        fs << "color red\n";
//        fs << "select atomno=" << setw(4) << SEP + ali[i] << endl;
//        fs << "color red\n";
//    }
//    fs << "select all\nexit\n";

//    // saved a
//    for (i = 0; i < a->len; i++)
//    {
//        fs << "ATOM  " << setw(5) << i << " " << " CA " << " " << setw(3)
//           << bio::aaConvert13(a->aa[i]) << " " << " " << setw(4) << i
//           << " " << "   ";
//        for (j = 0; j < 3; j++)
//            fs << setw(8) << a->ca[i][j];
//        fs << endl;
//    }
//    fs << "TER   " << endl;
//    for (i = 1; i < a->len; i++)
//        fs << "CONECT" << setw(5) << i - 1 << setw(5) << i << endl;

//    // new b
//    for (i = 0; i < b->len; i++)
//    {
//        fs << "ATOM  " << setw(5) << SEP + i << " " << " CA " << " " << setw(3)
//           << bio::aaConvert13(b->aa[i]) << " " << " " << setw(4) << i
//           << " " << "   ";
//        for (j = 0; j < 3; j++)
//            fs << setw(8) << mb[i][j];
//        fs << endl;
//    }
//    fs << "TER   " << endl;
//    for (i = 1; i < b->len; i++)
//        fs << "CONECT" << setw(5) << SEP + i - 1 << setw(5) << SEP + i << endl;

//    fs.close();
//}
//void hfali::move(double rot[][3], double trans[3]) const
//{
//    int i, j, k, n;
//    double ps[3], pq[3];

//    for (j = 0; j < 3; j++)
//        ps[j] = pq[j] = 0;

//    for (n = i = 0; i < a->len; i++)
//    {
//        k = ali[i];
//        if (k != -1)
//        {
//            for (j = 0; j < 3; j++)
//            {
//                ps[j] += a->ca[i][j];
//                pq[j] += b->ca[k][j];
//            }
//            n++;
//        }
//    }

//    for (j = 0; j < 3; j++)
//    {
//        ps[j] /= n;
//        pq[j] /= n;
//    }

//    cor.get_rotate_matrix(rot);
//    for (i = 0; i < 3; i++)
//    {
//        trans[i] = 0;
//        for (j = 0; j < 3; j++)
//        {
//            trans[i] += rot[i][j] * pq[j];
//        }
//    }
//    for (j = 0; j < 3; j++)
//        trans[j] = ps[j] - trans[j];
//}

//void hfali::adjacent_window_score(string a, string b, int width, int(*score)(
//                                      char, char), int lower_socre, vector<SFP> &q)
//{
//    int i, j, k, m, n, s, la = a.size(), lb = b.size();
//    vector<int> z(min(la, lb) + 1);
//    //vector<int> z(max(la, lb));

//    if (la < width || lb < width)
//    {
//        cerr << "code too short!" << endl;
//        exit(1);
//    }

//    q.clear();
//    // i, k ~ a
//    // j ~ b, fixed
//    for (k = 0; k <= la - width; k++) // stand on 'b', shift 'a'
//    {
//        z[0] = 0;
//        for (i = k, j = 0; i < la && j < lb; i++, j++)
//            z[j + 1] = z[j] + score(a[i], b[j]);

//        for (m = 0, n = width; n <= j; m++, n++)
//        {
//            s = z[n] - z[m];
//            if (s > lower_socre)
//                q.push_back(SFP(m + k, m, s));
//        }
//    }
//    // i, k ~ b
//    // j ~ a, fixed
//    for (k = 1; k <= lb - width; k++) // stand on 'a', shift 'b'. Set k=1 to exclude head-head
//    {
//        z[0] = 0;
//        for (i = k, j = 0; i < lb && j < la; i++, j++)
//            z[j + 1] = z[j] + score(b[i], a[j]);

//        for (m = 0, n = width; n <= j; m++, n++)
//        {
//            s = z[n] - z[m];
//            if (s > lower_socre)
//                q.push_back(SFP(m, m + k, s));
//        }
//    }
//}
//void hfali::window_score(string a, string b, int width,
//                         int(*score)(char, char), int lower_socre, vector<SFP> &q)
//{
//    int i, j, k, s, la = a.size(), lb = b.size();

//    if (la < width || lb < width)
//    {
//        cerr << "code too short!" << endl;
//        exit(1);
//    }
//    q.clear();

//    for (i = 0; i <= la - width; i++)
//        for (j = 0; j <= lb - width; j++)
//        {
//            for (s = k = 0; k < width; k++)
//                s += score(a[i + k], b[j + k]);
//            if (s > lower_socre)
//                q.push_back(SFP(i, j, s));
//        }
//}

//void hfali::shave_adjacent_window_score(vector<SFP> &q, int d)
//{
//    vector<SFP>::size_type i, j, ix, nn;
//    int t;

//    for (nn = i = 0; i != q.size(); i = j)
//    {
//        t = q[i].score;
//        ix = i;
//        for (j = i + 1; j != q.size(); ++j)
//        {
//            if (q[j].ia - q[i].ia >= d || (q[i].ia - q[i].ib != q[j].ia
//                                           - q[j].ib))
//                break; // get proper j
//            if (q[j].score > t) // find maximum in q[*].score
//            {
//                t = q[j].score;
//                ix = j;
//            }
//        }
//        q[nn++] = q[ix];
//    }
//    q.erase(q.begin() + nn, q.end());
//}
