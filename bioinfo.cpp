#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cstdlib>
#include <cstring>
#include <set>
#include <cassert>
#include "bioinfo.h"
using namespace std;

namespace
{
// PDB-style ATOM lines parser: readPDB()
// update to PDB format version 3.2
// COLUMNS        DATA  TYPE    FIELD        DEFINITION
// -------------------------------------------------------------------------------------
//  1 -  6        Record name   "ATOM  "
//  7 - 11        Integer       serial       Atom  serial number.
// 13 - 16        Atom          name         Atom name.
// 17             Character     altLoc       Alternate location indicator.
// 18 - 20        Residue name  resName      Residue name.
// 22             Character     chainID      Chain identifier.
// 23 - 26        Integer       resSeq       Residue sequence number.
// 27             AChar         iCode        Code for insertion of residues.
// 31 - 38        Real(8.3)     x            Orthogonal coordinates for X in Angstroms.
// 39 - 46        Real(8.3)     y            Orthogonal coordinates for Y in Angstroms.
// 47 - 54        Real(8.3)     z            Orthogonal coordinates for Z in Angstroms.
// 55 - 60        Real(6.2)     occupancy    Occupancy.
// 61 - 66        Real(6.2)     tempFactor   Temperature  factor.
// 77 - 78        LString(2)    element      Element symbol, right-justified.
// 79 - 80        LString(2)    charge       Charge  on the atom.

/* e.g.
 ATOM    173  CA  ASN A  21      27.678  17.546  -5.580  1.00 12.02           C
 */
struct pdb_atom
{
    char altLoc;
    char resName[4];
    char chainID;
    int resSeq;
    char iCode;
    string signature;
    double x[3];
    pdb_atom(const string &line);

};
pdb_atom::pdb_atom(const string &line)
{
    int i;
    altLoc = line[16];
    for (i = 0; i < 3; i++)
        resName[i] = line[17 + i];
    resName[i] = '\0';
    chainID = line[21];
    resSeq = atoi(line.substr(22, 4).c_str());
    iCode = line[26];

    signature = line.substr(12, 15);
    signature.erase(4, 1);

    for (i = 0; i < 3; i++)
        x[i] = atof(line.substr(30 + 8 * i, 8).c_str());
}
}

namespace bio
{
int readPDB(const char *fn, string &aa, double ca[][3], int capacity)
{

    // domain of a protein. Note: for altered residues, only the first one
    // is recorded, others are discarded.
    vector<pdb_atom> domain;

    // first appeared alternate residue
    set<string> alt_aa;

    // `insert_status.second' is false when an alternate residue is inserted repeatedly
    pair<set<string>::iterator, bool> insert_status;

    ifstream fs(fn);
    if (fs.fail())
    {
        cerr << "cannot open file " << endl;
        exit(1);
    }

    string line;
    bool done, caAtom;

    while (getline(fs, line))
    {
        done = line.compare(0, 3, "TER") == 0 || line.compare(0, 3, "END") == 0;

        if (done)
            break;

        caAtom = (line.compare(0, 6, "ATOM  ") == 0 || line.compare(0, 6,
                                                                    "HETATM") == 0) && line.compare(12, 4, " CA ") == 0;

        if (!caAtom)
            continue;

        pdb_atom cur(line);

        if (cur.altLoc == ' ')
        {
            domain.push_back(cur);
        }
        else
        {
            insert_status = alt_aa.insert(cur.signature);
            if (insert_status.second) // insert successfully
                domain.push_back(cur);
        }

    }

    fs.close();

    int length = static_cast<int> (domain.size());

    if (capacity < length)
    {
        cerr << "Array storage for coordinates is not enough!" << endl;
        exit(1);
    }
    aa.clear();
    for (int i = 0; i < length; i++)
    {
        for (int j = 0; j < 3; j++)
            ca[i][j] = domain[i].x[j];
        aa.push_back(aaConvert31(domain[i].resName));
    }
    return length;
}
// AA letters x and y are in upper case.
int blosum62Score(char x, char y)
{
    // BLOSUM Clustered Scoring Matrix in 1/2 Bit Units
    // Cluster Percentage: >= 62
    // @
    // http://www.ncbi.nlm.nih.gov/Class/FieldGuide/BLOSUM62.txt
    static const int BLOSUM62[24][24] =
    {
	//A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  *
        { 4, -1, -2, -2, 0, -1, -1, 0, -2, -1, -1, -1, -1, -2, -1, 1, 0,
          -3, -2, 0, -2, -1, 0, -4 },//A
        { -1, 5, 0, -2, -3, 1, 0, -2, 0, -3, -2, 2, -1, -3, -2, -1, -1, -3,
          -2, -3, -1, 0, -1, -4 },//R
        { -2, 0, 6, 1, -3, 0, 0, 0, 1, -3, -3, 0, -2, -3, -2, 1, 0, -4, -2,
          -3, 3, 0, -1, -4 },//N
        { -2, -2, 1, 6, -3, 0, 2, -1, -1, -3, -4, -1, -3, -3, -1, 0, -1,
          -4, -3, -3, 4, 1, -1, -4 },//D
        { 0, -3, -3, -3, 9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1,
          -2, -2, -1, -3, -3, -2, -4 },//C
        { -1, 1, 0, 0, -3, 5, 2, -2, 0, -3, -2, 1, 0, -3, -1, 0, -1, -2,
          -1, -2, 0, 3, -1, -4 },//Q
        { -1, 0, 0, 2, -4, 2, 5, -2, 0, -3, -3, 1, -2, -3, -1, 0, -1, -3,
          -2, -2, 1, 4, -1, -4 },//E
        { 0, -2, 0, -1, -3, -2, -2, 6, -2, -4, -4, -2, -3, -3, -2, 0, -2,
          -2, -3, -3, -1, -2, -1, -4 },//G
        { -2, 0, 1, -1, -3, 0, 0, -2, 8, -3, -3, -1, -2, -1, -2, -1, -2,
          -2, 2, -3, 0, 0, -1, -4 },//H
        { -1, -3, -3, -3, -1, -3, -3, -4, -3, 4, 2, -3, 1, 0, -3, -2, -1,
          -3, -1, 3, -3, -3, -1, -4 },//I
        { -1, -2, -3, -4, -1, -2, -3, -4, -3, 2, 4, -2, 2, 0, -3, -2, -1,
          -2, -1, 1, -4, -3, -1, -4 },//L
        { -1, 2, 0, -1, -3, 1, 1, -2, -1, -3, -2, 5, -1, -3, -1, 0, -1, -3,
          -2, -2, 0, 1, -1, -4 },//K
        { -1, -1, -2, -3, -1, 0, -2, -3, -2, 1, 2, -1, 5, 0, -2, -1, -1,
          -1, -1, 1, -3, -1, -1, -4 },//M
        { -2, -3, -3, -3, -2, -3, -3, -3, -1, 0, 0, -3, 0, 6, -4, -2, -2,
          1, 3, -1, -3, -3, -1, -4 },//F
        { -1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4, 7, -1,
          -1, -4, -3, -2, -2, -1, -2, -4 },//P
        { 1, -1, 1, 0, -1, 0, 0, 0, -1, -2, -2, 0, -1, -2, -1, 4, 1, -3,
          -2, -2, 0, 0, 0, -4 },//S
        { 0, -1, 0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1, 1, 5,
          -2, -2, 0, -1, -1, 0, -4 },//T
        { -3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1, 1, -4, -3,
          -2, 11, 2, -3, -4, -3, -2, -4 },//W
        { -2, -2, -2, -3, -2, -1, -2, -3, 2, -1, -1, -2, -1, 3, -3, -2, -2,
          2, 7, -1, -3, -2, -1, -4 },//Y
        { 0, -3, -3, -3, -1, -2, -2, -3, -3, 3, 1, -2, 1, -1, -2, -2, 0,
          -3, -1, 4, -3, -2, -1, -4 },//V
        { -2, -1, 3, 4, -3, 0, 1, -1, 0, -3, -4, 0, -3, -3, -2, 0, -1, -4,
          -3, -3, 4, 1, -1, -4 },//B
        { -1, 0, 0, 1, -3, 3, 4, -2, 0, -3, -3, 1, -1, -3, -1, 0, -1, -3,
          -2, -2, 1, 4, -1, -4 },//Z
        { 0, -1, -1, -1, -2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -2, 0, 0,
          -2, -1, -1, -1, -1, -1, -4 },//X
        { -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4,
          -4, -4, -4, -4, -4, -4, -4, 1 }, //*
    };
    // ARNDCQEGHILKMFPSTWYVBZX*
    //A, B, C, D, E,  F, G, H, I,  J,  K,  L,  M, N, O,  P, Q, R,  S,  T,  U,  V,  W,  X,  Y, Z,
    static const int INDEX[26] =
    { 0, 20, 4, 3, 6, 13, 7, 8, 9, -1, 11, 10, 12, 2, -1, 14, 5, 1, 15, 16, -1,
      19, 17, 22, 18, 21 };

    if (x < 'A' || x > 'Z' || y < 'A' || y > 'Z')
    {
        cerr << "Invalid AA letter" << endl;
        cerr << ":" << int(x) << ":" << int(y) << ":" << endl;
        cerr << ":" << char(x) << ":" << char(y) << ":" << endl;
        exit(1);
    }

    int p = INDEX[x - 'A'];
    int q = INDEX[y - 'A'];

    // Invalid AA letters are mapped to `-1'
    if (p == -1 || q == -1)
    {
        cerr << "Invalid AA letter" << endl;
        cerr << ":" << int(x) << ":" << int(y) << ":" << endl;
        cerr << ":" << char(x) << ":" << char(y) << ":" << endl;
        exit(1);
    }
    return BLOSUM62[p][q];
}

int cleScore(char cle_1, char cle_2)
{
    // 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17
    static const int CLESUM[18][18] =
    {
	{ 71, 19, 12, -19, -27, -22, -7, -43, -39, -22, -17, -12, -2, 10, 23, 34,
          13, -900, },
	{ 19, 49, 5, 11, 13, 4, 12, -94, -79, -55, -49, -13, -14, -11, -12, 41, 10,
          -900, },
        { 12, 5, 51, 19, 2, 19, -6, -75, -59, -42, -32, -2, -12, -5, 1, 4,
          27, -900, },
        { -19, 11, 19, 49, 21, 20, -32, -122, -105, -87, -81, -24, -50,
          -45, -43, -11, 13, -900, },
        { -27, 13, 2, 21, 34, 24, -24, -125, -112, -91, -83, -23, -47, -43,
          -49, -6, -8, -900, },
        { -22, 4, 19, 20, 24, 48, -6, -106, -95, -73, -67, -18, -34, -32,
          -31, -2, 0, -900, },
        { -7, 12, -6, -32, -24, -6, 68, -49, -40, -21, -11, 27, 21, -7, -2,
          3, 8, -900, },
        { -43, -94, -75, -122, -125, -106, -49, 25, 14, 15, 8, -61, -2,
          -32, -54, -58, -87, -900, },
        { -39, -79, -59, -105, -112, -95, -40, 14, 51, 12, 17, -48, -4,
          -16, -37, -55, -69, -900, },
        { -22, -55, -42, -87, -91, -73, -21, 15, 12, 38, 16, -31, 17, -1,
          -23, -33, -43, -900, },
        { -17, -49, -32, -81, -83, -67, -11, 8, 17, 16, 51, 0, 14, 28, 5,
          -35, -24, -900, },
        { -12, -13, -2, -24, -23, -18, 27, -61, -48, -31, 0, 71, 4, 5, -5,
          -14, 24, -900, },
        { -2, -14, -12, -50, -47, -34, 21, -2, -4, 17, 14, 4, 59, 8, 5, 7,
          -7, -900, },
        { 10, -11, -5, -45, -43, -32, -7, -32, -16, -1, 28, 5, 8, 89, 14,
          -4, 31, -900, },
        { 23, -12, 1, -43, -49, -31, -2, -54, -37, -23, 5, -5, 5, 14, 102,
          2, -13, -900, },
        { 34, 41, 4, -11, -6, -2, 3, -58, -55, -33, -35, -14, 7, -4, 2, 64,
          6, -900, },
        { 13, 10, 27, 13, -8, 0, 8, -87, -69, -43, -24, 24, -7, 31, -13, 6,
          88, -900, },
        { -900, -900, -900, -900, -900, -900, -900, -900, -900, -900, -900,
          -900, -900, -900, -900, -900, -900, -900, }, };
    if (cle_1 < 'A' || cle_1 > 'R' || cle_2 < 'A' || cle_2 > 'R')
    {
        cerr << "Invalid cle code" << endl;
        exit(1);
    }
    return CLESUM[cle_1 - 'A'][cle_2 - 'A'];
}

char aaConvert31(const char *resName)
{
    static const char *MAP[] =
    { "ALA", "A", "VAL", "V", "LEU", "L", "ILE", "I", "CYS", "C", "MET", "M",
      "PRO", "P", "PHE", "F", "TYR", "Y", "TRP", "W", "ASP", "D", "ASN",
      "N", "GLU", "E", "GLN", "Q", "HIS", "H", "SER", "S", "THR", "T",
      "ARG", "R", "LYS", "K", "GLY", "G", };
    int i, n = sizeof(MAP) / sizeof(MAP[0]) / 2;
    for (i = 0; i < n; i++)
        if (!strcmp(resName, MAP[2 * i]))
            return MAP[2 * i + 1][0];
    return 'X';
}

const char *aaConvert13(char c)
{
    static const char * MAP[] =
    { "ALA", "VAL", "LEU", "ILE", "CYS", "MET", "PRO", "PHE", "TYR", "TRP",
      "ASP", "ASN", "GLU", "GLN", "HIS", "SER", "THR", "ARG", "LYS",
      "GLY", };
    // AVLICMPFYWDNEQHSTRKG
    // A   B   C   D   E   F   G   H   I   J   K   L   M   N   O   P   Q   R   S   T   U   V   W   X   Y   Z
    static int INDEX[26] =
    { 0, -1, 4, 10, 12, 7, 19, 14, 3, -1, 18, 2, 5, 11, -1, 6, 13, 17, 15, 16,
      -1, 1, 9, -1, 8, -1, };
    if (c < 'A' || c > 'Z')
    {
        cerr << "Invalid AA letter" << endl;
        cerr << ":" << char(c) << ":" << int(c) << ":" << endl;
        exit(1);
    }

    int i = INDEX[c - 'A'];

    if (i != -1)
        return MAP[i];
    else
        return "XAA";
}

int readDALI(const char *fn, string &query, string & subject, int ali[],
             int capacity)
{
    query.clear();
    subject.clear();

    ifstream fs(fn);
    if (fs.fail())
    {
        cerr << "file open error" << endl;
        exit(1);
    }

    int i, j, k, qs_len, qi, si;
    string line, eff_line, *ps;

    for (i = 0; getline(fs, line); i++)
    {
        j = i % 7;

        if (j == 1)
            ps = &query;
        else if (j == 3)
            ps = &subject;
        else
            continue;

        eff_line = line.substr(6, 60);
        for (k = 0; k < 60; k++)
            if (eff_line[k] == ' ')
                break;

        if (k == 60)
            *ps += eff_line;
        else
            *ps += eff_line.substr(0, k);
    }

    fs.close();

    assert(query.size()==subject.size());

    qs_len = subject.size();

    //ali.assign(qs_len, -1);

    for (qi = si = i = 0; i < qs_len; i++)
    {
        if (subject[i] != '-' && query[i] != '-')
            ali[si] = qi;
        else
            ali[si] = -1;

        if (subject[i] != '-')
            si++;

        if (si > capacity)
        {
            cerr << "Array storage for alignment is not enough!" << endl;
            exit(1);
        }

        if (query[i] != '-')
            qi++;
    }
    //ali.erase(ali.begin()+si, ali.end());

    return si;
}
}

//Decomp_xyz::Decomp_xyz(const string line)
//{
//	if (line[13] == '!')
//	{
//		if (line[14] == '*')
//			line_status = new_chain;
//		else
//			line_status = chain_break;
//		return; // no data, nothing need constructed
//	}
//	else
//	{
//		line_status = normal;
//	}
//	int i;
//	for (i = 0; i < 5; i++)
//		pnum[i] = line[5 + i];
//	pnum[5] = '\0';
//	chainID = line[11];
//	aa = line[13] > 'Z' ? 'C' : line[13];
//	ss = line[16]; // 8 letters
//	elab = line[18];
//	cl = line[20];
//	for (i = 0; i < 3; i++)
//		ca[i] = atof(line.substr(22 + i * 7, 7).c_str()); // ca %7.1f
//	// reduce
//	if (ss == 'H' || ss == 'G' || ss == 'I')
//		ss = 'H';
//	else if (ss == 'E' || ss == 'B')
//		ss = 'E';
//	else
//		ss = 'C';
//}
