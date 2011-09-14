#ifndef  BIOINFO_H
#define  BIOINFO_H
#include <string>
#include <vector>

#include "gencle.h"

namespace bio
{
int blosum62Score(char aa_1, char aa_2);
int cleScore(char cle_1, char cle_2);

char aaConvert31(const char *resName);

const char *aaConvert13(char c);

// read PDB ATOM style (astral) file
int readPDB(const char *fn, std::string &aa, double ca[][3], int capacity);

// read DALI alignment file
int readDALI(const char *fn, std::string &query, std::string & subject,
             int ali[], int capacity);

// read xyz file
//struct Decomp_xyz;
}

//namespace bio
//{
//struct Decomp_xyz
//{
//	Decomp_xyz(const std::string line);
//	enum status
//	{
//		chain_break, new_chain, normal
//	};
//	status line_status;
//	char pnum[6], chainID, aa, ss, elab, cl;
//	double ca[3];
//};
//}

#endif  /*BIOINFO_H*/
