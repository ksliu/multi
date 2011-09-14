#ifndef  BIOINFO_H
#define  BIOINFO_H

#include <string>

namespace bio
{
// amino acid: character and integer codes
char aaIntToChar(int);
int aaCharToInt(char);

// amino acid: 1-letter and 3-letters codes
char aaConvert31(const char *);
const char *aaConvert13(char);

// distance from CB to CA atom
double sideChainLengthInt(int);
double sideChainLengthChar(char);

// BLOSUM score for different amino acid pair
int blosum62Score(char, char);

// CLE score for different amino acid pair
int cleScore(char, char);

// consensus CLE code for a given CLE set
char cleConsensus(const std::string &s);
}

#endif
