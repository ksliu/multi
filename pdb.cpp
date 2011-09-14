#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <fstream>
#include <set>
#include "pdb.h"
#include "dar.h"
#include "bioinfo.h"

using namespace std;

PDB::PDB(const std::string &fn, char c)
{
    if (c == '_')
        cacheFirst(fn);
    else
        cacheRange(fn, toupper(c));
}
void PDB::cacheFirst(const std::string & fn)
{
    ifstream fs(fn.c_str());
    if (fs.fail())
    {
        cerr << "cannot open file " << fn << endl;
        exit(1);
    }

    string line;
    set<string> altLoc;

    while (getline(fs, line))
    {
        if (line.compare(0, 3, "TER") == 0 || line.compare(0, 3, "END") == 0) break;
        if ((line.compare(0, 6, "ATOM  ") == 0|| line.compare(0, 6, "HETATM") == 0))
        {
            if (line[16] == ' ' || altLoc.insert(line.substr(21, 6) + line.substr(12, 4)).second)
                domain.push_back(line);
        }
    }

    fs.close();
}
void PDB::cacheRange(const std::string &fn, char chain)
{
    ifstream fs(fn.c_str());
    if (fs.fail())
    {
        cerr << "cannot open file " << fn << endl;
        exit(1);
    }

    string line;
    set<string> altLoc;

    while (getline(fs, line))
        if ((line.compare(0,4, "ATOM")==0 || line.compare(0, 6, "HETATM")==0)
                && line[21] == chain )
            break;

    if (!fs)
    {
        cerr << "Chain " << chain << "was not found in file " << fn << endl;
        exit(1);
    }

    do
    {
        if (line.compare(0, 3, "TER") == 0 || line.compare(0, 3, "END") == 0) break;
        if (line.compare(0, 6, "ATOM  ") == 0|| line.compare(0, 6, "HETATM") == 0)
        {
            if (line[16] == ' ' || altLoc.insert(line.substr(21, 6) + line.substr(12, 4)).second)
                domain.push_back(line);
        }
    }
    while (getline(fs, line));

    fs.close();
}

void PDB::getAA(string &aa) const
{
    aa.clear();
    for (vector<string>::size_type i = 0; i != domain.size(); ++i)
    {
        if (domain[i].compare(12, 4, " CA ") == 0)
            aa.push_back(bio::aaConvert31(domain[i].substr(17, 3).c_str()));
    }
}
int PDB::getCaNumber() const
{
    int n = 0;
    for (vector<string>::size_type i = 0; i != domain.size(); ++i)
    {
        if (domain[i].compare(12, 4, " CA ") == 0)
            n++;
    }
    return n;
}
void PDB::getCa(Dar &ca) const
{
    int caNum = getCaNumber();
    ca.reset(caNum);
    int iCa = 0;
    for (vector<string>::size_type i = 0; i != domain.size(); ++i)
    {
        if (domain[i].compare(12, 4, " CA ") == 0)
        {
            for (int j = 0; j < 3; j++)
                ca[iCa][j] = atof(domain[i].substr(30 + 8 * j, 8).c_str());
            iCa++;
        }
    }
    if (iCa != caNum)
    {
        cerr << "mismatch in CA" << endl;
        exit(1);
    }
}

void PDB::getCaCb(Dar &ca, Dar &cb) const
{
    int caNum = getCaNumber();
    ca.reset(caNum);
    cb.reset(caNum);

    int i, iCa, iCb, delta;
    string lastResidueID, currentResidueID;

    iCa = iCb = 0;
    for (vector<string>::size_type ix = 0; ix != domain.size(); ++ix)
    {
        const string & line = domain[ix];
        currentResidueID = line.substr(21, 6);
        if ( currentResidueID != lastResidueID )
        {
            delta = iCa - iCb;
            if (delta == 1)
            {
                for (i = 0; i < 3; i++)
                    cb[iCb][i] = ca[iCb][i];
                iCb++;
            }
            else if (delta == -1)
            {
                iCb--;
            }
        }
        lastResidueID = currentResidueID;

        if (line.compare(12, 4, " CA ") == 0)
        {
            if (iCa < caNum)
            {
                for (i = 0; i < 3; i++)
                    ca[iCa][i] = atof(line.substr(30 + 8 * i, 8).c_str());
                iCa++;
            }
        }
        else if (line.compare(12, 4, " CB ") == 0)
        {
            if (iCb < caNum)
            {
                for (i = 0; i < 3; i++)
                    cb[iCb][i] = atof(line.substr(30 + 8 * i, 8).c_str());
                iCb++;
            }
        }
    }
    delta = iCa - iCb;
    if (delta == 1)
    {
        for (i = 0; i < 3; i++)
            cb[iCb][i] = ca[iCb][i];
        iCb++;
    }
    else if (delta == -1)
    {
        iCb--;
    }

    if (iCa != caNum || iCb != caNum)
    {
        cerr << "mismatch in CA, CB" << endl;
        exit(1);
    }
}

int PDB::bbAtomType(const string & line) const
{
    if (line.compare(12, 4, " N  ") == 0)
        return 0;
    if (line.compare(12, 4, " CA ") == 0)
        return 1;
    if (line.compare(12, 4, " C  ") == 0)
        return 2;
    if (line.compare(12, 4, " O  ") == 0)
        return 3;
    if (line.compare(12, 4, " CB ") == 0)
        return 4;
    return -1;
}

void PDB::getBB(Dar bb[5]) const
{
    const int NBBTYPE = 5;
    int caNum = getCaNumber();
    for (int i=0; i<NBBTYPE; i++)
        bb[i].reset(caNum);

    int iAt[NBBTYPE] = {0};
    string lastResidueID, currentResidueID;

    for (vector<string>::size_type ix = 0; ix != domain.size(); ++ix)
    {
        const string & line = domain[ix];
        currentResidueID = line.substr(21, 6);
        if ( currentResidueID != lastResidueID )
        {
            for (int atom=0; atom<NBBTYPE; atom++)
            {
                switch (iAt[1] - iAt[atom])
                {
                case 0:
                    break;
                case 1:
                    for (int j = 0; j < 3; j++)
                        bb[atom][iAt[atom]][j] = bb[1][iAt[atom]][j];
                    iAt[atom]++;
                    break;
                case -1:
                    iAt[atom]--;
                    break;
                default:
                    cout << "error: impossible "<< endl;
                    exit(1);
                }
            }
        }
        lastResidueID = currentResidueID;

        int type = bbAtomType(line);
        if (type != -1)
        {
            // when last residue contain atoms other than "CA"
            // space of these atoms will be not correct.
            if (iAt[type] < caNum)
            {
                for (int j = 0; j < 3; j++)
                    bb[type][iAt[type]][j]  = atof(line.substr(30 + 8 * j, 8).c_str());
                iAt[type]++;
            }
        }
    }
    for (int atom=0; atom<NBBTYPE; atom++)
    {
        switch (iAt[1] - iAt[atom])
        {
        case 0:
            break;
        case 1:
            for (int j = 0; j < 3; j++)
                bb[atom][iAt[atom]][j] = bb[1][iAt[atom]][j];
            iAt[atom]++;
            break;
        case -1:
            iAt[atom]--;
            break;
        default:
            cout << "error: impossible "<< endl;
            exit(1);
        }
    }
    // check
    for (int atom=0; atom< NBBTYPE; atom++)
    {
        if (iAt[atom] != caNum)
        {
            cout << "mismatch in getBB()" << endl;
            exit(1);
        }
    }
}
