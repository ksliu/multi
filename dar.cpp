#include "dar.h"
#include <iostream>
#include <cstdlib>

#define DAR_RANGE_CHECK

Dar::Dar(int n): p(0), s(0)
{
    if(n > 0)
    {
        p = new double [n][3];
        s = n;
    }
}
Dar::~Dar()
{
    delete [] p;
}
int Dar::len() const
{
    return s;
}
void Dar::reset(int n)
{
    delete [] p;
    if (n>0)
    {
        p = new double [n][3];
        s = n;
    }
    else
    {
        p = 0;
        s = 0;
    }
}

double * Dar::operator[] (int i)
{
#ifdef DAR_RANGE_CHECK
    if(i <0 || i>=s)
    {
        std::cout << "Dar range error " << std::endl;
        exit(1);
    }
#endif
    return p[i];
}

const double * Dar::operator[] (int i) const
{
#ifdef DAR_RANGE_CHECK
    if(i <0 || i>=s)
    {
        std::cout << "Dar range error " << std::endl;
        exit(1);
    }
#endif
    return p[i];
}
double (*Dar::ptr()) [3]
{
    return p;
}
const double (*Dar::ptr() const) [3]
{
    return p;
}
