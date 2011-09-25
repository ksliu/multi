#ifndef TRANSFORM_H
#define TRANSFORM_H

#include <vector>

class transform
{
public:
    transform();
    void run(double s[][3], int ls, double q[][3], int lq, const std::vector<int> &);

    void getPartial(int n);
    //    ~transform();

private:
    typedef double PA [][3];
    PA subject, query, subjectAliged, queryAligned, queryShadow;
    int subjectLength, queryLength, sqLength;
    double subjectAligedCenter[3], queryAlignedCenter[3], rotation[3][3];
    int numAliged;
};

#endif // TRANSFORM_H
