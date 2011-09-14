#ifndef DAR_H
#define DAR_H

class Dar
{
public:
    Dar(int n=0);
    ~Dar();

    void reset(int n);
    int len() const;

    double * operator[] (int i);
    const double * operator[] (int i) const;

    double (*ptr()) [3];
    const double (*ptr() const)[3];

private:
    double (*p)[3];
    int s;

    Dar(const Dar &);
    Dar operator= (const Dar &);
};
#endif // DAR_H
