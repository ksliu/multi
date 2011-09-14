#ifndef  READ_CONF_H
#define  READ_CONF_H
#include <string>
#include <map>

class Para
{
public:
    static Para * instance(const char *s);
    static Para * getPointer();
    static void destorySelf();

    double get_value(const std::string &key) const;
private:
    Para(const char *s);
    Para(const Para &);
    Para operator=(const Para &);

    std::map<std::string, double> dict;
    static Para * data;
};

#endif  /*READ_CONF_H*/
