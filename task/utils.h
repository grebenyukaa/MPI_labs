#ifndef UTILS_H
#define UTILS_H

#include <fstream>
#include <assert.h>
#include "matrix.h"

class Logger
{
public:
    Logger(const Matrix::index_type nodeid);
    ~Logger();

    template<class T>
    std::ofstream& operator<<(const T& data)
    {
        *m_ofs << data;
        return *m_ofs;
    }
private:
    std::ofstream* m_ofs;
};

//

class APUtils
{
public:
    template<class T>
    static T sum(const T from, const T to, const T step = 1)
    {
        assert(to >= from);
        return ((to + from) * ((to - from) / step + 1) / step) / 2;
    }
};
#endif //UTILS_H
