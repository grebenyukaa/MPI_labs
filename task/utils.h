#ifndef UTILS_H
#define UTILS_H

#include <fstream>

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

//AP - arithmetic progression
template<class T>
class APUtils
{
    static T sum(const T count, const T from = 0, const T step = 1)
    {
        return (2*from + step*(count - 1)) * count / 2;
    }
};

#endif //UTILS_H
