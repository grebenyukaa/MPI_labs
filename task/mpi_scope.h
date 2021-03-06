#ifndef MPI_SCOPE_H
#define MPI_SCOPE_H

#include <string>

class MPI_Scope
{
public:
    MPI_Scope();
    ~MPI_Scope();
};

//

class MPI_Trace
{
public:
    static const std::string m_trace_file_name;
    enum Colors
    {
        ClRecv = 0xFF0000,
        ClSend = 0x00FF00,
        ClProbe = 0x0000FF
    };

    MPI_Trace(unsigned int color);
    ~MPI_Trace();
    static void Init(const int world_rank);
private:
    unsigned int m_color;
};

#endif //MPI_SCOPE_H
