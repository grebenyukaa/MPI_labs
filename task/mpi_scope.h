#ifndef MPI_SCOPE_H
#define MPI_SCOPE_H

#ifdef MPITV_ENABLED
#include<pcontrol.h>
#endif

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
    static const char* m_trace_file_name;
    enum Colors
    {
        ClRecv = 0xFF0000,
        ClSend = 0x00FF00,
        ClProbe = 0x0000FF
    };

    //checked on cluster -- was unable to find anything regarding mpicxx-TV
#ifndef MPITV_ENABLED
    enum PControlCodes
    {
        PC_TRACEEVENT = 8,
        PC_TRARELEVEL = 6,
        PC_TRACENODE  = 3,
        PC_TRACEFILES = 101
    };
#else
    enum PControlCodes
    {
        PC_TRACEEVENT = TRACEEVENT,
        PC_TRARELEVEL = TRACELEVEL,
        PC_TRACENODE  = TRACENODE,
        PC_TRACEFILES = TRACEFILES
    };
#endif

    MPI_Trace(unsigned int color);
    ~MPI_Trace();
    static void Init();
private:
    unsigned int m_color;
};

#endif //MPI_SCOPE_H
