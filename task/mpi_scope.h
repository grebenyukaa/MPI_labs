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

#ifdef MPITV_ENABLED
#include <pcontrol.h>
#define MPI_TRACEEVENT TRACEEVENT
#define MPI_TRACELEVEL TRACELEVEL
#define MPI_TRACENODE  TRACENODE
#define MPI_TRACEFILES TRACEFILES
#define MPI_TRACEFLUSH TRACEFLUSH
#else
#define MPI_TRACEEVENT 8
#define MPI_TRACELEVEL 6
#define MPI_TRACENODE  3
#define MPI_TRACEFILES 101
#define MPI_TRACEFLUSH 5
#endif

#ifdef _MPITV_ENABLED
#define MPI_TRACE_EVENT(x, color) \
MPI_Pcontrol(MPI_TRACEEVENT, "entry", color, 0, NULL);\
x;\
MPI_Pcontrol(MPI_TRACEEVENT, "exit", color, 0, NULL);
#else
#define MPI_TRACE_EVENT(x, color) x;
#endif

class MPI_Trace
{
public:
    static const std::string m_trace_file_name;
    enum Colors
    {
        ClRecv = 1,
        ClSend = 2
    };

    MPI_Trace(unsigned int color);
    ~MPI_Trace();
    static void Init(const int world_rank);
private:
    unsigned int m_color;
};

#endif //MPI_SCOPE_H
