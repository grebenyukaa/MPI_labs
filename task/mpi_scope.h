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
int MPI_TRACEEVENT = TRACEEVENT;
int MPI_TRACELEVEL = TRACELEVEL;
int MPI_TRACENODE  = TRACENODE;
int MPI_TRACEFILES = TRACEFILES;
#else
int MPI_TRACEEVENT = 8;
int MPI_TRACELEVEL = 6;
int MPI_TRACENODE  = 3;
int MPI_TRACEFILES = 101;
#endif

#define MPI_TRACE_EVENT(x, color) \
    MPI_Pcontrol(MPI_TRACEEVENT, "entry", color, 0, NULL);\
    x;\
    MPI_Pcontrol(MPI_TRACEEVENT, "exit", color, 0, NULL);

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
