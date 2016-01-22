#ifndef MPI_SCOPE_H
#define MPI_SCOPE_H

#define TRACE_FILE_NAME

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
    enum PControlCodes
    {
        TRACEEVENT = 8,
        TRARELEVEL = 6,
        TRACENODE  = 3,
        TRACEFILES = 101
    };

    MPI_Trace(unsigned int color);
    ~MPI_Trace();
private:
    unsigned int m_color;
};

#endif //MPI_SCOPE_H
