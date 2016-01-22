#ifndef MPI_TRACE_C_H
#define MPI_TRACE_C_H

#ifdef __cplusplus
extern "C"{
    extern int MPI_TRACEEVENT;
    extern int MPI_TRACELEVEL;
    extern int MPI_TRACENODE;
    extern int MPI_TRACEFILES;

    extern void MPI_Trace_Init(const char* filename);
    extern void MPI_Trace_Event(const unsigned int color, const char* name);
}
#else
extern int MPI_TRACEEVENT;
extern int MPI_TRACELEVEL;
extern int MPI_TRACENODE;
extern int MPI_TRACEFILES;

extern void MPI_Trace_Init(const char* filename);
extern void MPI_Trace_Event(const unsigned int color, const char* name);
#endif

#endif //MPI_TRACE_C_H
