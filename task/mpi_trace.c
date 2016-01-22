#include <mpi.h>
#include <stddef.h>

#include "mpi_trace.h"

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

void MPI_Trace_Init(const char* filename)
{
    MPI_Pcontrol(MPI_TRACEFILES, NULL, filename, 0);
    MPI_Pcontrol(MPI_TRACELEVEL, 1, 1, 1);
    MPI_Pcontrol(MPI_TRACENODE, 1000000, 1, 1);
}

void MPI_Trace_Event(const unsigned int color, const char* name)
{
    MPI_Pcontrol(MPI_TRACEEVENT, name, color, 0, NULL);
}
