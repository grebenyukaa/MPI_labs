#include <mpi.h>
#include <cstddef>
#include <iostream>

#include "mpi_scope.h"

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

const char* MPI_Trace::m_trace_file_name = "mpi_trace.trc";

MPI_Trace::MPI_Trace(unsigned int color)
    :
    m_color(color)
{
    std::cout << "MPI_TRACE entry" << std::endl;
    MPI_Pcontrol(MPI_TRACEEVENT, "entry", color, 0, NULL);
}

MPI_Trace::~MPI_Trace()
{
    MPI_Pcontrol(MPI_TRACEEVENT, "exit", m_color, 0, NULL);
}

void MPI_Trace::Init()
{
    std::cout << "MPI_TRACE init" << std::endl;
    MPI_Pcontrol(MPI_TRACEFILES, NULL, m_trace_file_name, 0);
    MPI_Pcontrol(MPI_TRACELEVEL, 1, 1, 1);
    MPI_Pcontrol(MPI_TRACENODE, 1000000, 1, 1);
    std::cout << "  done" << std::endl;
}
