#include <mpi.h>
#include <cstddef>

#include "mpi_scope.h"

#ifdef MPITV_ENABLED
#include<pcontrol.h>
#endif

const char* MPI_Trace::m_trace_file_name = "mpi_trace.trc";

MPI_Trace::MPI_Trace(unsigned int color)
{
    MPI::Pcontrol(MPI_Trace::PC_TRACEEVENT, "entry", color, 0, NULL);
}

MPI_Trace::~MPI_Trace()
{
    MPI::Pcontrol(MPI_Trace::PC_TRACEEVENT, "exit", m_color, 0, NULL);
}

void MPI_Trace::Init()
{
    MPI::Pcontrol(MPI_Trace::PC_TRACEFILES, NULL, m_trace_file_name, 0);
    MPI::Pcontrol(MPI_Trace::PC_TRARELEVEL, 1, 1, 1);
    MPI::Pcontrol(MPI_Trace::PC_TRACENODE, 1000000, 1, 1);
}
