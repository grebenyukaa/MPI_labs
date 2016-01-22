#include <mpi.h>
#include "mpi_scope.h"

const char* MPI_Trace::m_trace_file_name = "mpi_trace.trc";

MPI_Trace::MPI_Trace(unsigned int color)
{
    MPI::Pcontrol(MPI_Trace::TRACEEVENT, "entry", color, 0, NULL);
}

MPI_Trace::~MPI_Trace()
{
    MPI::Pcontrol(MPI_Trace::TRACEEVENT, "exit", m_color, 0, NULL);
}
