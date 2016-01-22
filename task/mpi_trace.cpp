#include <mpi.h>
#include <cstddef>

#include "mpi_scope.h"
#include "mpi_trace.h"

const char* MPI_Trace::m_trace_file_name = "mpi_trace.trc";

MPI_Trace::MPI_Trace(unsigned int color)
    :
    m_color(color)
{
    MPI_Trace_Event(color, "entry");
}

MPI_Trace::~MPI_Trace()
{
    MPI_Trace_Event(m_color, "exit");
}

void MPI_Trace::Init()
{
    MPI_Trace_Init(m_trace_file_name);
}
