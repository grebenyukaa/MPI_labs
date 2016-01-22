#include "mpi_scope.h"
#include <mpi.h>

MPI_Scope::MPI_Scope()
{
    MPI::Init();
    MPI::Pcontrol(MPI_Trace::TRACEFILES, NULL, MPI_Trace::m_trace_file_name, 0);
    MPI::Pcontrol(MPI_Trace::TRARELEVEL, 1, 1, 1);
    MPI::Pcontrol(MPI_Trace::TRACENODE, 1000000, 1, 1);
    //MPI::COMM_WORLD.Set_errhandler(MPI::ERRORS_THROW_EXCEPTIONS);
}

MPI_Scope::~MPI_Scope()
{
    MPI::COMM_WORLD.Abort(0);
}
