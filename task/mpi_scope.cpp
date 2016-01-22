#include "mpi_scope.h"
#include <mpi.h>

MPI_Scope::MPI_Scope()
{
    MPI::Init();
    //MPI::COMM_WORLD.Set_errhandler(MPI::ERRORS_THROW_EXCEPTIONS);
}

MPI_Scope::~MPI_Scope()
{
    MPI::COMM_WORLD.Abort(0);
}
