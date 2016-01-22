#include "mpi_scope.h"
#include <mpi.h>
#include <cstddef>

MPI_Scope::MPI_Scope()
{
    //MPI::Init();
    MPI_Init(NULL, NULL);
    //MPI::COMM_WORLD.Set_errhandler(MPI::ERRORS_THROW_EXCEPTIONS);
}

MPI_Scope::~MPI_Scope()
{
    MPI_Abort(MPI_COMM_WORLD, 0);
    //MPI::COMM_WORLD.Abort(0);
}
