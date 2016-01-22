#include <mpi.h>
#include <cstddef>

#include "mpi_scope.h"

MPI_Scope::MPI_Scope()
{
    //MPI::Init();
    MPI_Init(NULL, NULL);

    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    MPI_Trace::Init(world_rank);
}

MPI_Scope::~MPI_Scope()
{
    MPI_Abort(MPI_COMM_WORLD, 0);
    //MPI::COMM_WORLD.Abort(0);
}
