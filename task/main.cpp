#include <iostream>
#include <iomanip>
#include <cstdio>

#include <mpi.h>

#include "mpi_scope.h"
#include "matrix.h"

static const int msize = 8;
static const double precision = 1e-7;

int main()
{
    {
#if defined(MATRIX_MUL_MPI) || defined(MATRIX_MUL_MPI_OMP)
        MPI_Scope scope; (void)scope;


        int world_rank = MPI::COMM_WORLD.Get_rank();
        if (world_rank == 0)
#endif //MPI
        {
            Matrix m(msize, msize);
            double d = msize;
            for (int i = 0; i < m.getRowCount(); ++i)
            {
                m.at(i, i) = i + 1 + d;
                for (int j = i + 1; j < m.getRowCount(); ++j)
                {
                    m.at(i, j) = i + 1 + d;
                    m.at(j, i) = i + 1 + d;
                }
            }
            std::cout << m << std::endl;

            m.compute_eigenvalues(precision);
            const std::vector<double>& evs = m.getEigenValues();
            for (size_t i = 0; i < evs.size(); ++i)
                std::cout << std::setw(10) << evs[i];
            std::cout << std::endl;
        }
#if defined(MATRIX_MUL_MPI) || defined(MATRIX_MUL_MPI_OMP)
        else
        {
            Matrix::find_max_mpi_server(msize, msize);
        }
#endif //MPI

        //getchar();
    }

    return 0;
}
