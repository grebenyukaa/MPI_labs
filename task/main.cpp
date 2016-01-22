#include <iostream>
#include <iomanip>
#include <cstdio>
#include <ctime>

#include <mpi.h>

#include "mpi_scope.h"
#include "matrix.h"

static const int msize = 1000;
static const double precision = 1e-6;

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
            //std::cout << m << std::endl;
            std::cout << "Computation start" << std::endl;
            double before, after;
#if defined(MATRIX_MUL_MPI) || defined(MATRIX_MUL_MPI_OMP)
            before = MPI::Wtime();
#elif defined(MATRIX_MUL_PLAIN)
            before = clock();
#endif
            m.compute_eigenvalues(precision);

#if defined(MATRIX_MUL_MPI) || defined(MATRIX_MUL_MPI_OMP)
            after = MPI::Wtime();
#elif defined(MATRIX_MUL_PLAIN)
            after = clock();
#endif
            /*std::cout << "Computed eigenvalues:" << std::endl;
            const std::vector<double>& evs = m.getEigenValues();
            for (size_t i = 0; i < evs.size(); ++i)
                std::cout << std::setw(10) << evs[i];
            std::cout << std::endl;*/

            std::cout << std::endl;
            std::cout << "Time elapsed:";
#if defined(MATRIX_MUL_MPI) || defined(MATRIX_MUL_MPI_OMP)
            std::cout << std::setw(10) << after - before << "s" << std::endl;
#elif defined(MATRIX_MUL_PLAIN)
            std::cout << std::setw(10) << (after - before) / CLOCKS_PER_SEC << "s" << std::endl;
#endif
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
