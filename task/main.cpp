#include <mpi.h>

#include <iostream>
#include <iomanip>
#include <cstdio>
#include <ctime>

#include "mpi_scope.h"
#include "matrix.h"

static const int msize = 1000;
static const double precision = 1e-6;

int main()
{
    {
#if defined(COMPUTATION_MPI) || defined(COMPUTATION_MPI_OMP)
        MPI_Scope scope; (void)scope;


        int world_rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
        //int world_rank = MPI::COMM_WORLD.Get_rank();
        if (world_rank == 0)
#endif //MPI
        {
            Matrix m(msize, msize);
            double d = msize;
            for (Matrix::index_type i = 0; i < m.getRowCount(); ++i)
            {
                m.at(i, i) = i + 1 + d;
                for (Matrix::index_type j = i + 1; j < m.getRowCount(); ++j)
                {
                    m.at(i, j) = i + 1 + d;
                    m.at(j, i) = i + 1 + d;
                }
            }
            std::cout << "Dimentions: " << m.getRowCount() << " x " << m.getColCount() << std::endl;
            //std::cout << m << std::endl;
            std::cout << "Computation start" << std::endl;
            double before, after;
#if defined(COMPUTATION_MPI) || defined(COMPUTATION_MPI_OMP)
            before = MPI_Wtime();
            //before = MPI::Wtime();
#elif defined(COMPUTATION_PLAIN)
            before = clock();
#endif
            m.compute_eigenvalues(precision);

#if defined(COMPUTATION_MPI) || defined(COMPUTATION_MPI_OMP)
            after = MPI_Wtime();
            //after = MPI::Wtime();
#elif defined(COMPUTATION_PLAIN)
            after = clock();
#endif
            /*
            std::cout << "Computed eigenvalues:" << std::endl;
            const std::vector<double>& evs = m.getEigenValues();
            for (size_t i = 0; i < evs.size(); ++i)
                std::cout << std::setw(10) << evs[i];
            std::cout << std::endl;
            */

            std::cout << std::endl;
            std::cout << "Time elapsed:";
#if defined(COMPUTATION_MPI) || defined(COMPUTATION_MPI_OMP)
            std::cout << std::setw(10) << after - before << "s" << std::endl;
#elif defined(COMPUTATION_PLAIN)
            std::cout << std::setw(10) << (after - before) / CLOCKS_PER_SEC << "s" << std::endl;
#endif
        }
#if defined(COMPUTATION_MPI) || defined(COMPUTATION_MPI_OMP)
        else
        {
            Matrix::find_max_mpi_server(msize, msize);
        }
#endif //MPI

        //getchar();
    }

    return 0;
}
