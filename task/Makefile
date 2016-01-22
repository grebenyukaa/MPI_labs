mpi:
	mpicxx *.cpp *.h -D COMPUTATION_MPI -Wall -fopenmp -lm -o ../work/ppch_mpi

mpi_omp:
	mpicxx *.cpp *.h -D COMPUTATION_MPI_OMP -Wall -fopenmp -lm -o ../work/ppch_mpi_omp

plain:
	mpicxx *.cpp *.h -D COMPUTATION_PLAIN -Wall -fopenmp -lm -o ../work/ppch_plain

mpi-tv:
	mpicxx-TV *.cpp *.h -D COMPUTATION_MPI -Wall -fopenmp -lm -o ../work/ppch_mpi_tv

all:
	mpi mpi_omp plain mpi-tv

clean:
	@rm -rf ../work/*

.PHONY: all