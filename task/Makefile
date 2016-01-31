BUILD_DIR = $(realpath ../../work)
SOURCE_FILES = utils.cpp mpi_trace.cpp mpi_scope.cpp matrix.cpp main.cpp

CXX_FLAGS = $(SOURCE_FILES) -Wall -fopenmp -lm

mpi:
		@mkdir $(BUILD_DIR)/mpi
		mpicxx -D COMPUTATION_MPI $(CXX_FLAGS) -o $(BUILD_DIR)/ppch_mpi

mpi_omp:
		@mkdir $(BUILD_DIR)/mpi_omp
		mpicxx -D COMPUTATION_MPI_OMP $(CXX_FLAGS) -o $(BUILD_DIR)/ppch_mpi_omp

plain:
		@mkdir $(BUILD_DIR)/plain
		mpicxx -D COMPUTATION_PLAIN $(CXX_FLAGS) -o $(BUILD_DIR)/ppch_plain

mpi-tv:
		@mkdir $(BUILD_DIR)/mpi-tv
		mpicxx-TV -D"MPITV_ABORT_AFTER=100" -D COMPUTATION_MPI -D MPITV_ENABLED $(CXX_FLAGS) -o $(BUILD_DIR)/ppch_mpi-tv

all: mpi mpi_omp plain

all-tv: all mpi-tv


MATRIX_SIZES = 100 200 300 400 500
mpi-size-test:
		$(foreach size,$(MATRIX_SIZES),mpicxx -D"MATRIX_SIZE=$(size)" -DCOMPUTATION_MPI $(CXX_FLAGS) -o $(BUILD_DIR)/ppch_mpi_$(size);)

clean-res:
		rm -f $(BUILD_DIR)/mpi/*
		rm -f $(BUILD_DIR)/mpi_omp/*
		rm -f $(BUILD_DIR)/plain/*

clean:
		@rm -rf $(BUILD_DIR)/*

.PHONY: all
