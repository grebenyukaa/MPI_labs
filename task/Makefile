BUILD_DIR = $(realpath ../../work)
SOURCE_FILES = mpi_trace.cpp mpi_scope.cpp matrix.cpp main.cpp

TRACE_LIB = mpi-trace
TRACE_LIB_TV = mpi-trace-tv

CXX_FLAGS = $(SOURCE_FILES) -Wl,-rpath=$(BUILD_DIR) -L$(BUILD_DIR) -Wall -fopenmp -lm

mpi: trace_lib
		mpicxx -D COMPUTATION_MPI $(CXX_FLAGS_TRACE) -o $(BUILD_DIR)/ppch_mpi

mpi_omp: trace_lib
		mpicxx -D COMPUTATION_MPI_OMP $(CXX_FLAGS_TRACE) -o $(BUILD_DIR)/ppch_mpi_omp

plain: trace_lib
		mpicxx -D COMPUTATION_PLAIN $(CXX_FLAGS_TRACE) -o $(BUILD_DIR)/ppch_plain

mpi-tv:
		mpicxx-TV -D COMPUTATION_MPI -D MPITV_ENABLED $(CXX_FLAGS) -o $(BUILD_DIR)/ppch_mpi-tv

all: mpi mpi_omp plain

all-tv: all mpi-tv

clean:
		@rm -rf $(BUILD_DIR)/*

.PHONY: all
