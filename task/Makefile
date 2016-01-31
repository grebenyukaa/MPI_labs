BUILD_DIR = $(realpath ../../work)
SOURCE_FILES = utils.cpp mpi_trace.cpp mpi_scope.cpp matrix.cpp main.cpp

CXX_FLAGS = $(SOURCE_FILES) -Wall -fopenmp -lm

MATRIX_SIZE = 1000

mpi:
		@mkdir $(BUILD_DIR)/mpi
		mpicxx -D"MATRIX_SIZE=$(MATRIX_SIZE)" -DCOMPUTATION_MPI $(CXX_FLAGS) -o $(BUILD_DIR)/ppch_mpi

mpi_omp:
		@mkdir $(BUILD_DIR)/mpi_omp
		mpicxx -D"MATRIX_SIZE=$(MATRIX_SIZE)" -DCOMPUTATION_MPI_OMP $(CXX_FLAGS) -o $(BUILD_DIR)/ppch_mpi_omp

plain:
		@mkdir $(BUILD_DIR)/plain
		mpicxx -D"MATRIX_SIZE=$(MATRIX_SIZE)" -DCOMPUTATION_PLAIN $(CXX_FLAGS) -o $(BUILD_DIR)/ppch_plain

mpi-tv:
		@mkdir $(BUILD_DIR)/mpi-tv
		mpicxx-TV -D"MATRIX_SIZE=$(MATRIX_SIZE)" -D"MPITV_ABORT_AFTER=100" -DCOMPUTATION_MPI -D MPITV_ENABLED $(CXX_FLAGS) -o $(BUILD_DIR)/ppch_mpi-tv

all: mpi mpi_omp plain

all-tv: all mpi-tv


MATRIX_SIZES = 500 600 700 800 900 1000
test:
		@mkdir $(BUILD_DIR)/mpi/by_size
		@mkdir $(BUILD_DIR)/mpi/by_node
		@mkdir $(BUILD_DIR)/mpi_omp/by_size
		@mkdir $(BUILD_DIR)/mpi_omp/by_node
		$(foreach size,$(MATRIX_SIZES),mpicxx -D"MATRIX_SIZE=$(size)" -DCOMPUTATION_MPI $(CXX_FLAGS) -o $(BUILD_DIR)/ppch_mpi_$(size);)
		$(foreach size,$(MATRIX_SIZES),mpicxx -D"MATRIX_SIZE=$(size)" -DCOMPUTATION_MPI_OMP $(CXX_FLAGS) -o $(BUILD_DIR)/ppch_mpi_omp_$(size);)

clean-test:
		@rm -f $(BUILD_DIR)/mpi/*
		@rm -f $(BUILD_DIR)/mpi_omp/*
		@rm -f $(BUILD_DIR)/plain/*

clean:
		@rm -rf $(BUILD_DIR)/*

.PHONY: all
