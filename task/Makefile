BUILD_DIR = $(realpath ../../work)
SOURCE_FILES = mpi_trace.cpp mpi_scope.cpp matrix.cpp main.cpp

TRACE_LIB = mpi-trace
TRACE_LIB_TV = mpi-trace-tv

CXX_FLAGS = $(SOURCE_FILES) -Wl,-rpath=$(BUILD_DIR) -L$(BUILD_DIR) -Wall -fopenmp -lm

CXX_SO_FLAGS = -Wall -Werror -fPIC -shared
C_A_FLAGS = -Wall -Werror

CXX_FLAGS_TRACE = $(CXX_FLAGS) -l$(TRACE_LIB)
CXX_FLAGS_TRACE_TV = $(CXX_FLAGS) -l$(TRACE_LIB_TV)

trace_lib:
		mpicc -c mpi_trace.c $(C_A_FLAGS) -o $(BUILD_DIR)/$(TRACE_LIB).o
		ar rsc $(BUILD_DIR)/lib$(TRACE_LIB).a $(BUILD_DIR)/$(TRACE_LIB).o

trace_so:
		mpicc -c mpi_trace.c $(CXX_SO_FLAGS) -o $(BUILD_DIR)/lib$(TRACE_LIB).so

trace_lib-tv:
		mpicc-TV -D MPITV_ENABLED -c mpi_trace.c $(C_A_FLAGS) -o $(BUILD_DIR)/$(TRACE_LIB_TV).o
		ar rsc $(BUILD_DIR)/lib$(TRACE_LIB_TV).a $(BUILD_DIR)/$(TRACE_LIB_TV).o

trace_so-tv:
		mpicc-TV -D MPITV_ENABLED -c mpi_trace.c $(CXX_SO_FLAGS) -o $(BUILD_DIR)/lib$(TRACE_LIB_TV).so

mpi: trace_lib
		mpicxx -D COMPUTATION_MPI $(CXX_FLAGS_TRACE) -o $(BUILD_DIR)/ppch_mpi

mpi_omp: trace_lib
		mpicxx -D COMPUTATION_MPI_OMP $(CXX_FLAGS_TRACE) -o $(BUILD_DIR)/ppch_mpi_omp

plain: trace_lib
		mpicxx -D COMPUTATION_PLAIN $(CXX_FLAGS_TRACE) -o $(BUILD_DIR)/ppch_plain

mpi-so-tv: trace_so-tv
		mpicxx-TV -D COMPUTATION_MPI $(CXX_FLAGS_TRACE_TV) -o $(BUILD_DIR)/ppch_mpi-tv

mpi-lib-tv: trace_lib-tv
		mpicxx-TV -D COMPUTATION_MPI $(CXX_FLAGS_TRACE_TV) -o $(BUILD_DIR)/ppch_mpi-tv

mpi-tv:
		mpicxx-TV -D COMPUTATION_MPI -D MPITV_ENABLED mpi_trace.c $(CXX_FLAGS) -o $(BUILD_DIR)/ppch_mpi-tv

all: mpi mpi_omp plain

all-tv: all mpi-tv

clean:
		@rm -rf $(BUILD_DIR)/*

.PHONY: all
