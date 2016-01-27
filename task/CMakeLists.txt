project(task)

find_package(MPI)
include_directories(SYSTEM ${MPI_INCLUDE_PATH})

SET(CMAKE_CXX_COMPILER "mpicxx")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Werror")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DCMAKE_BUILD_TYPE=Debug")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -DCMAKE_BUILD_TYPE=Release")

cmake_minimum_required(VERSION 2.8)
aux_source_directory(. SRC_LIST)
add_executable(${PROJECT_NAME} ${SRC_LIST})
target_link_libraries(${PROJECT_NAME} ${MPI_C_LIBRARIES})
