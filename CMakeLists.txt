cmake_minimum_required(VERSION 3.10.0)

set(CMAKE_C_COMPILER icc)
set(CMAKE_CXX_COMPILER icpc)
set(CMAKE_Fortran_COMPILER ifort)
project(main0 LANGUAGES CXX Fortran)
set(CMAKE_CXX_STANDARD 17)

# build lib

file(REMOVE_RECURSE
        ${PROJECT_SOURCE_DIR}/run
        ${PROJECT_SOURCE_DIR}/build)
        # ${PROJECT_SOURCE_DIR}/run/restart)

file(GLOB_RECURSE SRCS src/*.f90)
add_library(integral SHARED src/ellint.cpp)
add_executable(${PROJECT_NAME} ${SRCS})
target_link_libraries(main0 PUBLIC integral)
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/run)
