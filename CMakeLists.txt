cmake_minimum_required (VERSION 2.6)

project(OmegaMaxEnt)

include_directories("armadillo-5.600.2/include" "/usr/local/opt/lapack/include")

add_definitions(-std=c++11 -O0 -DARMA_DONT_USE_WRAPPER)

link_directories("/usr/local/opt/lapack/lib")

add_executable(OmegaMaxEnt graph_2D.cpp generique.cpp OmegaMaxEnt_data.cpp OmegaMaxEnt_main.cpp)

target_link_libraries(OmegaMaxEnt blas lapack fftw3 gsl)
