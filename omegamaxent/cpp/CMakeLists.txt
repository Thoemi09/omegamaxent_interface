
add_executable(OmegaMaxEnt graph_2D.cpp graph_3D.cpp generique.cpp OmegaMaxEnt_data.cpp OmegaMaxEnt_main.cpp)

target_compile_options(OmegaMaxEnt PRIVATE -std=c++11 -O0)

target_compile_definitions(OmegaMaxEnt PRIVATE ARMA_DONT_USE_WRAPPER)

find_package(GSL REQUIRED)
target_link_libraries(OmegaMaxEnt GSL::gsl GSL::gslcblas)
target_link_libraries(OmegaMaxEnt nda::blas_lapack)
target_link_libraries(OmegaMaxEnt triqs::fftw)

install(TARGETS OmegaMaxEnt DESTINATION bin)
