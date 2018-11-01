cmake_minimum_required (VERSION 2.6)

project(OmegaMaxEnt)

add_executable(OmegaMaxEnt graph_2D.cpp generique.cpp OmegaMaxEnt_data.cpp OmegaMaxEnt_main.cpp)

if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  set(LAPACK_ROOT /usr/local/opt/lapack CACHE PATH "lapack install path")
  if (NOT IS_DIRECTORY "${LAPACK_ROOT}/include")
    message(FATAL_ERROR "Please install the brew version of LAPACK first with command 'brew install lapack' or set LAPACK_ROOT\n")
  else()
	message("Darwin build settings")
	target_include_directories(OmegaMaxEnt PRIVATE "${LAPACK_ROOT}/include")
	link_directories("${LAPACK_ROOT}/lib")
  endif()
endif()

target_compile_options(OmegaMaxEnt PRIVATE -std=c++11 -O0)

target_compile_definitions(OmegaMaxEnt PRIVATE ARMA_DONT_USE_WRAPPER)

target_link_libraries(OmegaMaxEnt blas lapack fftw3 gsl)

install(TARGETS OmegaMaxEnt DESTINATION bin)

