configure_file(omegamaxent_interface-config.cmake.in omegamaxent_interface-config.cmake @ONLY)
configure_file(omegamaxent_interface-config-version.cmake.in omegamaxent_interface-config-version.cmake @ONLY)
install(
  FILES
    ${CMAKE_CURRENT_BINARY_DIR}/omegamaxent_interface-config.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/omegamaxent_interface-config-version.cmake
  DESTINATION lib/cmake/omegamaxent_interface
)
