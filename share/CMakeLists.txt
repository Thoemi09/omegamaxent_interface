add_subdirectory(cmake)

if(NOT CMAKE_INSTALL_PREFIX STREQUAL TRIQS_ROOT)

  configure_file(omegamaxent_interface.modulefile.in omegamaxent_interface.modulefile @ONLY)
  configure_file(omegamaxent_interfacevars.sh.in omegamaxent_interfacevars.sh @ONLY)

  install(
    FILES
      ${CMAKE_CURRENT_BINARY_DIR}/omegamaxent_interface.modulefile
      ${CMAKE_CURRENT_BINARY_DIR}/omegamaxent_interfacevars.sh
    DESTINATION share
  )
 
  message(STATUS "***************************************************************")
  message(STATUS "* Custom install Location. Use:                                ")
  message(STATUS "*                                                              ")
  message(STATUS "*   source ${CMAKE_INSTALL_PREFIX}/share/omegamaxent_interfacevars.sh      ")
  message(STATUS "*                                                              ")
  message(STATUS "* to set up the environment variables                          ")
  message(STATUS "***************************************************************")

endif()
