# List of all tests
set(all_tests test_Freq test_Time test_matrix test_block test_Freq2 test_Freq_boson test_Freq_chi test_Time_boson test_Time_chi)

foreach(test ${all_tests})
  get_filename_component(test_name ${test} NAME_WE)
  get_filename_component(test_dir ${test} DIRECTORY)
  add_test(NAME Py_${test_name} COMMAND ${TRIQS_PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/${test_dir}/${test_name}.py WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${test_dir})
  set_property(TEST Py_${test_name} APPEND PROPERTY ENVIRONMENT
    PATH=${PROJECT_BINARY_DIR}/omegamaxent/cpp:$ENV{PATH}
    PYTHONPATH=${PROJECT_BINARY_DIR}/python:$ENV{PYTHONPATH}
  )
endforeach()
