set(all_tests test_Freq test_Time test_matrix test_block test_Freq2)

foreach(t ${all_tests})
  add_test(NAME ${t} COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/${t}.py)
endforeach()

# Set the PythonPath for all tests
set_property(TEST ${all_tests} PROPERTY ENVIRONMENT
  PATH=${CMAKE_BINARY_DIR}/omegamaxent/cpp:$ENV{PATH}
  PYTHONPATH=${CMAKE_SOURCE_DIR}/python:$ENV{PYTHONPATH}
)
