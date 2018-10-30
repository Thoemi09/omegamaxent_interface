set(all_tests test_Freq test_Time test_matrix test_block test_Freq2)

foreach(t ${all_tests})
  add_test(NAME ${t} COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/${t}.py)
endforeach()
