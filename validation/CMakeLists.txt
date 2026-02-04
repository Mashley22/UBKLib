function(add_unit_test TEST_NAME)

  add_executable(${TEST_NAME})

  target_link_libraries(${TEST_NAME}
    PRIVATE
      UBKLib
      Catch2::Catch2WithMain
  )

  target_compile_options(${TEST_NAME}
    PRIVATE
      ${PEDANTIC_COMPILE_OPTIONS}
  )

  catch_discover_tests(${TEST_NAME})

endfunction()

add_unit_test(test_utils)

target_sources(test_utils
  PRIVATE
    test_utils.cpp
)

add_unit_test(test_tracer)

target_sources(test_tracer
  PRIVATE
    test_tracer.cpp
)
