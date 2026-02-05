function(add_unit_test TEST_NAME)

  add_executable(${TEST_NAME})

  target_link_libraries(${TEST_NAME}
    PRIVATE
      UBKLib
  )

  target_compile_options(${TEST_NAME}
    PRIVATE
      ${PEDANTIC_COMPILE_OPTIONS}
  )

endfunction()

add_unit_test(sample_basic)

target_sources(sample_basic
  PRIVATE
    sample_basic.cpp
)
