function(add_sample SAMPLE_NAME)

  add_executable(${SAMPLE_NAME})

  target_link_libraries(${SAMPLE_NAME}
    PRIVATE
      UBKLib
  )

  target_compile_options(${SAMPLE_NAME}
    PRIVATE
      ${PEDANTIC_COMPILE_OPTIONS}
  )

endfunction()

add_sample(sample_basic)

target_sources(sample_basic
  PRIVATE
    sample_basic.cpp
)

add_sample(sample)

target_sources(sample
  PRIVATE
    sample.cpp
)

add_sample(sample_k0)

target_sources(sample_k0
  PRIVATE
    sample_k0.cpp
)

add_sample(sample_cli)

target_sources(sample_cli
  PRIVATE
    sample_cli.cpp
)

if(UBKLIB_PY_SAMPLES)

  pybind11_add_module(sample_py MODULE sample_py.cpp)

  target_link_libraries(sample_py
    PRIVATE
      UBKLib
  )

  target_compile_options(sample_py
    PRIVATE
      ${PEDANTIC_COMPILE_OPTIONS}
  )

  set_target_properties(sample_py PROPERTIES
    OUTPUT_NAME UBKLibpp
  )
endif()
