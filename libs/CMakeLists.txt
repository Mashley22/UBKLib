target_include_directories(UBKLib
  PUBLIC
    ${CMAKE_CURRENT_LIST_DIR}/include
)

if(UBKLIB_UNIT_TEST OR UBKLIB_SAMPLES)
  target_compile_options(UBKLib
    PRIVATE
      ${PEDANTIC_COMPILE_OPTIONS}
  )
endif()

add_subdirectory(utils)
add_subdirectory(field_models)
add_subdirectory(seed_points)
add_subdirectory(tracer)

target_sources(UBKLib
  PUBLIC
    FILE_SET public_cxx_modules TYPE CXX_MODULES
    BASE_DIRS ${CMAKE_CURRENT_LIST_DIR}
    FILES ubk.cppm 
)
