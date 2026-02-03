target_include_directories(UBKLib
  PUBLIC
    ${CMAKE_CURRENT_LIST_DIR}/include
)

add_subdirectory(utils)
add_subdirectory(field_models)
add_subdirectory(tracer)

target_sources(UBKLib
  PUBLIC
    FILE_SET public_cxx_modules TYPE CXX_MODULES
    BASE_DIRS ${CMAKE_CURRENT_LIST_DIR}
    FILES ubk.cppm 
)
