target_include_directories(UBKLib
  PUBLIC
    ${CMAKE_CURRENT_LIST_DIR}/include
)

add_subdirectory(utils)
add_subdirectory(field_models)
add_subdirectory(tracer)
