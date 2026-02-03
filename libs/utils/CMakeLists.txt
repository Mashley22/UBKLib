target_sources(UBKLib
  PUBLIC
    FILE_SET public_cxx_modules TYPE CXX_MODULES
    BASE_DIRS ${CMAKE_CURRENT_LIST_DIR}/..
    FILES check.cppm 
      utils.cppm
      math.cppm
)
