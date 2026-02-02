target_sources(UBKLib
  PUBLIC
    FILE_SET public_cxx_modules TYPE CXX_MODULES
    FILES traits.cppm
      electric_potentials/cross_tail.cppm
      field_models.cppm
)
