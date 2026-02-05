# since the c++ modules stuff
cmake_minimum_required(VERSION 3.28)

project(UBK VERSION 0.0.1 LANGUAGES CXX)

set(MINIMUM_CXX_STANDARD 23)

if(NOT DEFINED CMAKE_CXX_STANDARD)
  message(STATUS "No c++ standard found, defaulting to ${MINIMUM_CXX_STANDARD}")
  set(CMAKE_CXX_STANDARD ${MINIMUM_CXX_STANDARD})
elseif(CMAKE_CXX_STANDARD LESS ${MINIMUM_CXX_STANDARD})
  message(FATAL_ERROR "c++ standard set too low: ${CMAKE_CXX_STANDARD}. Need ${MINIMUM_CXX_STANDARD} atleast!")
endif()

if(UBKLIB_UNIT_TEST)
  add_compile_definitions(UBKLIB_UNIT_TEST UBKLIB_CHECK_FAIL_THROWS)
endif()

if(UBKLIB_SAMPLES)
  add_compile_options(-O3 -march=native)
  set(CMAKE_INTERPROCEDURAL_OPTIMIZATION ON)
endif()

set(PEDANTIC_COMPILE_OPTIONS
  -Wall 
  -Wextra 
  -Wpedantic 
  -Werror 
  -Wshadow 
  -Wconversion 
  -Wsign-conversion 
  -Wcast-align 
  -Wcast-qual 
  -Wformat=2 
  -Wundef 
  -Wnull-dereference 
  -Wdouble-promotion 
  -Wimplicit-fallthrough 
  -Wmissing-declarations 
  -Wredundant-decls 
  -Wstrict-overflow=5 
  -Wswitch-enum 
  -Wunreachable-code
)

add_library(UBKLib)

add_subdirectory(libs)

if(UBKLIB_UNIT_TEST)
  enable_testing()
  include(CTest)

  add_subdirectory(extern/catch2)
  list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/extern/catch2/extras")
  include(Catch)

  add_subdirectory(validation)

endif()

if(UBKLIB_SAMPLES)
  add_subdirectory(samples)
endif()
