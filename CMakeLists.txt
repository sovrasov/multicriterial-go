cmake_minimum_required(VERSION 3.1)

option(ENABLE_GCOV_SUPPORT "Turns on gcov-compatible compilation mode" OFF)

set(PROJECT_NAME_STR multicriterial_sample)
project(${PROJECT_NAME_STR})

#initial setup
enable_language(C)
enable_language(CXX)
set(CMAKE_CXX_STANDARD 11)
if((${ENABLE_GCOV_SUPPORT}) AND (CMAKE_CXX_COMPILER_ID MATCHES "GNU"))
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fprofile-arcs -ftest-coverage")
  message("Gcov support is enabled.")
endif()

#build gtest
set(GOOGLETEST_ROOT 3rd-party/googletest/googletest CACHE STRING "Google Test source root")
include_directories(
    ${PROJECT_SOURCE_DIR}/${GOOGLETEST_ROOT}
    ${PROJECT_SOURCE_DIR}/${GOOGLETEST_ROOT}/include
    )
set(GOOGLETEST_SOURCES
    ${PROJECT_SOURCE_DIR}/${GOOGLETEST_ROOT}/src/gtest-all.cc
    ${PROJECT_SOURCE_DIR}/${GOOGLETEST_ROOT}/src/gtest_main.cc
    )
foreach(_source ${GOOGLETEST_SOURCES})
    set_source_files_properties(${_source} PROPERTIES GENERATED 1)
endforeach()

add_library(googletest ${GOOGLETEST_SOURCES})

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
  target_link_libraries(googletest pthread)
endif()

#build optimizer
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/3rd-party)

file(GLOB SRC_FILES ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/*.hpp)

add_library(multi_objective_optimizer ${SRC_FILES})
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
  target_compile_options(multi_objective_optimizer INTERFACE -Wall -fopenmp -march=native)
elseif (CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
  add_definitions(/EHsc /openmp)
endif()

#build sample
add_executable(${PROJECT_NAME_STR} ${PROJECT_SOURCE_DIR}/sample/main.cpp)
target_link_libraries(${PROJECT_NAME_STR} multi_objective_optimizer)

#build tests
add_executable(unit_tests test/main.cpp test/solver_tests.cpp)
target_link_libraries(unit_tests
  multi_objective_optimizer googletest)

include(CTest)
enable_testing()
add_test(unit ${PROJECT_BINARY_DIR}/unit_tests)
