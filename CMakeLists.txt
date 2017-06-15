cmake_minimum_required(VERSION 3.1)

set(PROJECT_NAME_STR multicriterial_sample)

project (${PROJECT_NAME_STR})

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src)

file(GLOB SRC_FILES ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp
                    ${CMAKE_CURRENT_SOURCE_DIR}/src/*.hpp)

add_executable(${PROJECT_NAME_STR} ${SRC_FILES})

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
  set_target_properties(${PROJECT_NAME_STR} PROPERTIES LINK_FLAGS "-fopenmp")
  target_compile_options(${PROJECT_NAME_STR} PRIVATE -std=c++11 -Wall -fopenmp -march=native)
elseif (CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
  add_definitions(/EHsc /openmp)
endif()
