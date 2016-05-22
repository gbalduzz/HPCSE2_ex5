cmake_minimum_required(VERSION 3.1)
project(ex5)

set(CMAKE_CXX_FLAGS_DEBUG "-g3 -Wall -DDEBUG")
set(CMAKE_CXX_FLAGS_RELEASE "-DNDEBUG -O4 -fstrict-aliasing -ffast-math -mtune=native -march=native -funroll-loops ")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14  -Wstrict-aliasing=2 -fopenmp")

add_definitions(-DPROJECT_SOURCE_DIR=${PROJECT_SOURCE_DIR})
add_definitions(-D_GLIBCXX_PARALLEL)

include_directories(ex5 PRIVATE ${PROJECT_SOURCE_DIR}/include)
add_executable(ex5 main.cpp)

#include(CTest)
#add_executable(test_radius test_radius.cpp)
#add_test(radius test_radius)
add_executable(test_binomial tests/test_binomial.cpp)
add_executable(test_e2e tests/test_e2e.cpp)


