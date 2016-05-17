cmake_minimum_required(VERSION 3.1)
project(ex5)

set(CMAKE_CXX_FLAGS_DEBUG '-g3')
set(CMAKE_CXX_FLAGS_RELEASE "-O4 -fstrict-aliasing -ffast-math -mtune=native -march=native ")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14  -Wstrict-aliasing=2 ")

add_definitions(-DPROJECT_SOURCE_DIR=${PROJECT_SOURCE_DIR})

set(SOURCE_FILES main.cpp)
add_executable(ex5 ${SOURCE_FILES})
