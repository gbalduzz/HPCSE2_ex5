cmake_minimum_required(VERSION 3.1)
project(ex5)

set(CMAKE_CXX_FLAGS_DEBUG "-g3 -Wall -DDEBUG")
set(CMAKE_CXX_FLAGS_RELEASE "-DNDEBUG -O4 -fstrict-aliasing -ffast-math -mtune=native -march=native -funroll-loops ")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14  -Wstrict-aliasing=2 -fopenmp")

add_definitions(-DPROJECT_SOURCE_DIR=${PROJECT_SOURCE_DIR})
add_definitions(-D_GLIBCXX_PARALLEL)

#compile ispc files
include(ispc.cmake)
include(m4generate.cmake)
m4generate(src/p2e.ispc.m4)
m4generate(src/e2p.ispc.m4)
set(ISPC_SOURCES src/p2p.ispc ${M4DIR}/p2e.ispc ${M4DIR}/e2p.ispc)
create_ispc_library("${ISPC_SOURCES}")

include_directories(ex5 PRIVATE ${PROJECT_SOURCE_DIR}/include)
add_executable(ex5 main.cpp)
target_link_libraries(ex5 ispc)




