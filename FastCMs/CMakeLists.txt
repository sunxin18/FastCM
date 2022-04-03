
cmake_minimum_required (VERSION 2.8)
SET(CMAKE_BUILD_TYPE DEBUG)

project (Demo)

add_definitions(-std=c++11)
add_executable(fastcm core.cpp graph.cpp main.cpp)
