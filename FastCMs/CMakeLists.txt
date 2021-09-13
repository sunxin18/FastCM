# CMake 最低版本号要求

cmake_minimum_required (VERSION 2.8)

# 项目信息

project (Demo)

# 指定生成目标
add_definitions(-std=c++11)
add_executable(Demo core.hpp graph.h graph_IO.h main.cpp)
