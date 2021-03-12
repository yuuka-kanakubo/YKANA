cmake_minimum_required(VERSION 2.8)

project(YKSTAT, CXX)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -Ofast -Wall")


add_subdirectory(src)
