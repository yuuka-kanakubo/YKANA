cmake_minimum_required(VERSION 2.8)

project(YKANA, CXX)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -Ofast -Wall")
message(${CMAKE_CXX_FLAGS})

add_subdirectory(src)
add_subdirectory(main)
