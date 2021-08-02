set(
    CMAKE_RUNTIME_OUTPUT_DIRECTORY
    ${CMAKE_HOME_DIRECTORY}
   )
add_executable(YKSTAT main.cpp)
target_link_libraries(YKSTAT YKSTAT_LIB)

