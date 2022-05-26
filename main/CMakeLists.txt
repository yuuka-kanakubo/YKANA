set(
    CMAKE_RUNTIME_OUTPUT_DIRECTORY
    ${CMAKE_HOME_DIRECTORY}
   )
add_executable(YKANA main.cpp)
target_link_libraries(YKANA YKANA_LIB)

