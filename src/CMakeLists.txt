# Source files for library
add_library(YKSTAT_LIB STATIC
	          ReadIn.cpp Fill.cpp Stat.cpp Write.cpp
	   )
target_include_directories(YKSTAT_LIB PUBLIC ${PROJECT_SOURCE_DIR}/include)
