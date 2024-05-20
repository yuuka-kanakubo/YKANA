# Source files for library
add_library(YKANA_LIB STATIC
	          ReadIn.cpp Fill.cpp Stat.cpp Write.cpp BSTR.cpp
	   )
target_include_directories(YKANA_LIB PUBLIC ${PROJECT_SOURCE_DIR}/include)
