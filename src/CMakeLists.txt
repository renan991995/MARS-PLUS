cmake_minimum_required (VERSION 2.6)


#set(CMAKE_CXX_COMPILER g++)
set(CMAKE_CXX_COMPILER "/home/tom61212/bin/gcc920/bin/g++")
set(CMAKE_C_COMPILER "/home/tom61212/bin/gcc920/bin/gcc")


set(filename	
		CASES_NEU.cpp		CASES_NEU.h
		CASES_IL.cpp		CASES_IL.h
		CASES_IL_INDEPENDENT.cpp		CASES_IL_INDEPENDENT.h
		UTILITY.cpp			UTILITY.h
		MOLECULE.cpp		MOLECULE.h
		IL.cpp				IL.h
		ELEMENTS.cpp		ELEMENTS.h
		PARAMETER.h
		main.cpp			
	)

project(GA)
add_executable(MARS ${filename})
#target_link_libraries(MARS "/home/akitainu/bin/openbabel-install/lib/libopenbabel.so.5")
#include_directories("/home/akitainu/bin/openbabel-install/include/openbabel-2.0")
target_link_libraries(MARS "/home/tom61212/bin/obabel-install/lib/libopenbabel.so.7")
include_directories("/home/tom61212/bin/obabel-install/include/openbabel3")
include_directories("/home/tom61212/share/Eigen3/eigen-3.3.7/build/include/eigen3")
#include_directories("/home/akitainu/bin/openbabel-install/lib")

message("Default build type is ${CMAKE_BUILD_TYPE}")

message("CMAKE_C_FLAGS_DEBUG is ${CMAKE_C_FLAGS_DEBUG}")
message("CMAKE_C_FLAGS_RELEASE is ${CMAKE_C_FLAGS_RELEASE}")
message("CMAKE_C_FLAGS_RELWITHDEBINFO is ${CMAKE_C_FLAGS_RELWITHDEBINFO}")
message("CMAKE_C_FLAGS_MINSIZEREL is ${CMAKE_C_FLAGS_MINSIZEREL}")

message("CMAKE_CXX_FLAGS_DEBUG is ${CMAKE_CXX_FLAGS_DEBUG}")
message("CMAKE_CXX_FLAGS_RELEASE is ${CMAKE_CXX_FLAGS_RELEASE}")
message("CMAKE_CXX_FLAGS_RELWITHDEBINFO is ${CMAKE_CXX_FLAGS_RELWITHDEBINFO}")
message("CMAKE_CXX_FLAGS_MINSIZEREL is ${CMAKE_CXX_FLAGS_MINSIZEREL}")

if (NOT CMAKE_BUILD_TYPE)
    message(STATUS "No build type selected, default to Release")
    set(CMAKE_BUILD_TYPE "debug")
endif()


