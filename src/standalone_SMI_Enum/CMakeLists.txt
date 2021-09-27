cmake_minimum_required (VERSION 3.5)


#set(D_GLIBCXX_USE_CXX11_ABI 0)

set(filename	
		SMI_Enumerator.cpp
	)

set(exename
        SMI_Enumerator
    )


project(PROG)

# ENV setting
set(RDBASE "/home/tom61212/share/rdkit_py/rdkit")
set(MYBOOST "/home/tom61212/share/boost_1_75_0_py/build")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${RDBASE}/Code/cmake/Modules")


# Compiler setting
# note that if you haven't installed/built the toolkit with CoordGen, you'll have problems with this.
set(CMAKE_CXX_COMPILER "/home/tom61212/bin/gcc920/bin/g++")
set(CMAKE_C_COMPILER "/home/tom61212/bin/gcc920/bin/gcc")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall  -std=c++11 -DRDK_BUILD_COORDGEN_SUPPORT=ON" )


# Boost settings
#set(MYBOOST "/home/tom61212/share/boost_1_75_0_py/build")
set(Boost_USE_STATIC_LIBS ON)
set(Boost_USE_MULTITHREADED OFF)
set(Boost_USE_STATIC_RUNTIME ON)
include_directories( ${MYBOOST}/include )
link_directories( ${MYBOOST}/lib )

# RDKit setting
# specify where CMake can find the RDKit libraries
include_directories ( ${RDBASE}/Code ) #${CAIRO_INCLUDE_DIRS}
link_directories ( ${RDBASE}/lib )
set(RDKit_LIBS  RDKitChemReactions  RDKitFileParsers    RDKitSmilesParse        RDKitDepictor
                RDKitRDGeometryLib  RDKitRDGeneral      RDKitSubstructMatch     RDKitSubgraphs
                RDKitMolDraw2D      RDKitGraphMol       RDKitDistGeometry       RDKitDistGeomHelpers
                RDKitMolAlign       RDKitOptimizer      RDKitForceField         RDKitForceFieldHelpers
                RDKitAlignment      RDKitForceField     RDKitMolTransforms      RDKitEigenSolvers
    )
find_package (Threads)
set(RDKit_THREAD_LIBS Threads::Threads)


#OBabel setting
#include_directories("/home/akitainu/bin/openbabel-install/include/openbabel-2.0")   # .../openbabel-2.0 or .../openbabel3
#link_directories( "/home/akitainu/bin/openbabel-install/lib/libopenbabel.so.5" )
link_libraries( "/home/tom61212/bin/obabel-install/lib/libopenbabel.so.7" )
include_directories("/home/tom61212/bin/obabel-install/include/openbabel3")


#Eigen3 setting
include_directories("/home/tom61212/share/Eigen3/eigen-3.3.7/build/include/eigen3")


#GCC 9.2.0 setting
include_directories("/home/tom61212/bin/gcc920/lib64")


set( LIBS 
		${RDKIT_LIBRARIES} 
		${RDKit_THREAD_LIBS}  
		#"/home/akitainu/bin/openbabel-install/lib/libopenbabel.so.5"
		"/home/tom61212/bin/obabel-install/lib/libopenbabel.so.7"
	) #${CAIRO_LIBRARIES} Boost::iostreams


set(EXECUTABLE_OUTPUT_PATH ${CMAKE_SOURCE_DIR})
add_executable(${exename} ${filename})
target_link_libraries(${exename} ${LIBS} ${RDKit_LIBS}) #"/home/akitainu/bin/openbabel-install/lib/libopenbabel.so.5"
#set_target_properties(${exename} PROPERTIES
#                                        CXX_STANDARD 11  # C++11...
#                                        CXX_STANDARD_REQUIRED ON  #...is required...
#                                        CXX_EXTENSIONS OFF  #...without compiler extensions like gnu++11
#                    )

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


