cmake_minimum_required (VERSION 3.5.1)

# Set architecture to arm64
if(APPLE)
    set(CMAKE_OSX_ARCHITECTURES "arm64")
    set(OpenMP_CXX_FLAGS "-Xpreprocessor -fopenmp -I/opt/homebrew/opt/libomp/include")
    set(OpenMP_CXX_LIB_NAMES "omp")
    set(OpenMP_omp_LIBRARY "/opt/homebrew/opt/libomp/lib/libomp.dylib")
endif()

project(demo CXX)

# Find OpenMP
find_package(OpenMP REQUIRED)
if(OpenMP_CXX_FOUND)
    message(STATUS "Found OpenMP ${OpenMP_CXX_VERSION}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

set(CMAKE_CXX_STANDARD 17)
option( BUILD_STATIC "Set to ON to include static versions of the library" OFF)

if(NOT OpenFHE_DIR)
    set(OpenFHE_DIR "" CACHE PATH "Path to OpenFHE build (leave empty if installed system-wide)")
endif()
find_package(OpenFHE CONFIG REQUIRED PATHS ${OpenFHE_DIR})
if (OpenFHE_FOUND)
    message(STATUS "FOUND PACKAGE OpenFHE")
    message(STATUS "OpenFHE Version: ${BASE_OPENFHE_VERSION}")
    message(STATUS "OpenFHE installed as shared libraries: ${OpenFHE_SHARED}")
    message(STATUS "OpenFHE include files location: ${OpenFHE_INCLUDE}")
    message(STATUS "OpenFHE lib files location: ${OpenFHE_LIBDIR}")
    message(STATUS "OpenFHE Native Backend size: ${OpenFHE_NATIVE_SIZE}")
else()
    message(FATAL_ERROR "PACKAGE OpenFHE NOT FOUND")
endif ()

set( CMAKE_CXX_FLAGS ${OpenFHE_CXX_FLAGS} )

include_directories( ${OPENMP_INCLUDES} )
include_directories( ${OpenFHE_INCLUDE} )
include_directories( ${OpenFHE_INCLUDE}/third-party/include )
include_directories( ${OpenFHE_INCLUDE}/core )
include_directories( ${OpenFHE_INCLUDE}/pke )
include_directories( ${OpenFHE_INCLUDE}/binfhe )
### add directories for other OpenFHE modules as needed for your project

link_directories( ${OpenFHE_LIBDIR} )
link_directories( ${OPENMP_LIBRARIES} )
if(BUILD_STATIC)
    set( CMAKE_EXE_LINKER_FLAGS "${OpenFHE_EXE_LINKER_FLAGS} -static")
    link_libraries( ${OpenFHE_STATIC_LIBRARIES} OpenMP::OpenMP_CXX )
else()
    set( CMAKE_EXE_LINKER_FLAGS ${OpenFHE_EXE_LINKER_FLAGS} )
    link_libraries( ${OpenFHE_SHARED_LIBRARIES} OpenMP::OpenMP_CXX )
endif()

### ADD YOUR EXECUTABLE(s) HERE
### add_executable( EXECUTABLE-NAME SOURCES )
###
### EXAMPLE:
### add_executable( test demo-simple-example.cpp )

add_executable( ranking
					src/utils-basics.cpp
					src/utils-eval.cpp
					src/utils-matrices.cpp
					src/utils-ptxt.cpp
					src/ranking.cpp
					src/test-ranking.cpp
				)


add_executable( minimum
					src/utils-basics.cpp
					src/utils-eval.cpp
					src/utils-matrices.cpp
					src/utils-ptxt.cpp
					src/ranking.cpp
					src/minimum.cpp
					src/test-minimum.cpp
				)


add_executable( median
					src/utils-basics.cpp
					src/utils-eval.cpp
					src/utils-matrices.cpp
					src/utils-ptxt.cpp
					src/ranking.cpp
					src/median.cpp
					src/test-median.cpp
				)


add_executable( sorting
					src/utils-basics.cpp
					src/utils-eval.cpp
					src/utils-matrices.cpp
					src/utils-ptxt.cpp
					src/sorting.cpp
					src/test-sorting.cpp
				)


add_executable( statistics-test
                                        src/utils-basics.cpp
                                        src/utils-eval.cpp
                                        src/utils-matrices.cpp
                                        src/utils-ptxt.cpp
                                        src/ranking.cpp
                                        src/statistics.cpp
                                        src/cohort-io.cpp
                                        src/test-statistics.cpp
                                )

add_executable( mpc-statistics
                                        src/utils-basics.cpp
                                        src/utils-eval.cpp
                                        src/utils-matrices.cpp
                                        src/utils-ptxt.cpp
                                        src/ranking.cpp
                                        src/ranking-mpc.cpp
                                        src/statistics.cpp
                                        src/cohort-io.cpp
                                        src/mpc-statistics.cpp
                                )
