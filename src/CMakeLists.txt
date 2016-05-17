cmake_minimum_required(VERSION 2.6)


project(LevelSet_IMM C CXX)

set(CMAKE_BUILD_DIR "build")

set(MY_SOURCES  allocator.cpp ggLS.cpp mymath.cpp SCORE_Random.cpp grainhdl.cpp box.cpp marchingSquares.cpp Settings.cpp contourSector.cpp grahamScan.cpp junction.cpp minimalisticBoundary.cpp grainBoundary.cpp misorientation.cpp Quaternion.cpp IterativeGrainScheduler.cpp SquaresGrainScheduler.cpp)

#set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}")






set(CMAKE_CXX_FLAGS "-g -O3 -Wall -m64")
#set(CMAKE_CXX_FLAGS "-g -O3 -Wall -m64 -mkl -DMKL_ILP64")

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}") 

set(VORO++_FOUND true)
set(VORO++_INCLUDE_DIRS ${CMAKE_SOURCE_DIR}/voro++/)
set(VORO++_LIBRARIES ${CMAKE_SOURCE_DIR}/voro++/lib/libvoro++.a)
set(VORO++_DEFINITIONS " ")


# OPEN MP
FIND_PACKAGE(OpenMP REQUIRED)
if(OPENMP_FOUND)
message("OPENMP FOUND")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()


#find_package(Voro++ REQUIRED)


#INCLUDE(FindBoost.cmake)
#find_package(Boost COMPONENTS serialization mpi REQUIRED)
#find_package(Boost COMPONENTS serialization mpi REQUIRED)

#include_directories(${Boost_INCLUDE_PATH})
include_directories(${VORO++_INCLUDE_DIRS})
#include_directories(${MPI_INCLUDE_PATH})

link_directories(${Boost_LIBRARY_DIR} ${CMAKE_SOURCE_DIR})

#find_package(MPI REQUIRED)

ADD_DEFINITIONS("-std=c++0x")
ADD_DEFINITIONS("-DUSE_FFTW")
#ADD_DEFINITIONS("-DUSE_MKL")


add_executable(LevelSet_IMM ${MY_SOURCES})

target_link_libraries(LevelSet_IMM ${MPI_LIBRARIES}  ${VORO++_LIBRARIES} 
${Boost_LIBRARIES} -lfftw3 -lfftw3f -lm -lnuma -ljemalloc)


ADD_CUSTOM_TARGET (clearAll
	COMMAND ${CMAKE_COMMAND} -E remove *.png *.gnu *.ori Level* *.txt *.gif
	COMMAND ${CMAKE_COMMAND} -E echo "rm *.png *.gnu *.ori *.gif PARA* NET*"
)

ADD_CUSTOM_TARGET (rm_cmake
	COMMAND ${CMAKE_COMMAND} -E remove CMakeFiles/
)

ADD_CUSTOM_TARGET (doc
	COMMAND doxygen
)
