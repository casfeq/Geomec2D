cmake_minimum_required(VERSION 3.10)

if($ENV{DESKTOP_SESSION} STREQUAL "gnome")
    set(CMAKE_C_COMPILER mpicc)
    set(CMAKE_CXX_COMPILER mpicxx)
endif()

project(FVM-2D)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/CMakeModules/)

if($ENV{DESKTOP_SESSION} STREQUAL "ubuntu")
    set(PETSC_VERSION "3.12.2")
    set(PETSC_DIR "/home/carlos/Libraries/petsc-3.12.2")
    set(PETSC_ARCH "arch-linux2-c-opt")
elseif($ENV{DESKTOP_SESSION} STREQUAL "gnome")
    set(PETSC_VERSION "3.12.0")
    set(PETSC_DIR "/home/carlos/Libraries/petsc-3.12.0/release/shared")
    set(PETSC_ARCH "")
endif()
find_package(PETSc ${PETSC_VERSION} REQUIRED)
if(PETSC_FOUND)
	include_directories(${PETSC_INCLUDES})
endif()

include_directories(${PROJECT_SOURCE_DIR}/include/)
add_executable($ENV{sourceName} ${PROJECT_SOURCE_DIR}/source/$ENV{sourceName}.cpp)
target_link_libraries($ENV{sourceName} ${PETSC_LIBRARIES})