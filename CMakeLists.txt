# ---------- Project Setup ----------

cmake_minimum_required(VERSION 3.21)
 
project(VFEDemo
	LANGUAGES CXX
	HOMEPAGE_URL https://vc.tf.fau.de
	)

set(CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)

# ---------- Dependencies ----------

# ----- VTK -----
find_package(VTK REQUIRED)
include_directories(SYSTEM ${VTK_INCLUDE_DIRS})

# ----- Eigen -----
include(FetchContent)
FetchContent_Declare(eigen3
        GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
        GIT_TAG 9441d94dccccd5db8d64179516fdc5b53994a047
        )
FetchContent_Populate(eigen3)
add_library(eigen INTERFACE)
add_library (eigen3::eigen ALIAS eigen)
target_include_directories(eigen INTERFACE ${eigen3_SOURCE_DIR})

# ----- OpenMP -----

find_package(OpenMP)
if(OPENMP_FOUND)
  message(STATUS "Using OpenMP parallelization")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
else()
  message(STATUS "Not using OpenMP parallelization")
endif()

# ---------- Subdirectories ----------

add_subdirectory(variational)

# ---------- Executables ----------

# ----- Ocean -----
add_executable(Ocean vtkCriticalLinesExtractor.cxx vtkCriticalLinesExtractor.h Ocean.cxx)
target_link_libraries(Ocean ${VTK_LIBRARIES} eigen variational)
vtk_module_autoinit(TARGETS Ocean MODULES ${VTK_LIBRARIES})

# ----- EarthMantle -----
add_executable(EarthMantle vtkIsoLineExtractor.cxx vtkIsoLineExtractor.h EarthMantle.cxx)
target_link_libraries(EarthMantle ${VTK_LIBRARIES} eigen variational)
vtk_module_autoinit(TARGETS EarthMantle MODULES ${VTK_LIBRARIES})