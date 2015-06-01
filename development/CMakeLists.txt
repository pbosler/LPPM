CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

project (LPM)
set ( LPM_VERSION_MAJOR 1)
set ( LPM_VERSION_MINOR 0)

FIND_PACKAGE(MPI REQUIRED)

enable_language(Fortran)

# set default build type to release
if ( NOT CMAKE_BUILD_TYPE)
	set (CMAKE_BUILD_TYPE RELEASE CACHE STRING "Choose type of build. Options are : NONE, DEBUG, RELEASE." FORCE)
endif ( NOT CMAKE_BUILD_TYPE)

# default installation
get_filename_component (default_prefix "." ABSOLUTE)
set (CMAKE_INSTALL_PREFIX ${default_prefix}/install CACHE STRING
      "Choose the installation directory; by default it installs in the install directory."
      FORCE)
      
# FFLAGS depend on the compiler
if ( ${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")
  # gfortran
  set (CMAKE_Fortran_FLAGS_RELEASE "-O2 -fopenmp")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g -fopenmp")
elseif ( ${CMAKE_Fortran_COMPILER_ID} STREQUAL "Intel")
  # ifort
  set (CMAKE_Fortran_FLAGS_RELEASE "-O2 -openmp")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g -openmp -debug all -debug-parameters -warn all -traceback")
else ( ${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")
  message ("CMAKE_Fortran_COMPILER full path: " ${CMAKE_Fortran_COMPILER})
  message ("Fortran compiler: " ${Fortran_COMPILER_NAME})
  message ("No optimized Fortran compiler flags are known, we just try -O2...")
  set (CMAKE_Fortran_FLAGS_RELEASE "-O2")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g")
endif ( ${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")

add_library( lpmFortran NumberKinds.f90 OutputWriter.f90 Logger.f90 IntegerList.f90 STDIntVector.f90
						SphereGeometry.f90 PlaneGeometry.f90 
						Particles.f90 Field.f90
						Edges.f90 Faces.f90 PolyMesh2d.f90 )

configure_file(cubedSphereSeed.dat ${CMAKE_BINARY_DIR}/cubedSphereSeed.dat COPYONLY)
configure_file(icosTriSeed.dat ${CMAKE_BINARY_DIR}/icosTriSeed.dat COPYONLY)
configure_file(quadRectSeed.dat ${CMAKE_BINARY_DIR}/quadRectSeed.dat COPYONLY)
configure_file(triHexSeed.dat ${CMAKE_BINARY_DIR}/triHexSeed.dat COPYONLY)

# UNIT TEST EXECUTABLES
ADD_EXECUTABLE(particlesAndFieldTest.exe ParticlesAndFieldTest.f90 )
TARGET_LINK_LIBRARIES(particlesAndFieldTest.exe lpmFortran)

ADD_EXECUTABLE(planarMeshTest.exe PlanarMeshTest.f90)
TARGET_LINK_LIBRARIES(planarMeshTest.exe lpmFortran)


