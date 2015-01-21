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
  set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g -openmp")
else ( ${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")
  message ("CMAKE_Fortran_COMPILER full path: " ${CMAKE_Fortran_COMPILER})
  message ("Fortran compiler: " ${Fortran_COMPILER_NAME})
  message ("No optimized Fortran compiler flags are known, we just try -O2...")
  set (CMAKE_Fortran_FLAGS_RELEASE "-O2")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g")
endif ( ${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")

add_library( base NumberKinds3.f90 OutputWriter2.f90 IntegerList2.f90 Logger2.f90 SphereGeom3.f90 PlaneGeometry.f90)

add_library( mesh Particles.f90 Edges.f90 Panels.f90 SphereMesh2.f90 PlaneMesh.f90)
target_link_libraries( mesh base )

add_library( interp ssrfpack.f stripack.f STRIPACKInterface2.f90 SSRFPACKInterface2.f90 bivar.f90 BIVARInterface.f90)
target_link_libraries( interp base mesh )

add_library( remesh SphereRemesh.f90 PlaneRemesh.f90 PlaneVorticity.f90 PlaneTracer.f90 BVEVorticity.f90 Tracers.f90 )
TARGET_LINK_LIBRARIES( remesh base mesh interp)

add_library( problemSet Tracers.f90 BVEVorticity.f90 SWEVorticityAndDivergence.f90 
						Advection2.f90 BVEDirectSum.f90 SWEDirectSum.f90 
						PlaneDirectSum.f90 PlaneVorticity.f90 PlaneTracer.f90  )
TARGET_LINK_LIBRARIES( problemSet base mesh remesh interp )

add_library( output VTKOutput.f90 PlaneOutput.f90 LatLonOutput.f90 ) 
TARGET_LINK_LIBRARIES( output base mesh interp)

ADD_EXECUTABLE( sweTestCase1 TestCase1.f90 )
TARGET_LINK_LIBRARIES(sweTestCase1 base mesh interp remesh problemSet output )

ADD_EXECUTABLE( advectGaussHills AdvectGaussHills.f90 )
TARGET_LINK_LIBRARIES(advectGaussHills base mesh interp remesh problemSet output )

ADD_EXECUTABLE( advectGaussHillsDirect AdvectGaussHillsDirect.f90 )
TARGET_LINK_LIBRARIES(advectGaussHillsDirect base mesh interp remesh problemSet output )

ADD_EXECUTABLE( advectMovingVortices AdvectMovingVortices.f90 )
TARGET_LINK_LIBRARIES(advectMovingVortices base mesh interp remesh problemSet output )

#ADD_EXECUTABLE( advectMovingVorticesWithVorticityRefinement AdvectMovingVortsRefineVorticity.f90 )
#TARGET_LINK_LIBRARIES(advectMovingVorticesWithVorticityRefinement base mesh interp remesh problemSet output )

ADD_EXECUTABLE( advectMovingVorticesDirect AdvectMovingVortsDirect.f90 )
TARGET_LINK_LIBRARIES(advectMovingVorticesDirect base mesh interp remesh problemSet output )

ADD_EXECUTABLE( advectSlottedCylinders AdvectSlottedCylinders.f90  )
TARGET_LINK_LIBRARIES(advectSlottedCylinders base mesh interp remesh problemSet output )

ADD_EXECUTABLE( advectSlottedCylindersDirect AdvectSlottedCylindersDirect.f90 ) 
TARGET_LINK_LIBRARIES(advectSlottedCylindersDirect base mesh interp remesh problemSet output )

ADD_EXECUTABLE( advectCorrelatedTracers AdvectCorrelated.f90 ) 
TARGET_LINK_LIBRARIES(advectCorrelatedTracers base mesh interp remesh problemSet output )

ADD_EXECUTABLE( advectCorrelatedTracersDirect AdvectCorrelatedDirect.f90 )
TARGET_LINK_LIBRARIES(advectCorrelatedTracersDirect base mesh interp remesh problemSet output )

#ADD_EXECUTABLE( solidBodyRotationBVE BVESolidBodyRotation.f90 ) 
#TARGET_LINK_LIBRARIES(solidBodyRotationBVE base mesh interp remesh problemSet output )

#ADD_EXECUTABLE( singleGaussianVortexBVE BVESingleGaussianVortex.f90 )
#TARGET_LINK_LIBRARIES(singleGaussianVortexBVE base mesh interp remesh problemSet output )

#ADD_EXECUTABLE( rossbyHaurwitz4Wave BVERH4.f90 ) 
#TARGET_LINK_LIBRARIES(rossbyHaurwitz4Wave base mesh interp remesh problemSet output )

ADD_EXECUTABLE(sweTestCase2 SWETestCase2.f90 ) 
TARGET_LINK_LIBRARIES(sweTestCase2 base mesh interp remesh problemSet output )

ADD_EXECUTABLE( lambDipole LambDipole.f90 )
TARGET_LINK_LIBRARIES(lambDipole base mesh interp remesh problemSet output )

ADD_EXECUTABLE( twoDipoles TwoDipoles.f90 )
TARGET_LINK_LIBRARIES(twoDipoles base mesh interp remesh problemSet output )

#ADD_EXECUTABLE( convertVTK2NCL ConvertVTKtoNativeNCL.f90 )
#TARGET_LINK_LIBRARIES(convertVTK2NCL base )

#ADD_EXECUTABLE( interpVTK2NCL InterpDataFromVTKtoNCL.f90 )
#TARGET_LINK_LIBRARIES(interpVTK2NCL base interp )

ADD_EXECUTABLE( sweDivergenceTerms SWEDivergenceEqnTerms.f90 )
TARGET_LINK_LIBRARIES(sweDivergenceTerms base mesh interp remesh problemSet output )

ADD_EXECUTABLE( planeTest planeTester.f90 )
TARGET_LINK_LIBRARIES(planeTest base mesh interp remesh problemSet output )

ADD_EXECUTABLE( icosTriTest IcosTriTest2.f90 )
TARGET_LINK_LIBRARIES(icosTriTest base mesh interp remesh problemSet output )

ADD_EXECUTABLE( cubedSphereTest CubedSphereTest2.f90 )
TARGET_LINK_LIBRARIES(cubedSphereTest base mesh interp remesh problemSet output )



