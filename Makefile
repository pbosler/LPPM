SHELL = /bin/bash

## MAKEFILE FOR Lagrangian Particle/Panel Method on an Earth-Sized Sphere

# MODIFY THIS SECTION FOR EACH PLATFORM #

#----------------#
# FERRARI LAPTOP #
#FF=ifort
#FF_FLAGS=-O0 -g -check bounds -check pointers -check uninit -traceback -warn all -debug extended -openmp
#FF_FLAGS= -O2 -openmp -warn all
#VTK_INCLUDE=/usr/local/include/vtk-5.8
#VTK_LIB_DIR=/usr/local/lib/vtk-5.8
#MKLROOT=/opt/intel/mkl
#MKL_LINK=-L$(MKLROOT)/lib $(MKLROOT)/lib/libmkl_lapack95_lp64.a -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lpthread -lm
#MKL_COMPILE=-openmp -I$(MKLROOT)/include/intel64/lp64 -I$(MKLROOT)/include
#----------------#

#----------------#
# vortex.math.lsa.umich.edu
FF = ifort
#FF_FLAGS = -g -O0 -check bounds -check pointer -check uninit -traceback -warn all -debug extended -openmp
FF_FLAGS=-O2 -warn all -opt-report 1
MKL_ROOT=/usr/local/intel/Compiler/11.1/056/mkl
MKL_LINK=-L$(MKL_ROOT)/lib $(MKL_ROOT)/lib/libmkl_lapack95_lp64.a -lmkl_intel_lp64 -lmk_intel_thread -lmkl_core -lpthread -lm
MKL_COMPILE=-openmp -I$(MKL_ROOT)/include/intel64/lp64 -I$(MKL_ROOT)/include


#--------------#
# TANK DESKTOP #

#FF=ifort
#FF_FLAGS=-g -traceback -warn all -debug extended
#FF_FLAGS=-O2 -warn all -opt-report 1
#VTK_INCLUDE=/usr/local/include/vtk-5.10
#VTK_LIB_DIR=/usr/local/lib/vtk-5.10
#VTK_LIBS=-lvtkCommon -lvtkGraphics -lvtkRendering -lvtkViews -lvtkWidgets -lvtkImaging -lvtkHybrid -lvtkIO -lvtkFiltering
#MKLROOT=/opt/intel/mkl
#MKL_THREADING_LAYER=intel
#MKL_INTERFACE_LAYER=lp64
#MKL_LINK=-L$(MKLROOT)/lib $(MKLROOT)/lib/libmkl_lapack95_lp64.a \
#-lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -openmp -lpthread -lm
#MKL_COMPILE=-openmp -I/opt/intel/mkl/include/intel64/lp64 -I/opt/intel/mkl/include 

#--------------#

#############################################################
## MAKE RULES
#############################################################
%.o : %.f90
	$(FF) -c $(FF_FLAGS) $< `mpif90 -showme:compile` $(MKL_COMPILE)
%.o : %.f
	$(FF) -c $(FF_FLAGS) $<
clean:
	rm *.o *.mod *genmod.f90
cleanx:
	rm *.exe		

#############################################################
## LOCAL MAKE VARIABLES
#############################################################

BASE_OBJS = NumberKinds3.o OutputWriter2.o Logger2.o SphereGeom3.o IntegerList2.o 
base_objs: NumberKinds3.o OutputWriter2.o Logger2.o SphereGeom3.o IntegerList2.o
MESH_OBJS = $(BASE_OBJS) Particles.o Edges.o Panels.o SphereMesh2.o
mesh_objs: $(BASE_OBJS) Particles.o Edges.o Panels.o SphereMesh2.o 
INTERP_OBJS = $(BASE_OBJS) $(MESH_OBJS) ssrfpack.o stripack.o STRIPACKInterface2.o SSRFPACKInterface2.o
interp_objs: $(BASE_OBJS) $(MESH_OBJS) ssrfpack.o stripack.o STRIPACKInterface2.o SSRFPACKInterface2.o
OUTPUT_OBJS = VTKOutput.o $(BASE_OBJS) $(MESH_OBJS)
TEST_CASE_OBJS = Tracers.o BVEVorticity.o
ADVECTION_OBJS = $(BASE_OBJS) $(MESH_OBJS) $(INTERP_OBJS) $(OUTPUT_OBJS) $(TEST_CASE_OBJS) Advection2.o RefineRemesh2.o
BVE_OBJS = $(BASE_OBJS) $(MESH_OBJS) $(INTERP_OBJS) $(OUTPUT_OBJS) $(TEST_CASE_OBJS) BVEDirectSum.o RefineRemesh2.o


#############################################################
## LPPM MODEL RUNS
#############################################################
testCase1MPI.exe: TestCase1.o $(ADVECTION_OBJS)
	$(FF) $(FF_FLAGS) -o $@ $^ `mpif90 -showme:link` $(MKL_COMPILE)
solidBodyRotationMPI.exe: BVESolidBodyRotation.o $(BVE_OBJS)	
	$(FF) $(FF_FLAGS) -o $@ $^ `mpif90 -showme:link` $(MKL_COMPILE)
singleGaussianVortexMPI.exe: BVESingleGaussianVortex.o $(BVE_OBJS)	
	$(FF) $(FF_FLAGS) -o $@ $^ `mpif90 -showme:link` $(MKL_COMPILE)
	
#############################################################
## LPPM MODEL OBJECT FILES
#############################################################
TestCase1.o: TestCase1.f90 $(ADVECTION_OBJS)
BVESolidBodyRotation.o: BVESolidBodyRotation.f90 $(BVE_OBJS)
BVESingleGaussianVortex.o: BVESingleGaussianVortex.f90 $(BVE_OBJS)

#############################################################
## UNIT TEST EXECUTABLES
#############################################################
cubedSphereTestSerial.exe: CubedSphereTest2.o $(BASE_OBJS) $(MESH_OBJS) $(INTERP_OBJS) $(OUTPUT_OBJS) Advection2.o
	$(FF) $(FF_FLAGS) -o $@  $^ `mpif90 -showme:link` $(MKL_COMPILE)	

#############################################################
## UNIT TEST OBJECT FILES
#############################################################
CubedSphereTest2.o: CubedSphereTest2.f90 $(BASE_OBJS) $(MESH_OBJS) $(INTERP_OBJS) $(OUTPUT_OBJS)

#############################################################
## MODULES
#############################################################
NumberKinds3.o: NumberKinds3.f90
IntegerList2.o: NumberKinds3.o IntegerList2.f90
OutputWriter2.o: NumberKinds3.o OutputWriter2.f90
SphereGeom3.o: NumberKinds3.o SphereGeom3.f90
Logger2.o: NumberKinds3.o OutputWriter2.o Logger2.f90
Particles.o: Particles.f90 $(BASE_OBJS)
Edges.o: Edges.f90 $(BASE_OBJS)
Panels.o: Panels.f90 $(BASE_OBJS)
SphereMesh2.o: SphereMesh2.f90 Particles.o Edges.o Panels.o $(BASE_OBJS)
ssrfpack.o: ssrfpack.f
stripack.o: stripack.f
STRIPACKInterface2.o: STRIPACKInterface2.f90 $(BASE_OBJS) $(MESH_OBJS) stripack.o 
SSRFPACKInterface2.o: SSRFPACKInterface2.f90 $(BASE_OBJS) $(MESH_OBJS) STRIPACKInterface2.o ssrfpack.o
VTKOutput.o: VTKOutput.f90 $(BASE_OBJS) $(MESH_OBJS)
RefineRemesh2.o: RefineRemesh2.f90 $(BASE_OBJS) $(MESH_OBJS) $(INTERP_OBJS) Tracers.o BVEVorticity.o
Tracers.o: Tracers.f90 $(BASE_OBJS) $(MESH_OBJS)
BVEVorticity.o: BVEVorticity.f90 $(BASE_OBJS) $(MESH_OBJS)
Advection2.o: Advection2.f90 $(BASE_OBJS) $(MESH_OBJS)
BVEDirectSum.o: BVEDirectSum.f90 $(BASE_OBJS) $(MESH_OBJS)

#############################################################
## VTK EXECUTABLES
#############################################################

#############################################################
## VTK OBJECTS
#############################################################
VTKLookupTables.o: VTKLookupTables.cpp 
	g++ -c -O2 $< -I$(VTK_INCLUDE)
