SHELL = /bin/bash

#MACHINE='FERRARI'
MACHINE='TANK'
#MACHINE='VORTEX'

## MAKEFILE FOR Lagrangian Particle/Panel Method on an Earth-Sized Sphere

#-----------------------------------------------------------------------------#
# 			MODIFY THIS SECTION FOR MACHINE-SPECIFIC PATHS					  #
#																			  #
ifeq ($(MACHINE),'FERRARI') 
# FERRARI LAPTOP #	
  FF = ifort
  #FF_FLAGS = -O0 -g -check bounds -check pointers -check uninit -traceback -warn all -debug extended -openmp
  FF_FLAGS= -O2 -openmp -warn all -opt_report 1
  VTK_INCLUDE=/usr/local/include/vtk-5.8
  VTK_LIB_DIR=/usr/local/lib/vtk-5.8
  MKLROOT=/opt/intel/mkl
  MKL_LINK=-L$(MKLROOT)/lib $(MKLROOT)/lib/libmkl_lapack95_lp64.a -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lpthread -lm
  MKL_COMPILE=-openmp -I$(MKLROOT)/include/intel64/lp64 -I$(MKLROOT)/include
else ifeq ($(MACHINE),'VORTEX')
# vortex.math.lsa.umich.edu
  FF = ifort
  #FF_FLAGS = -O0 -g -check bounds -check pointers -check uninit -traceback -warn all -debug extended -openmp
  FF_FLAGS= -O2 -openmp -warn all -opt_report 1
  MKL_ROOT=/usr/local/intel/Compiler/11.1/056/mkl
  MKL_LINK=-L$(MKL_ROOT)/lib $(MKL_ROOT)/lib/libmkl_lapack95_lp64.a -lmkl_intel_lp64 -lmk_intel_thread -lmkl_core -lpthread -lm
  MKL_COMPILE=-openmp -I$(MKL_ROOT)/include/intel64/lp64 -I$(MKL_ROOT)/include
else ifeq ($(MACHINE),'TANK')
# TANK DESKTOP #
  FF = ifort
  #FF_FLAGS = -O0 -g -check bounds -check pointers -check uninit -traceback -warn all -debug extended -openmp
  FF_FLAGS= -O2 -openmp -warn all #-opt_report 1
  VTK_INCLUDE=/usr/local/include/vtk-5.10
  VTK_LIB_DIR=/usr/local/lib/vtk-5.10
  VTK_LIBS=-lvtkCommon -lvtkGraphics -lvtkRendering -lvtkViews -lvtkWidgets -lvtkImaging -lvtkHybrid -lvtkIO -lvtkFiltering
  MKLROOT=/opt/intel/mkl
  MKL_THREADING_LAYER=intel
  MKL_INTERFACE_LAYER=lp64
  MKL_LINK=-L$(MKLROOT)/lib $(MKLROOT)/lib/libmkl_lapack95_lp64.a \
  -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -openmp -lpthread -lm
  MKL_COMPILE=-openmp -I/opt/intel/mkl/include/intel64/lp64 -I/opt/intel/mkl/include 
endif
#-----------------------------------------------------------------------------#

#############################################################
## MAKE RULES
#############################################################
%.o : %.f90
	$(FF) -c $(FF_FLAGS) $< `mpif90 -showme:compile`
%.o : %.f
	$(FF) -c $(FF_FLAGS) $<
clean:
	rm *.o *.mod *genmod.f90
cleanx:
	rm *.exe		

#############################################################
## LOCAL MAKE VARIABLES
#############################################################

BASE_OBJS = NumberKinds3.o OutputWriter2.o Logger2.o SphereGeom3.o IntegerList2.o PlaneGeometry.o
base_objs: NumberKinds3.o OutputWriter2.o Logger2.o SphereGeom3.o IntegerList2.o PlaneGeometry.o

MESH_OBJS = $(BASE_OBJS) Particles.o Edges.o Panels.o SphereMesh2.o PlaneMesh.o
mesh_objs: $(BASE_OBJS) Particles.o Edges.o Panels.o SphereMesh2.o PlaneMesh.o

INTERP_OBJS = $(BASE_OBJS) $(MESH_OBJS) ssrfpack.o stripack.o STRIPACKInterface2.o SSRFPACKInterface2.o
interp_objs: $(BASE_OBJS) $(MESH_OBJS) ssrfpack.o stripack.o STRIPACKInterface2.o SSRFPACKInterface2.o

OUTPUT_OBJS = VTKOutput.o LatLonOutput.o PlaneOutput.o $(BASE_OBJS) $(MESH_OBJS) $(INTERP_OBJS)

TEST_CASE_OBJS = Tracers.o BVEVorticity.o SWEVorticityAndDivergence.o PlaneVorticity.o

ADVECTION_OBJS = $(BASE_OBJS) $(MESH_OBJS) $(INTERP_OBJS) $(OUTPUT_OBJS) $(TEST_CASE_OBJS) Advection2.o RefineRemesh2.o

BVE_OBJS = $(BASE_OBJS) $(MESH_OBJS) $(INTERP_OBJS) $(OUTPUT_OBJS) $(TEST_CASE_OBJS) BVEDirectSum.o RefineRemesh2.o

PLOTTING = ModelLookupTables.o ModelLookupTables.h

SWE_OBJS = $(BASE_OBJS) $(MESH_OBJS) $(INTERP_OBJS) $(OUTPUT_OBJS) SWEDirectSum.o RefineRemesh2.o $(TEST_CASE_OBJS)


#############################################################
## LPPM MODEL RUNS
#############################################################

testCase1MPI.exe: TestCase1.o $(ADVECTION_OBJS)
	$(FF) $(FF_FLAGS) -o $@ $^ `mpif90 -showme:link` $(MKL_COMPILE)
solidBodyRotationMPI.exe: BVESolidBodyRotation.o $(BVE_OBJS)	
	$(FF) $(FF_FLAGS) -o $@ $^ `mpif90 -showme:link` $(MKL_COMPILE)
singleGaussianVortexMPI.exe: BVESingleGaussianVortex.o $(BVE_OBJS) ReferenceSphere.o
	$(FF) $(FF_FLAGS) -o $@ $^ `mpif90 -showme:link` $(MKL_COMPILE)
rossbyHaurwitz4waveMPI.exe: BVERH4.o $(BVE_OBJS) ReferenceSphere.o	 
	$(FF) $(FF_FLAGS) -o $@ $^ `mpif90 -showme:link` $(MKL_COMPILE)
sweTestCase2MPI.exe: SWETestCase2.o $(SWE_OBJS)
	$(FF) $(FF_FLAGS) -o $@ $^ `mpif90 -showme:link` $(MKL_COMPILE)	
	
#############################################################
## LPPM MODEL OBJECT FILES
#############################################################

TestCase1.o: TestCase1.f90 $(ADVECTION_OBJS)
BVESolidBodyRotation.o: BVESolidBodyRotation.f90 $(BVE_OBJS)
BVESingleGaussianVortex.o: BVESingleGaussianVortex.f90 $(BVE_OBJS) ReferenceSphere.o
BVERH4.o: BVERH4.f90 $(BVE_OBJS) ReferenceSphere.o
SWETestCase2.o: SWETestCase2.f90 $(SWE_OBJS)

#############################################################
## UNIT TEST EXECUTABLES
#############################################################

cubedSphereTestSerial.exe: CubedSphereTest2.o $(BASE_OBJS) $(MESH_OBJS) $(INTERP_OBJS) $(OUTPUT_OBJS) Advection2.o
	$(FF) $(FF_FLAGS) -o $@  $^ `mpif90 -showme:link` $(MKL_COMPILE)	
planeTestSerial.exe: planeTester.o $(BASE_OBJS) $(MESH_OBJS) $(OUTPUT_OBJS) $(TEST_CASE_OBJS)
	$(FF) $(FF_FLAGS) -o $@ $^ 
sweDivergenceTerms.exe: SWEDivergenceEqnTerms.o $(BASE_OBJS) $(MESH_OBJS) $(OUTPUT_OBJS) 
	$(FF) $(FF_FLAGS) -o $@ $^ $(MKL_LINK) `mpif90 -showme:link`

#############################################################
## UNIT TEST OBJECT FILES
#############################################################

CubedSphereTest2.o: CubedSphereTest2.f90 $(BASE_OBJS) $(MESH_OBJS) $(INTERP_OBJS) $(OUTPUT_OBJS)
planeTester.o: planeTester.f90 $(BASE_OBJS) $(MESH_OBJS) $(OUTPUT_OBJS) $(TEST_CASE_OBJS)
SWEDivergenceEqnTerms.o: SWEDivergenceEqnTerms.f90 $(BASE_OBJS) $(MESH_OBJS) $(OUTPUT_OBJS)
	$(FF) $(FF_FLAGS) -c $< $(MKL_COMPILE) `mpif90 -showme:compile`

#############################################################
## MODULES
#############################################################

NumberKinds3.o: NumberKinds3.f90
IntegerList2.o: IntegerList2.f90 NumberKinds3.o 
OutputWriter2.o: OutputWriter2.f90 NumberKinds3.o 
SphereGeom3.o: SphereGeom3.f90 NumberKinds3.o 
PlaneGeometry.o: PlaneGeometry.f90 NumberKinds3.o
Logger2.o: Logger2.f90 NumberKinds3.o OutputWriter2.o 
Particles.o: Particles.f90 $(BASE_OBJS)
Edges.o: Edges.f90 $(BASE_OBJS)
Panels.o: Panels.f90 $(BASE_OBJS)
SphereMesh2.o: SphereMesh2.f90 Particles.o Edges.o Panels.o $(BASE_OBJS)
PlaneMesh.o: PlaneMesh.f90 Particles.o Edges.o Panels.o $(BASE_OBJS)
PlaneOutput.o: PlaneOutput.f90 PlaneMesh.o Particles.o Edges.o Panels.o $(BASE_OBJS)
ssrfpack.o: ssrfpack.f
stripack.o: stripack.f
STRIPACKInterface2.o: STRIPACKInterface2.f90 $(BASE_OBJS) $(MESH_OBJS) stripack.o 
SSRFPACKInterface2.o: SSRFPACKInterface2.f90 $(BASE_OBJS) $(MESH_OBJS) STRIPACKInterface2.o ssrfpack.o
VTKOutput.o: VTKOutput.f90 $(BASE_OBJS) $(MESH_OBJS)
RefineRemesh2.o: RefineRemesh2.f90 $(BASE_OBJS) $(MESH_OBJS) $(INTERP_OBJS) Tracers.o BVEVorticity.o
Tracers.o: Tracers.f90 $(BASE_OBJS) $(MESH_OBJS)
BVEVorticity.o: BVEVorticity.f90 $(BASE_OBJS) $(MESH_OBJS)
PlaneVorticity.o: PlaneVorticity.f90 $(BASE_OBJS) $(MESH_OBJS)
Advection2.o: Advection2.f90 $(BASE_OBJS) $(MESH_OBJS)
BVEDirectSum.o: BVEDirectSum.f90 $(BASE_OBJS) $(MESH_OBJS)
ReferenceSphere.o: ReferenceSphere.f90 $(BASE_OBJS) $(MESH_OBJS) $(INTERP_OBJS) RefineRemesh2.o
LatLonOutput.o: LatLonOutput.f90 $(BASE_OBJS) $(MESH_OBJS) $(INTERP_OBJS)
SWEDirectSum.o: SWEDirectSum.f90 $(BASE_OBJS) $(MESH_OBJS) $(INTERP_OBJS)
SWEVorticityAndDivergence.o: SWEVorticityAndDivergence.f90 $(BASE_OBJS) $(MESH_OBJS)

#############################################################
## VTK EXECUTABLES
#############################################################

plotNewRemeshingFrames.exe: vtkPlotNewRemeshingFrames.o
	g++ -O2 -Wno-deprecated -o $@ $< -I$(VTK_INCLUDE) -L$(VTK_LIB_DIR) $(VTK_LIBS)
plotRH4WaveFrames.exe: vtkPlotBigRHWave.o $(PLOTTING)
	g++ -O2 -Wno-deprecated -o $@ $^ -I$(VTK_INCLUDE) -L$(VTK_LIB_DIR) $(VTK_LIBS)

#############################################################
## VTK OBJECTS
#############################################################

ModelLookupTables.o: ModelLookupTables.cpp ModelLookupTables.h 
	g++ -c -Wno-deprecated -O2 $< -I$(VTK_INCLUDE)
vtkPlotNewRemeshingFrames.o: vtkPlotNewRemeshingFrames.cpp
	g++ -c -Wno-deprecated -O2 $< -I$(VTK_INCLUDE)
vtkPlotBigRHWave.o: vtkPlotBigRHWave.cpp $(PLOTTING)
	g++ -c -Wno-deprecated -O2 $^ -I$(VTK_INCLUDE)	