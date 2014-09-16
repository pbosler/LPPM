SHELL = /bin/bash

#MACHINE='FERRARI'
#MACHINE='TANK'
#MACHINE='VORTEX'
#MACHINE='LING'
MACHINE='WORK-LAPTOP'

## MAKEFILE FOR Lagrangian Particle/Panel Method on an Earth-Sized Sphere

#-----------------------------------------------------------------------------#
# 			MODIFY THIS SECTION FOR MACHINE-SPECIFIC PATHS					  #
#																			  #
ifeq ($(MACHINE),'FERRARI') 
# FERRARI LAPTOP #	
  FF = ifort
  FF_FLAGS = -O0 -g -check bounds -check pointers -check uninit -traceback -warn all -debug extended -openmp
  #FF_FLAGS= -O2 -openmp -warn all #-opt_report 1
  VTK_INCLUDE=/usr/local/include/vtk-5.8
  VTK_LIB_DIR=/usr/local/lib/vtk-5.8
  MKLROOT=/opt/intel/mkl
  MKL_LINK=-L$(MKLROOT)/lib $(MKLROOT)/lib/libmkl_lapack95_lp64.a -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lpthread -lm
  MKL_COMPILE=-openmp -I$(MKLROOT)/include/intel64/lp64 -I$(MKLROOT)/include
else ifeq ($(MACHINE),'VORTEX')
# vortex.math.lsa.umich.edu
  FF = ifort
  #FF_FLAGS = -O0 -g -check bounds -check pointers -check uninit -traceback -warn all -debug extended -openmp
  FF_FLAGS= -O2 -openmp -warn all #-opt_report 1
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
else ifeq ($(MACHINE),'LING')  
  FF=gfortran
  FF_FLAGS=-O2 -fopenmp -ffree-line-length-none
else ifeq ($(MACHINE),'WORK-LAPTOP')
  FF = ifort
  #FF_FLAGS = -O0 -g -check bounds -check pointers -check uninit -traceback -warn all -debug extended -openmp
  FF_FLAGS= -O2 -openmp -warn all #-opt_report 1
  MKL_ROOT=/opt/intel/mkl
  MKL_COMPILE=-openmp -I$(MKLROOT)/include/lp64 -mkl=parallel 
  MKL_LINK=$(MKLROOT)/lib/libmkl_blas95_lp64 $(MKLROOT)/lib/libmkl_lapack95_lp64 -lpthread -lm
  NETCDF_INCLUDE_DIR=-I/opt/local/include
  NETCDF_LIB_DIR=-L/opt/local/lib
  NETCDF_LIBS=-lnetcdff
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

MESH_OBJS = $(BASE_OBJS) Particles.o Edges.o Panels.o SphereMesh2.o 
mesh_objs: $(BASE_OBJS) Particles.o Edges.o Panels.o SphereMesh2.o 

INTERP_OBJS = $(BASE_OBJS) $(MESH_OBJS) ssrfpack.o stripack.o STRIPACKInterface2.o SSRFPACKInterface2.o
interp_objs: $(BASE_OBJS) $(MESH_OBJS) ssrfpack.o stripack.o STRIPACKInterface2.o SSRFPACKInterface2.o

OUTPUT_OBJS = VTKOutput.o LatLonOutput.o $(BASE_OBJS) $(MESH_OBJS) $(INTERP_OBJS)

TEST_CASE_OBJS = Tracers.o BVEVorticity.o SWEVorticityAndDivergence.o 

ADVECTION_OBJS = $(BASE_OBJS) $(MESH_OBJS) $(INTERP_OBJS) $(OUTPUT_OBJS) $(TEST_CASE_OBJS) Advection2.o SphereRemesh.o

BVE_OBJS = $(BASE_OBJS) $(MESH_OBJS) $(INTERP_OBJS) $(OUTPUT_OBJS) $(TEST_CASE_OBJS) BVEDirectSum.o SphereRemesh.o

PLANE_MESH = $(BASE_OBJS) Particles.o Edges.o Panels.o PlaneMesh.o
PLANE_REMESH = $(BASE_OBJS) $(PLANE_MESH) bivar.o BIVARInterface.o PlaneRemesh.o
PLANE_OUTPUT = $(BASE_OBJS) $(PLANE_MESH) PlaneOutput.o
PLANE_RUNS = $(BASE_OBJS) $(PLANE_MESH) $(PLANE_REMESH) $(PLANE_OUTPUT) PlaneDirectSum.o PlaneTracer.o PlaneVorticity.o

SWE_OBJS = $(BASE_OBJS) $(MESH_OBJS) $(INTERP_OBJS) $(OUTPUT_OBJS) SWEDirectSum.o SphereRemesh.o $(TEST_CASE_OBJS)


#############################################################
## LPPM MODEL RUNS
#############################################################

testCase1MPI.exe: TestCase1.o $(ADVECTION_OBJS)
	$(FF) $(FF_FLAGS) -o $@ $^ `mpif90 -showme:link`
advectGaussHillsMPI.exe: AdvectGaussHills.o $(ADVECTION_OBJS)
	$(FF) $(FF_FLAGS) -o $@ $^ `mpif90 -showme:link`
advectGaussHillsDirectMPI.exe: AdvectGaussHillsDirect.o $(ADVECTION_OBJS)
	$(FF) $(FF_FLAGS) -o $@ $^ `mpif90 -showme:link`	
advectMovingVorticesMPI.exe: AdvectMovingVortices.o $(ADVECTION_OBJS)
	$(FF) $(FF_FLAGS) -o $@ $^ `mpif90 -showme:link`
advectMovingVorticesWithVorticityRefinementMPI.exe: AdvectMovingVortsRefineVorticity.o $(ADVECTION_OBJS)
	$(FF) $(FF_FLAGS) -o $@ $^ `mpif90 -showme:link`
advectSlottedCylindersMPI.exe: AdvectSlottedCylinders.o $(ADVECTION_OBJS)
	$(FF) $(FF_FLAGS) -o $@ $^ `mpif90 -showme:link`
advectSlottedCylindersDirectMPI.exe: AdvectSlottedCylindersDirect.o $(ADVECTION_OBJS)
	$(FF) $(FF_FLAGS) -o $@ $^ `mpif90 -showme:link`		
advectCorrelatedMPI.exe: AdvectCorrelated.o $(ADVECTION_OBJS)
	$(FF) $(FF_FLAGS) -o $@ $^ `mpif90 -showme:link`
advectCorrelatedDirectMPI.exe: AdvectCorrelatedDirect.o $(ADVECTION_OBJS)
	$(FF) $(FF_FLAGS) -o $@ $^ `mpif90 -showme:link`		
advectCascadeMPI.exe: AdvectionCascade.o $(ADVECTION_OBJS)
	$(FF) $(FF_FLAGS) -o $@ $^ `mpif90 -showme:link`	
solidBodyRotationMPI.exe: BVESolidBodyRotation.o $(BVE_OBJS)	
	$(FF) $(FF_FLAGS) -o $@ $^ `mpif90 -showme:link` 
singleGaussianVortexMPI.exe: BVESingleGaussianVortex.o $(BVE_OBJS)
	$(FF) $(FF_FLAGS) -o $@ $^ `mpif90 -showme:link` 
rossbyHaurwitz4waveMPI.exe: BVERH4.o $(BVE_OBJS) 
	$(FF) $(FF_FLAGS) -o $@ $^ `mpif90 -showme:link` $(MKL_LINK)
sweTestCase2MPI.exe: SWETestCase2.o $(SWE_OBJS)
	$(FF) $(FF_FLAGS) -o $@ $^ `mpif90 -showme:link` $(MKL_COMPILE) $(MKL_LINK)	
lambDipoleMPI.exe: LambDipole.o $(PLANE_RUNS)
	$(FF) $(FF_FLAGS) -o $@ $^ `mpif90 -showme:link` 
twoDipolesMPI.exe: TwoDipoles.o $(PLANE_RUNS)
	$(FF) $(FF_FLAGS) -o $@ $^ `mpif90 -showme:link` 
reversibleDipolesMPI.exe: ReversibleDipoles.o $(PLANE_RUNS)
	$(FF) $(FF_FLAGS) -o $@ $^ `mpif90 -showme:link`	
planeAdvectRotationMPI.exe: PlaneRotationalAdvection.o $(PLANE_RUNS)
	$(FF) $(FF_FLAGS) -o $@ $^ `mpif90 -showme:link`			

#############################################################
## Utilities
#############################################################
convertVTKtoNCL.exe: ConvertVTKtoNativeNCL.o NumberKinds3.o SphereGeom3.o
	$(FF) -O3 -o $@ $^
ConvertVTKtoNativeNCL.o: ConvertVTKtoNativeNCL.f90 NumberKinds3.o SphereGeom3.o
	$(FF) -c -O3 $<
interpVTKtoNCL.exe: InterpDataFromVTKtoNCL.o NumberKinds3.o SphereGeom3.o $(INTERP_OBJS)
	$(FF) -O3 -openmp -o $@ $^ 
InterpDataFromVTKToNCL.o: InterpDataFromVTKtoNCL.f90 NumberKinds3.o SphereGeom3.o $(INTERP_OBJS)
	$(FF) -O3 -c InterpDataFromVTKtoNCL.f90

#############################################################
## LPPM MODEL OBJECT FILES
#############################################################

AdvectGaussHills.o: AdvectGaussHills.f90 $(ADVECTION_OBJS)
AdvectGaussHillsDirect.o: AdvectGaussHillsDirect.f90 $(ADVECTION_OBJS)
AdvectMovingVortices.o: AdvectMovingVortices.f90 $(ADVECTION_OBJS)
AdvectMovingVortsRefineVorticity.o: AdvectMovingVortsRefineVorticity.f90 $(ADVECTION_OBJS)
AdvectSlottedCylinders.o: AdvectSlottedCylinders.f90 $(ADVECTION_OBJS)
AdvectSlottedCylindersDirect.o: AdvectSlottedCylindersDirect.f90 $(ADVECTION_OBJS)
AdvectCorrelated.o: AdvectCorrelated.f90 $(ADVECTION_OBJS)
AdvectCorrelatedDirect.o: AdvectCorrelated.f90 $(ADVECTION_OBJS)
AdvectionCascade.o: AdvectionCascade.f90 $(ADVECTION_OBJS)
TestCase1.o: TestCase1.f90 $(ADVECTION_OBJS)
BVESolidBodyRotation.o: BVESolidBodyRotation.f90 $(BVE_OBJS)
BVESingleGaussianVortex.o: BVESingleGaussianVortex.f90 $(BVE_OBJS) 
BVERH4.o: BVERH4.f90 $(BVE_OBJS) 
SWETestCase2.o: SWETestCase2.f90 $(SWE_OBJS)
LambDipole.o: LambDipole.f90 $(PLANE_RUNS)
TwoDipoles.o: TwoDipoles.f90 $(PLANE_RUNS)
PlaneRotationalAdvection.o: PlaneRotationalAdvection.f90 $(PLANE_RUNS)
ReversibleDipoles.o: ReversibleDipoles.f90 $(PLANE_RUNS)

#############################################################
## UNIT TEST EXECUTABLES
#############################################################

cubedSphereTestSerial.exe: CubedSphereTest2.o $(BASE_OBJS) $(MESH_OBJS) $(INTERP_OBJS) $(OUTPUT_OBJS) Advection2.o
	$(FF) $(FF_FLAGS) -o $@  $^ `mpif90 -showme:link` $(MKL_COMPILE)
icosTriTestMPI.exe: IcosTriTest2.o	$(BASE_OBJS) $(MESH_OBJS) $(INTERP_OBJS) $(OUTPUT_OBJS)
	$(FF) $(FF_FLAGS) -o $@  $^ `mpif90 -showme:link` $(MKL_COMPILE)	
planeTestMPI.exe: planeTester.o $(BASE_OBJS) $(MESH_OBJS) $(OUTPUT_OBJS) $(TEST_CASE_OBJS) PlaneDirectSum.o PlaneRemesh.o bivar.o BIVARInterface.o
	$(FF) $(FF_FLAGS) -o $@ $^ `mpif90 -showme:link`
sweDivergenceTerms.exe: SWEDivergenceEqnTerms.o $(BASE_OBJS) $(MESH_OBJS) $(OUTPUT_OBJS) 
	$(FF) $(FF_FLAGS) -o $@ $^ $(MKL_LINK) `mpif90 -showme:link`

#############################################################
## UNIT TEST OBJECT FILES
#############################################################

CubedSphereTest2.o: CubedSphereTest2.f90 $(BASE_OBJS) $(MESH_OBJS) $(INTERP_OBJS) $(OUTPUT_OBJS)
IcosTriTest2.o: IcosTriTest2.f90 $(BASE_OBJS) $(MESH_OBJS) $(INTERP_OBJS) $(OUTPUT_OBJS)
planeTester.o: planeTester.f90 $(PLANE_RUNS)
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
bivar.o: bivar.f90 NumberKinds3.o
BIVARInterface.o: BIVARInterface.f90 $(BASE_OBJS) $(PLANE_MESH)
STRIPACKInterface2.o: STRIPACKInterface2.f90 $(BASE_OBJS) $(MESH_OBJS) stripack.o 
SSRFPACKInterface2.o: SSRFPACKInterface2.f90 $(BASE_OBJS) $(MESH_OBJS) STRIPACKInterface2.o ssrfpack.o
VTKOutput.o: VTKOutput.f90 $(BASE_OBJS) $(MESH_OBJS)
#RefineRemesh2.o: RefineRemesh2.f90 $(BASE_OBJS) $(MESH_OBJS) $(INTERP_OBJS) Tracers.o BVEVorticity.o
Tracers.o: Tracers.f90 $(BASE_OBJS) $(MESH_OBJS)
BVEVorticity.o: BVEVorticity.f90 $(BASE_OBJS) $(MESH_OBJS)
PlaneVorticity.o: PlaneVorticity.f90 $(BASE_OBJS) $(PLANE_MESH)
PlaneTracer.o: PlaneTracer.f90 $(BASE_OBJS) $(PLANE_MESH)
Advection2.o: Advection2.f90 $(BASE_OBJS) $(MESH_OBJS)
BVEDirectSum.o: BVEDirectSum.f90 $(BASE_OBJS) $(MESH_OBJS)
PlaneDirectSum.o: PlaneDirectSum.f90 $(BASE_OBJS) $(PLANE_MESH)
PlaneRemesh.o: PlaneRemesh.f90 $(BASE_OBJS) $(PLANE_MESH) PlaneTracer.o PlaneVorticity.o bivar.o BIVARInterface.o
#ReferenceSphere.o: ReferenceSphere.f90 $(BASE_OBJS) $(MESH_OBJS) $(INTERP_OBJS) RefineRemesh2.o
LatLonOutput.o: LatLonOutput.f90 $(BASE_OBJS) $(MESH_OBJS) $(INTERP_OBJS)
SWEDirectSum.o: SWEDirectSum.f90 $(BASE_OBJS) $(MESH_OBJS) $(INTERP_OBJS)
SWEVorticityAndDivergence.o: SWEVorticityAndDivergence.f90 $(BASE_OBJS) $(MESH_OBJS)
SphereRemesh.o: SphereRemesh.f90 $(BASE_OBJS) $(MESH_OBJS) Tracers.o BVEVorticity.o $(INTERP_OBJS)

#############################################################
## VTK EXECUTABLES
#############################################################

plotNewRemeshingFrames.exe: vtkPlotNewRemeshingFrames.o
	g++ -O2 -Wno-deprecated -o $@ $< -I$(VTK_INCLUDE) -L$(VTK_LIB_DIR) $(VTK_LIBS)
plotRH4WaveFrames.exe: vtkPlotRHWaveFrames.o ModelLookupTables.o
	g++ -O2 -Wno-deprecated -o $@ $^ -I$(VTK_INCLUDE) -L$(VTK_LIB_DIR) $(VTK_LIBS)
plotRH4WaveWithScalar.exe: vtkPlotRHWaveWithLatScalar.o ModelLookupTables.o	
	g++ -O2 -Wno-deprecated -o $@ $^ -I$(VTK_INCLUDE) -L$(VTK_LIB_DIR) $(VTK_LIBS)
plotPolarVortex.exe: vtkPlotPolarVortexFrames.o ModelLookupTables.o
	g++ -O2 -Wno-deprecated -o $@ $^ -I$(VTK_INCLUDE) -L$(VTK_LIB_DIR) $(VTK_LIBS) -fopenmp
#############################################################
## VTK OBJECTS
#############################################################

ModelLookupTables.o: ModelLookupTables.cpp ModelLookupTables.h 
	g++ -c -Wno-deprecated -O2 $< -I$(VTK_INCLUDE)
vtkPlotNewRemeshingFrames.o: vtkPlotNewRemeshingFrames.cpp
	g++ -c -Wno-deprecated -O2 $< -I$(VTK_INCLUDE)
vtkPlotBigRHWave.o: vtkPlotBigRHWave.cpp 
	g++ -c -Wno-deprecated -O2 $^ -I$(VTK_INCLUDE)	
vtkPlotRHWaveFrames.o: vtkPlotRHWaveFrames.cpp
	g++ -c -Wno-deprecated -O2 $< -I$(VTK_INCLUDE)
vtkPlotRHWaveWithLatScalar.o: vtkPlotRHWaveWithLatScalar.cpp
	g++ -c -Wno-deprecated -O2 $< -I$(VTK_INCLUDE)
vtkPlotPolarVortexFrames.o: vtkPlotPolarVortexFrames.cpp
	g++ -c -Wno-deprecated -O2 $< -I$(VTK_INCLUDE) -fopenmp