LPPM : Lagrangian Particle / Panel Method
=========

Code associated with 

P. Bosler,  L. Wang,  R. Krasny, and C. Jablonowski,  2014.  
	"A Lagrangian particle/panel method for the barotropic vorticity equation on a rotating sphere," Fluid Dynamics Research,  46 : 031406.

Bosler, P. A., "Particle Methods for Geophysical Flow on the Sphere," PhD Thesis, University of Michigan, 2013.

This software provides a Lagrangian Particle-Panel Method for solving the barotropic vorticity equation
and the advection equation on the sphere and in the plane.

A shallow water equation solver is in development.

Design / Use
=========

Although it is written in modern Fortran rather than C++, the code is organized in an object-oriented style.
This means that for each task the code performs, there is a derived data type that corresponds to that task (like a C++ class).  
Each module file contains one derived data type and its associated methods.  
For example, PlaneMesh.f90 contains the data type for moving meshes in the plane and interfaces to all the functions that are 
implemented on PlaneMesh objects.  
The PlaneDirectSum.f90 module defines the Runge-Kutta time stepping data type for plane meshes, and contains all of the timestepping functions.

Each object must be created in memory using its New method, many objects must then be initialized using a call to an Init or Initialize subroutine.  
For clean code, every time a New subroutine is called on an object, a corresponding call to that object's Delete subroutine should be used when the object is no longer needed or goes out of scope.

Example: 
---------
To define a Rossby-Haurwitz wave vorticity distribution for a sphere mesh, a user would first allocate the VorticitySetup object by calling its New subroutine.
Then the wave would be initialized with user-defined parameters (most commonly read from a namelist file) by calling the VorticitySetup module's InitRHWave routine.
Finally, the wave is defined and users can define the vorticity on a sphere mesh by calling the SetRHWaveOnMesh subroutine.
  
	call New(rhWave, nIntegers, nReals)
	call InitRH4Wave(rhWave, backgoundWindSpeed, amplitude)
	call SetRH4WaveOnMesh( sphere, rhWave)
	
Code runs...

	call Delete(rhWave)	

To create a new test case:
---------------
The first step is to create the object that will control your test case, either a tracer or a vorticity distribution, in the Tracers.f90 or BVEVorticity.f90 files.
Test cases are defined by BVESetup or TracerSetup objects.  
These objects hold the integers and real numbers that define your test case; numbers that, for example, define the amplitude of a Rossby-Haurwitz wave or the maximum value of an advected tracer.
To create the object the first step is to allocate the memory to hold these parameters using the "New" subroutine.  

	call New(rhWave, RHWAVE_N_INT, RHWAVE_N_REAL)

RHWAVE_N_INIT and RHWAVE_N_REAL are integer constants defined in the BVESetup.f90 module; their values are 0 and 2 for the RH4 test case,
which indicates that this test case is defined by 2 real numbers (the ampilitude and phase speeed of the wave).  
Now that memory is allocated to hold these parameters, the test case must be initialized from user input (typically read from a namelist file).
If the variables backgroundWindSpeed and amplitude are set by the user, the test case is defind by calling the "Init" subroutine on the test case object.

	call InitRH4Wave(rhWave, backgroundWindSpeed, amplitude)

Finally and most importantly is the subroutine which translates a test case into a tracer or vorticity distribution over a mesh.
These subroutine must all have the same interface -- only 2 arguments, the first being a mesh object and second being the test case definition.
This interface is required in order for the LagrangianRemeshing subroutines to work properly.

	call SetRH4WaveOnMesh( aMesh, rhWave)

SetRH4WaveOnMesh calls the mathematical functions that define the vorticity distribution of a Rossby-Haurwitz wave, given the parameters defined by the rhWave test case object.

To create a new test case then, users must first decide how many integers and reals are required to define their test case, and allocate the appropriate amount of memory with the "New" subroutine.
Then the test case object must be initialized with an "Init" subroutine to define the parameters in the test case object's memory.
Finally, users must write a SetTestCaseOnMesh subroutine and the associated mathematical functions that define their test case.

	
Requirements
=========
The code requires a modern Fortran compiler with an OpenMPI (<http://www.open-mpi.org>) distribution for parallel computing.  
Without OpenMPI serial computations are still possible, but the code is designed to use distributed memory parallelism.

Currently the software is written for Intel compilers and makes use of the Intel MKL software library (LAPACK, specifically).  
It is possible to use other compilers (e.g. gfortran) provided OpenMPI and LAPACK are also available.  

LPPM uses interpolation software provided by ACM TOMS algorithm 526 written by Hiroshi Akima and implemented by John Burkhardt (<https://orion.math.iastate.edu/burkardt/f_src/bivar/bivar.html>),
and ACM TOMS algorithms 772 and 773 written by Robert Renka.  All other code is original and authored by Peter Bosler.
	
	H. Akima, 1978, A method of bivariate interpolation and smooth surface fitting for irregularly distributed data points.
		ACM TOMS, 4.
	H. Akima, 1984, On estimating partial derivatives for bivariate interpolation of scattered data.
		Rocky Mountain Journal of Mathematics, 14:41-52.
	R. Renka, 1997, Algorithm 772 : STRIPACK : Delaunay triangulation and Voronoi diagram on the surface of a sphere.
		ACM TOMS, 23:416-434.
	R. Renka, 1997, Algorithm 773 : SSRFPACK : Interpolation of scattered data on the surface of a sphere with a surface under thension.
		ACM TOMS, 23:435-442.
	


Plotting and graphics
=========

LPPM includes subroutines for outputting data into formats readable by the VTK C++ graphics library and ParaView, available 
from www.kitware.org, as well as routines for the NCAR Command Language (NCL) <http://www.ncl.ucar.edu> and Matlab.  



