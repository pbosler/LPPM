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
	
	--- run the problem --
	
	call Delete(rhWave)

Requirements
=========
The code requires a modern Fortran compiler with an OpenMPI distribution for parallel computing.  
Without OpenMPI serial computations are still possible, but the code is designed to use distributed memory parallelism.

Currently the software is written for Intel compilers and makes use of the Intel MKL software library (LAPACK, specifically).  
It is possible to use other compilers (e.g. gfortran) provided OpenMPI and LAPACK are also available.  

LPPM uses interpolation software provided by ACM TOMS algorithm 526 written by Hiroshi Akima and implemented by John Burkhardt,
and ACM TOMS algorithms 772 and 773 written by Robert Renka.  All other code is authored by Peter Bosler.


Plotting and graphics
=========

LPPM includes subroutines for outputting data into formats readable by the VTK C++ graphics library and ParaView, available 
from www.kitware.org, as well as routines for the NCAR Command Language (NCL) and Matlab.  



