LPPM : Lagrangian Particle / Panel Method
=========

Code associated with 

P. Bosler,  L. Wang,  R. Krasny, and C. Jablonowski,  2014.  
	"A Lagrangian particle/panel method for the barotropic vorticity equation on a rotating sphere," Fluid Dynamics Research,  46 : 031406.

Bosler, P. A., "Particle Methods for Geophysical Flow on the Sphere," PhD Thesis, University of Michigan, 2013.

This software provides a Lagrangian Particle-Panel Method for solving the barotropic vorticity equation
and the advection equation on the sphere and in the plane.

A shallow water equation solver is in development.

Requirements
=========
The code requires a modern Fortran compiler with an OpenMPI distribution for parallel computing.  

Currently the software is written for Intel compilers and makes use of the Intel MKL software library (LAPACK, specifically).  
It is possible to use other compilers (e.g. gfortran) provided OpenMPI and LAPACK are also available.  

LPPM uses interpolation software provided by ACM TOMS algorithm 526 written by Hiroshi Akima and implemented by John Burkhardt,
and ACM TOMS algorithms 772 and 773 written by Robert Renka.  All other code is authored by Peter Bosler.


Plotting and graphics
=========

LPPM includes subroutines for outputting data into formats readable by the VTK C++ graphics library and ParaView, available 
from www.kitware.org, as well as routines for the NCAR Command Language (NCL) and Matlab.  


BigSphere
