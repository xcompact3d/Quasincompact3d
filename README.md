# README #

This is the repository for the Low Mach Number (LMN) version of incompact3d.

* Quick summary
* Road map
	* v1.0 Impose divergence of momentum
	* v1.1 Solve for divergence of velocity using extrapolated pressure field
	* v1.2 Solve for divergence of velocity using defect correction
* Version: v0.0.T3a

### How do I get set up? ###

* To create a case edit the initialisation/boundary routines in navier.f90 and set mesh size in module_param.f90. Other values (e.g. dt) can be set at runtime by editting incompact3d.prm
* To build: make clean && make
* Dependencies:
	* A Fortran90 compatible compiler (gfortran, ifort, ftn)
	* MPI wrappers: OpenMPI/vendorMPI are recommended (people have reported problems using MPICH)
* To run: mpiexec -np n ./incompact3d


### Who do I talk to? ###

* Paul Bartholomew <paul.bartholomew08@imperial.ac.uk>