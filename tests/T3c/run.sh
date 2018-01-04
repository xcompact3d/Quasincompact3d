#!/bin/bash

CWD=${PWD}
N=( 8 16 32 64 )
T=( 200 400 800 1600 )
NP=8

#
# Run the tests
#

# Loop over meshes
for i in "${N[@]}"
do
		# Build and copy executable to case directory
		cp ./module_param${i}.bckp ./module_param.f90 && \
				make clean && make && cp ./incompact3d DATA/ && \

				# Loop over timesteps
				for j in "${T[@]}"
				do
						# Copy runfile to case directory and run case
						cp ./incompact3d${j}.prm DATA/incompact3d.prm && \
								cd DATA && \
								mpiexec -np ${NP} ./incompact3d | tee OUTPUT${i}-${j}.log

						# Generate data for post-processing
						grep Error OUTPUT${i}-${j}.log > err.dat
						
						grep RHO err.dat > err_rho${i}-${j}.dat
						grep UVEL err.dat > err_u${i}-${j}.dat
						grep VVEL err.dat > err_v${i}-${j}.dat
						grep WVEL err.dat > err_w${i}-${j}.dat

						rm err.dat
						
						# Return to build directory
						cd ${CWD}
				done
done

# Run postprocessing
./postproc.py 
