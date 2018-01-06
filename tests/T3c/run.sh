#!/bin/bash

CWD=${PWD}
N=( 8 16 32 64 128 )
T=( 200 400 800 1600 )

# Mesh
NX=8
NY=8
NZ=8

# Processor decomposition
PROW=4
PCOL=2
NP=$((${PROW} * ${PCOL}))

echo "Running on ${NP} processors in ${PROW} X ${PCOL} arrangement"

VALGR=""
#VALGR="valgrind --tool=memcheck --log-file=memcheck%p.log --track-origins=yes"

#
# Run the tests
#

# Loop over meshes
for i in "${N[@]}"
do
		# Build and copy executable to case directory
		cp ./module_param${i}.bckp ./module_param.f90 && \
				sed -i -e "s/p_row=.*,/p_row=${PROW},/g" ./module_param.f90 && \
				sed -i -e "s/p_col=.*/p_col=${PCOL}/g" ./module_param.f90 && \
				make clean && make && cp ./incompact3d DATA/ && \
				
				# Loop over timesteps
				for j in "${T[@]}"
				do
						echo "Running mesh N=${i} for T=${j}"
						
						# Copy runfile to case directory and run case
						cp ./incompact3d${j}.prm DATA/incompact3d.prm && \
								cd DATA && \
								mpiexec -np ${NP} ${VALGR} ./incompact3d > OUTPUT${i}-${j}.log

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
