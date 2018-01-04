#!/bin/bash

CWD=${PWD}
N=( 8 16 32 64 )
NP=8

for i in "${N[@]}"
do
	cp ./module_param${i}.bckp ./module_param.f90 && \
		make clean && make && cp ./incompact3d DATA/ && \
		cd DATA && \
		mpiexec -np ${NP} ./incompact3d | tee OUTPUT${i}.log
	cd ${CWD}
done
