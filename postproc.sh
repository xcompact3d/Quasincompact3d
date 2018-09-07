#! /bin/sh

echo "Updating .dat files"
grep ENSTROPHY OUTPUT.log > ENSTROPHY.dat
grep KINETIC OUTPUT.log > KINETIC.dat
grep "PHI min" OUTPUT.log > PHIMIN.dat
grep "PHI max" OUTPUT.log > PHIMAX.dat
grep "RHO min" OUTPUT.log > RHOMIN.dat
grep "RHO max" OUTPUT.log > RHOMAX.dat

echo "Plotting..."
python plot.py
okular plot.eps 
