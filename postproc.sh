#! /bin/sh

echo "Updating .dat files"
grep ENSTROPHY OUTPUT.log > ENSTROPHY.dat
grep KINETIC OUTPUT.log > KINETIC.dat

echo "Plotting..."
python plot.py
okular plot.eps 
