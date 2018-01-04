#!/usr/bin/env python

import math
import matplotlib.pyplot as plt

N = [8,
		 16,
		 32,
		 64
] # Mesh size definition
T = [200,
		 400,
		 800,
		 1600,
] # Timestep definition

DATADIR = "./DATA/"
DATATYPE = ".dat"
DATAFILE = "err_"

VARS = ["rho",
				"u",
				"v",
				#"w"
]

LINES = ["-",
				 "--",
				 "-.",
				 ":"
]

MARKERS = ["",
					 "o",
					 "s",
					 "^"
]

def main():
	""" Main function. """

	errs = {} # Error dictionary
	for var in VARS:
		errs[var] = {}
		for n in N:
			errs[var][n] = {}
			for t in T:
				errs[var][n][t] = 0.0

				# Read error data
				with open(DATADIR + DATAFILE + var + str(n) + "-" + str(t) + DATATYPE, "r") as errfile:
					for row in errfile:
						errs[var][n][t] += float(row.split()[4])

				# Compute avg error
				errs[var][n][t] /= float(t)
				print var + str(n) + "-" + str(t) + ": " + str(errs[var][n][t])

	# Log errors
	with open (DATADIR + "err.log", "w") as logfile:
		for var in VARS:
			for n in N:
				for t in T:
					logfile.write(var + str(n) + "-" + str(t) + ": " + str(errs[var][n][t]))

	# Compute orders of convergence
	pn = {}
	pt = {}
	for var in VARS:
		pn[var] = {}
		pt[var] = {}
		nctr = 0
		for n in N:
			pn[var][n] = {}
			pt[var][n] = {}
			tctr = 0
			for t in T:
				pn[var][n][t] = 0
				pt[var][n][t] = 0

				if nctr:
					pn[var][n][t] = math.log(errs[var][N[nctr]][T[tctr]] \
																	 / errs[var][N[nctr - 1]][T[tctr]]) \
																	 / math.log(N[nctr - 1] / float(N[nctr]))
				if tctr:
					pt[var][n][t] = math.log(errs[var][N[nctr]][T[tctr]] \
																	 / errs[var][N[nctr]][T[tctr - 1]]) \
																	 / math.log(T[tctr] / float(T[tctr - 1]))
				tctr += 1
			nctr += 1

		for n in N:
			for t in T:
				print var + str(n) + "-" + str(t) + ": " + str(pt[var][n][t]) + " (dt)"
		for t in T:
			for n in N:
				print var + str(n) + "-" + str(t) + ": " + str(pn[var][n][t]) + " (dx)"

	# Error plots
	pltvar = {}
	varctr = 0
	plt.figure(figsize=(5.0, 3.5))
	for var in VARS:
		pltvar[var] = {}
		for t in T:
			pltvar[var][t] = []
			for n in N:
				pltvar[var][t].append(errs[var][n][t])
			for n in range(1, len(N) - 1):
				pltvar[var][t][n] /= pltvar[var][t][0]
			pltvar[var][t][0] = 1.0
			label = var + str(t)
			plt.plot(N, pltvar[var][t], label=label,
							 marker=MARKERS[varctr],
							 ls=LINES[varctr]
			)
		varctr += 1
	plt.title("Spatial Convergence")
	plt.xlabel("Ncells")
	plt.ylabel(r"$\varepsilon$")
	plt.xscale("log")
	plt.yscale("log")
	plt.legend()
	plt.show()
	plt.close()

	varctr = 0
	plt.figure(figsize=(5.0, 3.5))
	for var in VARS:
		pltvar[var] = {}
		for n in N:
			pltvar[var][n] = []
			for t in T:
				pltvar[var][n].append(errs[var][n][t])
			for t in range(1, len(T) - 1):
				pltvar[var][n][t] /= pltvar[var][n][0]
			pltvar[var][n][0] = 1.0
			label = var + str(n)
			plt.plot(T, pltvar[var][n], label=label,
							 marker=MARKERS[varctr],
							 ls=LINES[varctr]
			)
		varctr += 1
	plt.title("Temporal Convergence")
	plt.xlabel("Nsteps")
	plt.ylabel(r"$\varepsilon$")
	plt.xscale("log")
	plt.yscale("log")
	plt.legend()
	plt.show()
	plt.close()

if __name__ == "__main__":
	main()
