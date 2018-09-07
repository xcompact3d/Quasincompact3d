"""
.. module:: plot.py
    :synopsis: Plots the kinetic energy in file KINETIC.dat
.. moduleauthor:: Paul Bartholomew
"""

import matplotlib.pyplot as plt

def main():
	"""Main function, does the plotting."""

	t = []
	k = []
	with open("KINETIC.dat", "r") as data:
		for row in data:
			words = row.split()
			t.append(float(words[2]))
			k.append(float(words[3]))

	plt.subplot(2, 2, 1)
	plt.plot(t, k, color="black")
	plt.xlabel("t")
	plt.ylabel(r"$k$")
	plt.xlim(xmax=max(t))

	t = []
	enstrophy = []
	with open("ENSTROPHY.dat", "r") as data:
		for row in data:
			words = row.split()
			t.append(float(words[1]))
			enstrophy.append(float(words[2]))

	plt.subplot(2, 2, 3)
	plt.plot(t, enstrophy, color="black")
	plt.xlabel("t")
	plt.ylabel("Enstrophy")
	plt.xlim(xmax=max(t))

	# Scalar error
	phimin = []
	with open("PHIMIN.dat", "r") as data:
		for row in data:
			words = row.split()
			phimin.append(float(words[2]))
	for i in range(len(phimin)):
		phimin[i] = 100.0 * ((1.0 - phimin[i]) / 1.0)
	plt.subplot(2, 2, 2)
	plt.plot(phimin, ls="--", color="black")
	
	phimax = []
	with open("PHIMAX.dat", "r") as data:
		for row in data:
			words = row.split()
			phimax.append(float(words[2]))
	for i in range(len(phimax)):
		phimax[i] = 100.0 * ((2.0 - phimax[i]) / 2.0)
	plt.plot(phimax, color="black")
	plt.ylabel(r"$\varepsilon_{\phi}$ [%]")

	# Density error
	rhomin = []
	with open("RHOMIN.dat", "r") as data:
		for row in data:
			words = row.split()
			rhomin.append(float(words[2]))
	for i in range(len(rhomin)):
		rhomin[i] = 100.0 * ((1.0 - rhomin[i]) / 1.0)
	plt.subplot(2, 2, 4)
	plt.plot(rhomin, ls="--", color="black")
	
	rhomax = []
	with open("RHOMAX.dat", "r") as data:
		for row in data:
			words = row.split()
			rhomax.append(float(words[2]))
	for i in range(len(rhomax)):
		rhomax[i] = 100.0 * ((2.0 - rhomax[i]) / 2.0)
	plt.plot(rhomax, color="black")
	plt.ylabel(r"$\varepsilon_{\rho}$ [%]")
	
	# plt.show()
	plt.savefig("plot.eps")

if __name__ == "__main__":
	main()
