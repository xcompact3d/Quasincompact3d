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

	plt.plot(t, k)
	plt.xlabel("t")
	plt.ylabel("k")
	plt.xlim(xmax=max(t))
	plt.show()

	t = []
	enstrophy = []
	with open("ENSTROPHY.dat", "r") as data:
		for row in data:
			words = row.split()
			t.append(float(words[1]))
			enstrophy.append(float(words[2]))

	plt.plot(t, enstrophy)
	plt.xlabel("t")
	plt.ylabel("Enstrophy")
	plt.xlim(xmax=max(t))
	plt.show()

if __name__ == "__main__":
	main()
