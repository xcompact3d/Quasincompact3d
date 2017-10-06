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

	plt.subplot(2, 1, 1)
	plt.plot(t, k, color="black")
	plt.xlabel("t")
	plt.ylabel("k")
	plt.xlim(xmax=max(t))

	t = []
	enstrophy = []
	with open("ENSTROPHY.dat", "r") as data:
		for row in data:
			words = row.split()
			t.append(float(words[1]))
			enstrophy.append(float(words[2]))

	plt.subplot(2, 1, 2)
	plt.plot(t, enstrophy, color="black")
	plt.xlabel("t")
	plt.ylabel("Enstrophy")
	plt.xlim(xmax=max(t))
	# plt.show()
	plt.savefig("plot.eps")

if __name__ == "__main__":
	main()
