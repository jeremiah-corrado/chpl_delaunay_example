
import sys
import numpy as np
import matplotlib.pyplot as plt
import math

npoints = int(sys.argv[1])
x = np.zeros(npoints)
y = np.zeros(npoints)

with open("results/points.txt") as f:
    for (i, line) in enumerate(f):
        xy = line.split(" ")
        x[i] = xy[0]
        y[i] = xy[1]

plt.scatter(x, y)

with open("results/edges.txt") as f:
    for line in f.readlines():
        edge = line.split(" ")
        aID, bID = int(edge[0]), int(edge[1])
        plt.plot([x[aID], x[bID]], [y[aID], y[bID]], 'r-')

plt.savefig("tri.png")
