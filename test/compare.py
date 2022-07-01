#! /bin/python
import sys

file1 = []
file2 = []
with open(sys.argv[1]) as fp:
    line = fp.readline()
    while line:
        line2 = line.rstrip(" \n")
        tmp = line2.split(" ")
        tmp2 = [float(x) for x in tmp]
        file1.append(tmp2)
        line = fp.readline()

with open(sys.argv[1]) as fp:
    line = fp.readline()
    while line:
        line2 = line.rstrip(" \n")
        tmp = line2.split(" ")
        tmp2 = [float(x) for x in tmp]
        file2.append(tmp2)
        line = fp.readline()




import numpy as np
file1 = np.genfromtxt(sys.argv[1], delimiter=' ', dtype=np.double)
file2 = np.genfromtxt(sys.argv[2], delimiter=' ', dtype=np.double)

diff = np.abs(file1 - file2)
equal = not np.any(diff > 1e-14)
if equal:
    sys.exit(0)
else:
    sys.exit(1)
