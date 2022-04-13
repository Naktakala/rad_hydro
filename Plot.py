import sys
import numpy as np
import matplotlib.pyplot as plt

num_args = len(sys.argv)
if num_args == 1:
    sys.exit(0)

n = int(sys.argv[1])

ofile = open("Output.txt", "r")

u=np.zeros(100)
lines = ofile.readlines()
for line in lines:
    words = line.split()
    if len(words) == 0:
        continue

    if words[0] == "Timestep:":
        if int(words[1]) == n:
            for w in range(0,100):
                u[w] = float(words[5+w])


ofile.close()

plt.figure()
plt.plot(u,"k-o", markersize=3)
plt.savefig("Plot.png")

print("Done")
