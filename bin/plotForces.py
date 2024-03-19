#!/usr/bin/python3
import matplotlib.pyplot as plt
import sys

forcesfilename = "forces/0/forces.dat"
moment = False
if len(sys.argv) != 1:
    if sys.argv[1] == "moment":
        moment = True

results = [[] for _ in range(13)]
with open(forcesfilename, "r") as file:
    for line in file:
        if not line.startswith("#"):
            numbers = line.replace("(", "").replace(")", "").split()
            [results[i].append(float(num)) for i, num in enumerate(numbers)]

# combine pressure and viscous forces and moments
t = results[0]
fx = [fxp + fxv for fxp, fxv in zip(results[1], results[4])]
fy = [fyp + fyv for fyp, fyv in zip(results[2], results[5])]
fz = [fzp + fzv for fzp, fzv in zip(results[3], results[6])]
mx = [mxp + mxv for mxp, mxv in zip(results[7], results[10])]
my = [myp + myv for myp, myv in zip(results[8], results[11])]
mz = [mzp + mzv for mzp, mzv in zip(results[9], results[12])]

# write clean data file
with open("forces.dat", "w") as outForces:
    for data in zip(t, fx, fy, fz):
        outForces.write(" ".join([str(d) for d in data]) + "\n")

with open("moments.dat", "w") as outMoments:
    for data in zip(t, mx, my, mz):
        outMoments.write(" ".join([str(d) for d in data]) + "\n")

if moment:
    plt.plot(t, mx, label="mx")
    plt.plot(t, my, label="my")
    plt.plot(t, mz, label="mz")
    plt.ylabel("moment (Nm)")
else:
    plt.plot(t, fx, label="fx")
    plt.plot(t, fy, label="fy")
    plt.plot(t, fz, label="fz")
    plt.ylabel("force (N)")

plt.xlabel("iteration")
plt.grid(True)
plt.legend()
plt.show()
