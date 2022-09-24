#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import os
import math

order = 3
cfl_parameter = 0.4 / (2.0 * (order - 1.0) + 1)
test_time = 0.01
wave_velocity = 1.0
resolutions = [64, 128, 256, 512, 1024, 2048]
l1s = []


def compute_l1(x, prim, time):
    from numpy.polynomial.legendre import leggauss

    xg, w = leggauss(order)
    rho0 = np.zeros(len(prim))
    drho = np.zeros(len(prim))

    for i in range(order, (len(x) - order)):
        rho0[i] = 1.0 + 0.5 * np.sin(2.0 * math.pi * (x[i] - wave_velocity * time))
        drho[i] = np.abs(prim[i, 0] - rho0[i])

    nx = int(len(prim) / order)
    l1 = 0.0
    for i in range(1, nx - 1):
        for r in range(order):
            l1 += drho[i * order + r] * w[r]

    return l1 / float(nx)


for i in range(len(resolutions)):
    run_command = (
        "./euler num_zones="
        + str(resolutions[i])
        + " order="
        + str(order)
        + " cfl="
        + str(cfl_parameter)
        + " tmax="
        + str(test_time)
        + " gamma=1.4 rk=3 bc_type=1 run terminal=grid.dat grid:print terminal=data.dat prim:print"
    )
    os.system(run_command)
    f = open("grid.dat", "r")
    x = np.genfromtxt(f)
    f.close()
    f = open("data.dat", "r")
    prim = np.genfromtxt(f)
    f.close()
    l1i = compute_l1(x, prim, test_time)
    l1s.append(l1i)

for i in range(len(resolutions)):
    print(resolutions[i], l1s[i])

fig, ax = plt.subplots()
ax.set_yscale("log")
# ax.set(xlim=(30, 1100))
ax.set_title(r"Density Wave Convergence - Order " + str(order))
ax.set_xlabel(r"N")
ax.set_ylabel(r"L1 ($\rho$)")
ax.loglog(resolutions, l1s, "o")
ax.plot(
    [resolutions[0], resolutions[-1]],
    [l1s[0], (resolutions[0] / resolutions[-1]) ** order * l1s[0]],
    color="b",
    linestyle="dashed",
    linewidth=2,
)
plt.show()
