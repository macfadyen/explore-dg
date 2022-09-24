#!/usr/bin/env python3
import numpy as np
import math
import matplotlib.pyplot as plt
import os
from numpy.polynomial.legendre import leggauss

from euler_solvers import euler_exact_prim

NUM_FIELDS = 3

test = 2
order = 2
cfl_parameter = 0.4 / (2.0 * (order - 1.0) + 1)
tci_method = 3
tci = 0.0
use_analytic = False
floor = 0.0
tmax = 0.1

if test == 1:
    gamma = 5.0 / 3.0
    tmax = 8.0e-3
    resolution = 100
    bc = 0
    use_analytic = False
elif test == 2:
    gamma = 5.0 / 3.0
    tmax = 1e-1
    resolution = 100
    bc = 0
    use_analytic = False

tmax_an = 10.0 * tmax

xmin = 0.0
xmax = 1.0
dx = (xmax - xmin) / resolution

run_command = (
    "./srhd num_zones="
    + str(resolution)
    + " order="
    + str(order)
    + " tci="
    + str(tci)
    + " cfl="
    + str(cfl_parameter)
    + " gamma="
    + str(gamma)
    + " tmax="
    + str(tmax)
    + " init="
    + str(test)
    + " floor="
    + str(floor)
    + " tci_type="
    + str(tci_method)
    + " bc_type="
    + str(bc)
    + " solver_type=1 rk=1 run terminal=grid.dat grid:print terminal=prim.dat prim:print terminal=wgts.dat wgts:print"
)
print("run_command = ", run_command)
os.system(run_command)
f = open("grid.dat", "r")

xn = np.genfromtxt(f)
f.close()
f = open("prim.dat", "r")
prim = np.genfromtxt(f)
f.close()
f = open("wgts.dat", "r")
wgts = np.genfromtxt(f)
f.close()

# trim the ghost zones
xn = xn[order : (len(xn) - order)]
prim = prim[order : (np.shape(prim)[0] - order), :]
wgts = wgts[NUM_FIELDS : (np.shape(wgts)[0] - NUM_FIELDS), :]

# number of active node points (excluding ghost zones on each end)
nn = len(xn)
# number of zones (excluding ghost zones)
nx = int(nn / order)  # = order -2

# analytic solutions


def burgers_soln(x, t):
    import math
    from scipy import optimize

    """
    Analytic solution from initial condition from Burgers setup in
    sailfish/setups/simple1d.py
    u(x,t) > 0 for this setup
    Use root finder to find xsi such that xsi - x + f(xsi) * t = 0
    """

    a = 0.2
    k = 2.0 * math.pi
    average_wavespeed = 1.0
    # solution for the Burgers eqn in Example 3.1 includes this factor
    z = 2.0 * math.sqrt(3.0)

    def f(xsi):
        return xsi - x + t * (1.0 + a * math.sin(k * xsi)) * z

    def fder(xsi):
        return 1.0 + t * k * a * math.cos(k * xsi) * z

    def fder2(xsi):
        return -t * k * k * a * math.sin(k * xsi) * z

    # xsi0 is initial guess for xsi
    xsi0 = x - average_wavespeed * t * z
    xsi = optimize.newton(f, xsi0, fprime=fder, fprime2=fder2)
    return 1.0 + a * math.sin(k * xsi)


def prim_l1(x, prim, t):
    from numpy.polynomial.legendre import leggauss

    nx = int(len(xn) / order)

    xg, w = leggauss(order)
    drho = np.zeros(len(x))

    for i in range(len(x)):
        if test == 1:
            rho_an = burgers_soln(x[i], t)
        else:
            s = (x[i] - x0) / t
            rho_an = euler_exact_prim(priml, primr, s)[0]
        drho[i] = abs(prim[i, 0] - rho_an)
        # print("i, x[i], t, rho_an, drho[i] = ", i, x[i], t, rho_an, drho[i])

    l1 = 0.0
    for iz in range(nx):
        for r in range(order):
            l1 += drho[iz * order + r] * w[r]

    return 0.5 * l1 / nx


l1 = 0.0
if use_analytic:
    xan = np.linspace(xmin, xmax, num=4000)
    rho_an = np.zeros_like(xan)

    for i in range(len(xan)):
        if test == 1:
            rho_an[i] = burgers_soln(xan[i], tmax)
        elif test > 1 and test < 6:
            s = (xan[i] - x0) / tmax
            rho_an[i] = euler_exact_prim(priml, primr, s)[0]

    l1 = prim_l1(xn, prim, tmax)

rho_zone = np.zeros(resolution)
xz = xmin + dx * (np.arange(resolution) + 0.5)

iz = 0
for i in range(np.shape(wgts)[0]):
    if (i % NUM_FIELDS) == 0:
        rho_zone[iz] = wgts[i, 0]
        iz += 1

fig = plt.figure(1)

# fig, axs = plt.subplots(2,1)

ax1 = fig.add_subplot(211)
ax1.plot(xz, rho_zone, "ko", mfc="none")
if use_analytic or (test == 6 or test == 7):
    ax1.plot(xan, rho_an, "k-")
# ax1.plot(xn, prim[:, 2] / prim[:, 0])
# plt.xlabel(r"$x$")
plt.ylabel(r"$\rho$")

ax1.set_xlim([xmin, xmax])
# ax1.xlabel(r"$x$")
# ax1.ylabel(r"$\rho$")
# plt.show()

f = open("trzn.dat", "r")
trzn = np.genfromtxt(f)
f.close()
# nxt = np.shape(trzn)[1]
# trzn = trzn[:, 1 : (nxt - 1)]  # trim ghost zones

# ax2 = fig.add_subplot(211)

nt, nx = np.shape(trzn)
ntt = np.count_nonzero(trzn > tci)
print("number of troubled zones ntt = ", ntt)
pt = np.count_nonzero(trzn > tci) / np.size(trzn)
print("average percent troubled = ", 100 * pt)
mpt = max(np.count_nonzero(trzn > tci, axis=1)) / np.shape(trzn)[1]
print("maximum percent troubled = ", 100 * mpt)
a = np.where(trzn > tci, 1, 0)
ax2 = fig.add_subplot(212)
plt.xlabel(r"$x$")
# ax2.imshow(np.log10(np.flipud(trzn)), aspect="auto")
ax2.pcolor(a, cmap="binary")
fig.suptitle(
    "Test={} Order={} TCI_type={} TCI={:2.2} CFL={:2.2}\n L1= {:2.2} Trouble: ave= {:2.2%} max= {:2.2%}".format(
        test, order, tci_method, tci, cfl_parameter, l1, pt, mpt
    )
)
# plt.show()
plt.savefig("FS3.{}_{}_TCI{}_N{}.png".format(test, order, tci_method, resolution))

plt.show()
os.system("rm trzn.dat")
