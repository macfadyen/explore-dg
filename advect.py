"""
Solves the advection equation, nth order (DG) in space, 4th order in time

Author: Jonathan Zrake
"""

from timeit import default_timer
from numpy.polynomial.legendre import Legendre, leggauss
from numpy import (
    allclose,
    array,
    cos,
    einsum,
    exp,
    linspace,
    ndenumerate,
    sin,
    zeros,
    zeros_like,
)


class Timer:
    def __enter__(self):
        self.start = default_timer()
        return self

    def __exit__(self, *args):
        self.elapsed_secs = default_timer() - self.start


ADVECT_WAVESPEED = 1.0
CFL_NUMBER = 0.02


def leg(x, n, m=0):
    """
    Return the value of the scaled Legendre polynomial.

    The point x is in the range [-1, 1]. m is the derivative order.
    """
    c = [(2 * n + 1) ** 0.5 if i is n else 0.0 for i in range(n + 1)]
    return Legendre(c).deriv(m)(x)


def gauss_grid(x0, x1, num_zones, order, retstep=False):
    """
    Return a 2d array of Gauss quadrature points.

    The first axis is over zones, and the second axis is over Gauss points
    within the zone.
    """
    xsi, _ = leggauss(order)
    dx = (x1 - x0) / num_zones

    x = zeros([num_zones, order])

    for i in range(num_zones):
        x[i] = [x0 + dx * (i + 0.5 * (1.0 + y)) for y in xsi]

    if retstep:
        return dx, x
    else:
        return x


def field_to_weights(field):
    order = field.shape[-1]
    x, w = leggauss(order)
    phi = array([leg(x, n) for n in range(order)])
    return einsum("iq,nq,q,...", field, phi, w, 0.5)


def weights_to_field(weights):
    order = weights.shape[-1]
    x, _ = leggauss(order)
    phi = array([leg(x, n) for n in range(order)])
    return einsum("in,nq", weights, phi)


def flux_function(u):
    return u * ADVECT_WAVESPEED


def riemann_solver(ul, ur):
    if ADVECT_WAVESPEED < 0.0:
        return ur * ADVECT_WAVESPEED
    if ADVECT_WAVESPEED > 0.0:
        return ul * ADVECT_WAVESPEED


def time_derivative(w, dx):
    order = w.shape[-1]
    x_vol, gauss_weight = leggauss(order)
    x_srf = array([-1.0, 1.0])
    phi_grd = array([leg(x_vol, n, m=1) for n in range(order)])
    phi_vol = array([leg(x_vol, n, m=0) for n in range(order)])
    phi_srf = array([leg(x_srf, n, m=0) for n in range(order)])
    u_vol = einsum("in,nq", w, phi_vol)
    u_srf = einsum("in,nq", w, phi_srf)
    f_vol = flux_function(u_vol)
    f_hat = riemann_solver(u_srf[:-1, 1], u_srf[+1:, 0])
    w_dot_vol = einsum("iq,nq,q", f_vol, phi_grd, gauss_weight)
    w_dot_srf = zeros_like(w_dot_vol)
    w_dot_srf[+1:-1] += einsum("i,n", f_hat[:-1], phi_srf[:, 0])
    w_dot_srf[+1:-1] -= einsum("i,n", f_hat[+1:], phi_srf[:, 1])
    res = (w_dot_srf + w_dot_vol) / dx
    set_bc(res)
    return res


def set_bc(w):
    w[+0] = w[-2]
    w[-1] = w[+1]


def initial_condition(x, t=0.0, x0=0.0, x1=1.0):
    x = x - ADVECT_WAVESPEED * t
    for i, _ in ndenumerate(x):
        while x[i] < x0:
            x[i] += x1 - x0
        while x[i] > x1:
            x[i] -= x1 - x0
    return exp(-((x - 0.5) ** 2) / 0.01)


def main():
    import argparse
    import matplotlib.pyplot as plt

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--order",
        default=7,
        type=int,
        help="polynomial order to use",
    )
    parser.add_argument(
        "--zones",
        default=5,
        type=int,
        help="number of zones on the grid",
    )
    parser.add_argument(
        "--tfinal",
        "-e",
        default=0.1,
        type=float,
    )
    args = parser.parse_args()

    dx, grid = gauss_grid(0.0, 1.0, args.zones, order=args.order, retstep=True)
    dt = dx / abs(ADVECT_WAVESPEED) * CFL_NUMBER
    x_fine = linspace(0.0, 1.0, 2000)
    y_fine = initial_condition(x_fine)

    x = grid
    y = initial_condition(x)
    w = field_to_weights(y)
    t = 0.0
    n = 0
    L = lambda w: time_derivative(w, dx=dx)

    while t < args.tfinal:
        with Timer() as target:
            k1 = dt * L(w)
            k2 = dt * L(w + 0.5 * k1)
            k3 = dt * L(w + 0.5 * k2)
            k4 = dt * L(w + 1.0 * k3)
            w += (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0
            t += dt
            n += 1

        zps = args.zones / target.elapsed_secs
        if n % 1 == 0:
            print(f"[{n:04d}]: t={t:.3f} zps={zps:.3e}")

    y1 = weights_to_field(w)

    y1_true = initial_condition(x_fine, t=t, x0=dx, x1=1.0 - dx)
    y1_comp = initial_condition(x, t=t, x0=dx, x1=1.0 - dx)
    print(f"L2 error = {sum(((y1 - y1_comp) ** 2).flat) ** 0.5:.3e}")

    plt.plot(x_fine, y1_true, lw=0.5, ls="--", c="b", label=r"$f(x, t)$")

    for i in range(x.shape[0]):
        plt.plot(x[i, :], y1[i, :], "o", c="rk"[i % 2], mfc="none")

    plt.legend()
    plt.show()


if __name__ == "__main__":
    main()
