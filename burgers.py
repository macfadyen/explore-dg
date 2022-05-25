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
    diff,
    einsum,
    exp,
    linspace,
    maximum,
    ndenumerate,
    sin,
    where,
    zeros,
    zeros_like,
)


class Timer:
    def __enter__(self):
        self.start = default_timer()
        return self

    def __exit__(self, *args):
        self.elapsed_secs = default_timer() - self.start


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
    return 0.5 * u * u


def riemann_solver(ul, ur):
    al = ul
    ar = ur
    ap = maximum(al, ar)
    am = zeros_like(ap)
    fl = flux_function(ul)
    fr = flux_function(ur)
    return (fl * ap - fr * am - (ul - ur) * ap * am) / (ap - am)


def time_derivative(w, x_grid, dx, ck=0.01):
    order = w.shape[-1]
    faces = 0.5 * (x_grid.flat[:-1] + x_grid.flat[+1:])
    x_vol, gauss_weight = leggauss(order)
    x_srf = array([-1.0, 1.0])
    phi_grd = array([leg(x_vol, n, m=1) for n in range(order)])
    phi_vol = array([leg(x_vol, n, m=0) for n in range(order)])
    phi_srf = array([leg(x_srf, n, m=0) for n in range(order)])
    u_vol = einsum("in,nq", w, phi_vol)
    u_srf = einsum("in,nq", w, phi_srf)
    f_vol = flux_function(u_vol)
    f_hat = riemann_solver(u_srf[:-1, 1], u_srf[+1:, 0])
    f_god = riemann_solver(u_vol.flat[:-1], u_vol.flat[+1:])

    u_dot_god = zeros(w.size)
    u_dot_god[+1:-1] -= diff(f_god) / diff(faces)
    u_dot_god.shape = w.shape
    w_dot_god = field_to_weights(u_dot_god)

    w_dot_vol = einsum("iq,nq,q", f_vol, phi_grd, gauss_weight)
    w_dot_srf = zeros_like(w_dot_vol)
    w_dot_srf[+1:-1] += einsum("i,n", f_hat[:-1], phi_srf[:, 0])
    w_dot_srf[+1:-1] -= einsum("i,n", f_hat[+1:], phi_srf[:, 1])
    w_dot = (w_dot_srf + w_dot_vol) / dx

    troubled_zones = where(abs(w[:, -1] / w[:, 0]) > ck)[0]

    for zone in troubled_zones:
        w_dot[zone] = w_dot_god[zone]

    time_derivative.troubled_zones = troubled_zones

    set_bc(w_dot)
    return w_dot


def set_bc(w):
    w[+0] = w[-2]
    w[-1] = w[+1]


def initial_condition(x):
    return exp(-((x - 0.5) ** 2) / 0.01)


def plot_soln(x, w, filename=None):
    import matplotlib.pyplot as plt

    y1 = weights_to_field(w)
    plt.figure(figsize=[11, 8.5])
    for i in range(x.shape[0]):
        plt.plot(x[i, :], y1[i, :], "-o", c="rk"[i % 2], mfc="none")

    try:
        for zone in time_derivative.troubled_zones:
            plt.plot(x[zone], zeros_like(x[zone]), c="k", lw=8.0)
    except AttributeError:
        pass

    plt.ylim(0.0, 1.2)
    plt.ylabel(r"$u$")
    plt.xlabel(r"$x$")
    plt.tight_layout(pad=0.05)
    if filename is None:
        plt.show()
    else:
        # plt.xlim(0.635, 0.665)
        plt.savefig(filename)
    plt.close()


def main():
    import argparse

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
    parser.add_argument(
        "--rk-order",
        choices=[1, 2, 3, 4],
        default=2,
        type=int,
    )
    parser.add_argument(
        "--cfl",
        type=float,
        default=0.05,
    )
    parser.add_argument(
        "--limiter-param",
        "-c",
        type=float,
        default=0.01,
    )
    parser.add_argument(
        "--plot",
        type=str,
        default="end",
        choices=["end", "frames", "none"],
    )
    args = parser.parse_args()

    dx, grid = gauss_grid(0.0, 1.0, args.zones, order=args.order, retstep=True)
    dt = dx / 1.0 * args.cfl
    x = grid
    y = initial_condition(x)
    w = field_to_weights(y)
    t = 0.0
    n = 0
    L = lambda w: time_derivative(w, grid, dx=dx, ck=args.limiter_param)

    while t < args.tfinal:
        with Timer() as target:

            if args.rk_order == 1:
                w += dt * L(w)

            if args.rk_order == 2:
                w0 = w
                w1 = w0 + dt * L(w0)
                w2 = 1 / 2 * w0 + 1 / 2 * (w1 + dt * L(w1))
                w = w2

            if args.rk_order == 3:
                w0 = w
                w1 = w0 + dt * L(w0)
                w2 = 3 / 4 * w0 + 1 / 4 * (w1 + dt * L(w1))
                w3 = 1 / 3 * w0 + 2 / 3 * (w2 + dt * L(w2))
                w = w3

            if args.rk_order == 4:
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

        if args.plot == "frames":
            plot_soln(x, w, filename=f"frame.{n:04d}.png")

    if args.plot == "end":
        plot_soln(x, w)


if __name__ == "__main__":
    main()
