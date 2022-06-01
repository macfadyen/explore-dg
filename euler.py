"""
Solves the advection equation, nth order (DG) in space, 4th order in time

Author: Jonathan Zrake
"""

from timeit import default_timer
from numpy.polynomial.legendre import Legendre, leggauss
from numpy import (
    allclose,
    apply_along_axis,
    array,
    cos,
    diff,
    einsum,
    exp,
    linspace,
    maximum,
    ndenumerate,
    ones_like,
    pi,
    sin,
    transpose,
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


ADIABATIC_GAMMA = 5 / 3


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
    """
    q: field
    i, j, k: zone
    l, m, n: poly
    r, s, t: gauss point

    weights[i, q, l] = field[i, q, r] * phi[l, r] * w[r] * 0.5
    """
    order = field.shape[-1]
    x, w = leggauss(order)
    phi = array([leg(x, n) for n in range(order)])
    return einsum("iqr,lr,r,...->iql", field, phi, w, 0.5)


def weights_to_field(weights):
    """
    q: field
    i, j, k: zone
    l, m, n: poly
    r, s, t: gauss point

    field[i, q, r] = weights[i, q, l] * phi[l, r]
    """
    order = weights.shape[-1]
    x, _ = leggauss(order)
    phi = array([leg(x, n) for n in range(order)])
    return einsum("iql,lr->iqr", weights, phi)


def prim_to_cons(p):
    rho = p[0]
    vel = p[1]
    pre = p[2]
    return rho, rho * vel, 0.5 * rho * vel**2 + pre / (ADIABATIC_GAMMA - 1.0)


def cons_to_prim(u):
    d = u[0]
    s = u[1]
    E = u[2]
    pre = (E - 0.5 * s**2 / d) * (ADIABATIC_GAMMA - 1.0)
    if d <= 0.0:
        raise RuntimeError("negative density")
    if pre <= 0.0:
        raise RuntimeError("negative pressure")
    return d, s / d, pre


def sound_speed(p):
    rho, vel, pre = p
    return (ADIABATIC_GAMMA * pre / rho) ** 0.5


def flux_function(p):
    rho = p[0]
    vel = p[1]
    pre = p[2]
    return (
        rho * vel,
        rho * vel**2 + pre,
        (0.5 * rho * vel**2 + pre / (ADIABATIC_GAMMA - 1) + pre) * vel,
    )


def riemann(ul, ur):
    rho_l, vel_l, pre_l = pl = cons_to_prim(ul)
    rho_r, vel_r, pre_r = pr = cons_to_prim(ur)
    fl = flux_function(pl)
    fr = flux_function(pr)
    cs_l = sound_speed(pl)
    cs_r = sound_speed(pr)
    lam_pl = vel_l + cs_l
    lam_ml = vel_l - cs_l
    lam_pr = vel_r + cs_r
    lam_mr = vel_r - cs_r
    ap = max(0.0, lam_pl, lam_pr)
    am = min(0.0, lam_ml, lam_mr)
    return [
        (ap * fl[i] - am * fr[i] + ap * am * (ur[i] - ul[i])) / (ap - am)
        for i in range(3)
    ]


def time_derivative(w, x_grid, dx, setup, ck=0.01):
    order = w.shape[-1]
    faces = 0.5 * (x_grid.flat[:-1] + x_grid.flat[+1:])
    x_vol, gauss_weight = leggauss(order)
    x_srf = array([-1.0, 1.0])
    phi_grd = array([leg(x_vol, n, m=1) for n in range(order)])
    phi_vol = array([leg(x_vol, n, m=0) for n in range(order)])
    phi_srf = array([leg(x_srf, n, m=0) for n in range(order)])

    u_vol = einsum("iql,lr->iqr", w, phi_vol)
    u_srf = einsum("iql,lr->iqr", w, phi_srf)
    p_vol = apply_along_axis(cons_to_prim, axis=1, arr=u_vol)
    f_vol = apply_along_axis(flux_function, axis=1, arr=p_vol)

    ul, ur = u_srf[:-1, :, 1], u_srf[+1:, :, 0]
    f_hat = array([riemann(ul[i], ur[i]) for i in range(w.shape[0] - 1)])

    w_dot_vol = einsum("iqr,lr,r->iql", f_vol, phi_grd, gauss_weight)
    w_dot_srf = zeros_like(w_dot_vol)
    w_dot_srf[+1:-1] += einsum("iq,l->iql", f_hat[:-1], phi_srf[:, 0])
    w_dot_srf[+1:-1] -= einsum("iq,l->iql", f_hat[+1:], phi_srf[:, 1])
    w_dot = (w_dot_srf + w_dot_vol) / dx

    u = transpose(u_vol, (0, 2, 1)).reshape(-1, w.shape[1])
    ul, ur = u[:-1], u[+1:]
    f_god = array([riemann(ul[i], ur[i]) for i in range(u.shape[0] - 1)])
    u_dot_god = zeros_like(u)
    u_dot_god[order + 1 : -order - 1] -= (
        diff(f_god[order:-order], axis=0) / diff(faces)[:, None]
    )

    u_dot_god.shape = (w.shape[0], order, w.shape[1])
    u_dot_god = u_dot_god.transpose((0, 2, 1))
    w_dot_god = field_to_weights(u_dot_god)

    troubled_zones = where(abs(w[1:-1, 0, -1] / w[1:-1, 0, 0]) > ck)[0]

    # for zone in troubled_zones:
    #     w_dot[zone + 1] = w_dot_god[zone + 1]

    time_derivative.troubled_zones = troubled_zones

    set_bc(w_dot, setup)
    return w_dot


def add_guard_zones(w, setup):
    v = zeros((w.shape[0] + 2,) + w.shape[1:])
    v[1:-1] = w
    set_bc(v, setup)
    return v


def set_bc(w, setup):
    if setup == "density_wave":
        # periodic BC
        w[+0] = w[-2]
        w[-1] = w[+1]
    if setup == "sod":
        # zero-gradient BC
        w[+0] = w[+1]
        w[-1] = w[-2]


def initial_condition(x, setup):
    if setup == "density_wave":
        rho = 1.0 + 0.5 * sin(2 * pi * x)
        vel = 1.0 + zeros_like(x)
        pre = ones_like(x)

    if setup == "sod":
        rho = 1.0 * (x < 0.5) + 0.100 * (x >= 0.5)
        vel = 0.0 * (x < 0.5) + 0.000 * (x >= 0.5)
        pre = 1.0 * (x < 0.5) + 0.125 * (x >= 0.5)

    return transpose(array([rho, vel, pre]), (1, 0, 2))


def plot_soln(x, w, filename=None):
    import matplotlib.pyplot as plt

    y1 = weights_to_field(w)
    plt.figure(figsize=[11, 8.5])

    for i in range(x.shape[0]):
        plt.plot(x[i, :], y1[i, 0, :], "-o", c="rk"[i % 2], mfc="none")

    try:
        for zone in time_derivative.troubled_zones:
            plt.plot(x[zone], zeros_like(x[zone]) + 0.02, c="k", lw=1.0)
    except AttributeError:
        pass

    plt.ylim(0.0, 1.2)
    plt.ylabel(r"$\rho$")
    plt.xlabel(r"$x$")
    plt.tight_layout(pad=0.05)
    if filename is None:
        plt.show()
    else:
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
    parser.add_argument(
        "--setup",
        type=str,
        default="sod",
        choices=["sod", "density_wave"],
    )
    args = parser.parse_args()

    max_wavespeed = 2.0
    dx, grid = gauss_grid(0.0, 1.0, args.zones, order=args.order, retstep=True)
    dt = dx / max_wavespeed * args.cfl
    x = grid
    p = initial_condition(x, args.setup)
    u = apply_along_axis(prim_to_cons, axis=1, arr=p)
    w = field_to_weights(u)
    w = add_guard_zones(w, args.setup)
    t = 0.0
    n = 0
    L = lambda w: time_derivative(
        w, grid, dx=dx, setup=args.setup, ck=args.limiter_param
    )

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
            plot_soln(x, w[1:-1], filename=f"frame.{n:04d}.png")

    if args.plot == "end":
        plot_soln(x, w[1:-1])


if __name__ == "__main__":
    main()
