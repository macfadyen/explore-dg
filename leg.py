"""
Illustrates use of legendre polynomials to model a function

This program provides an option to disable the use of the numpy.einsum to
transform the modeled function between spatial and spectral representations.

Author: Jonathan Zrake
"""

from numpy.polynomial.legendre import Legendre, leggauss
from numpy import (
    allclose,
    array,
    cos,
    einsum,
    exp,
    linspace,
    sin,
    zeros,
    zeros_like,
)


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


def field_to_weights(field, use_einsum):
    """
    Convert from an array of field values.

    The conversion is done at the sub-cell Gauss quadrature points, to an
    array of weights.
    """
    order = field.shape[-1]
    x, w = leggauss(order)
    phi = array([leg(x, n) for n in range(order)])
    weights = zeros_like(field)

    if use_einsum:
        return einsum("iq,nq,q,...", field, phi, w, 0.5)

    else:
        for i in range(field.shape[0]):
            for n in range(order):
                for q in range(order):
                    weights[i, n] += field[i, q] * phi[n, q] * w[q] * 0.5

        return weights


def weights_to_field(weights, use_einsum):
    """
    Convert from an array of weights to an array of field values.

    The conversion yields function values at the sub-cell Gauss quadrature
    points.
    """
    order = weights.shape[-1]
    x, w = leggauss(order)
    phi = array([leg(x, n) for n in range(order)])

    if use_einsum:
        return einsum("in,nq", weights, phi)

    else:
        field = zeros_like(weights)

        for i in range(field.shape[0]):
            for n in range(order):
                for q in range(order):
                    field[i, q] += weights[i, n] * phi[n, q]

        return field


def function_to_model(x):
    return cos(5 * x) * sin(17 * x) * exp(-((x - 0.5) ** 2) / 0.1)


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
        "--no-einsum",
        action="store_true",
        help="do projections manually, not with numpy.einsum",
    )
    args = parser.parse_args()

    dx, grid = gauss_grid(0.0, 1.0, args.zones, order=args.order, retstep=True)
    x_fine = linspace(0.0, 1.0, 1000)
    y_fine = function_to_model(x_fine)

    x = grid
    y = function_to_model(x)
    w = field_to_weights(y, use_einsum=not args.no_einsum)
    y1 = weights_to_field(w, use_einsum=not args.no_einsum)

    assert allclose(y1, y, rtol=1e-13, atol=1e-13)

    plt.plot(x_fine, y_fine, "-", label=r"$f(x)$")
    plt.plot(x.flat, y.flat, "-o", mfc="none", label=r"$f(x)$ at Gauss points")
    plt.legend()
    plt.show()


if __name__ == "__main__":
    main()
