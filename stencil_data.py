import argparse


def generate(data, max_order):
    from numpy.polynomial.legendre import leggauss, Legendre
    from numpy import zeros

    if data == "gauss-nodes":
        print("double gauss_quadrature_node(int order, int index)")
    if data == "gauss-weights":
        print("double gauss_quadrature_weight(int order, int index)")
    if data == "lobatto-points":
        print("double lobatto_point(int order, int index)")

    print("{")
    print("    switch (order) {")
    for order in range(1, max_order + 1):
        x, w = leggauss(order)

        # A reference for Lobatto quadrature nodes and weights:
        # https://keisan.casio.com/exec/system/1280801905#
        lobatto = zeros(order + 1)
        lobatto[1:-1] = Legendre([0] * order + [1]).deriv().roots()
        lobatto[+0] = -1.0
        lobatto[-1] = +1.0
        lobatto[abs(lobatto) < 1e-15] = 0.0

        print(f"    case {order}:")
        print("        switch (index) {")

        for n in range(order + (data == "lobatto-points")):
            if data == "gauss-nodes":
                print(f"        case {n}: return {x[n]:+.14f};")
            if data == "gauss-weights":
                print(f"        case {n}: return {w[n]:+.14f};")
            if data == "lobatto-points":
                print(f"        case {n}: return {lobatto[n]:+.14f}")

        print("        }")
        print("        break;")
    print("    }")
    print("    return 0.0;")
    print("}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-order", default=3, type=int)
    args = parser.parse_args()

    generate("gauss-nodes", args.max_order)
    print()
    generate("gauss-weights", args.max_order)
    print()
    generate("lobatto-points", args.max_order)
