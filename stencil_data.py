import argparse


def generate(data, max_order):
    from numpy.polynomial.legendre import leggauss

    if data == "gauss-nodes":
        print("double gauss_quadrature_node(int order, int index)")
    if data == "gauss-weights":
        print("double gauss_quadrature_weight(int order, int index)")
    print("{")
    print("    switch (order) {")
    for order in range(1, max_order + 1):
        x, w = leggauss(order)
        print(f"    case {order}:")
        print("        switch (index) {")
        for n in range(order):
            if data == "gauss-nodes":
                print(f"        case {n}: return {x[n]:+.14f};")
            if data == "gauss-weights":
                print(f"        case {n}: return {w[n]:+.14f};")
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
