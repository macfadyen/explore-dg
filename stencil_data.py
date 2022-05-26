import argparse


def generate(data, max_order):
    from numpy.polynomial.legendre import leggauss

    print("switch (order) {")
    for order in range(1, max_order + 1):
        x, w = leggauss(order)
        print(f"case {order}:")
        print("    switch (index) {")
        for n in range(order):
            print(f"    case {n}:")
            if data == "gauss-nodes":
                print(f"        return {x[n]:+.14f};")
            if data == "gauss-weights":
                print(f"        return {w[n]:+.14f};")
        print("    }")
    print("}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--data", choices=["gauss-nodes", "gauss-weights"], default="gauss-nodes"
    )
    parser.add_argument("--max-order", default=3, type=int)
    args = parser.parse_args()

    generate(args.data, args.max_order)
