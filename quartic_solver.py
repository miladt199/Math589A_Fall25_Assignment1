import numpy as np

def solve_quartic(a, b, c, d, e):
    """Solve a*x^4 + b*x^3 + c*x^2 + d*x + e = 0.
    Returns a list of 1..4 roots (real numbers or complex numbers).
    If the leading coefficients are (near) zero the function will
    handle lower-degree polynomials automatically.
    tol: tolerance for treating small coefficients/imaginary parts as zero.
    """
    coeffs = [a, b, c, d, e]
    roots = []

    # ################################################################
    # Your solution goes here
    # ################################################################

    # Sample return statement
    return roots


def main():
    tests = [
        (1, 0, 0, 0, -1),  # roots of x^4 - 1 = 0 (2 real and 2 complex roots)
        (1, 0, 1, 0, -1),  # roots of x^4 + x^2 - 1 = 0 (2 real and 2 complex roots)
        (0, 1, -3, 2, 0),  # a=0 => cubic: x^3 - 3x^2 + 2x = 0  (roots 0,1,2)
    ]
    for a, b, c, d, e in tests:
        roots = solve_quartic(a, b, c, d, e)
        print(f"solve_quartic({a}, {b}, {c}, {d}, {e}) -> {roots}")

if __name__ == "__main__":
    main()
