import math, cmath

def _cbrt(z):
    # Real cube root for real inputs, complex principal cube root otherwise
    if isinstance(z, complex) or (isinstance(z, float) and z != z):  # NaN check pass-through
        return z**(1/3)
    if z >= 0:
        return z**(1/3)
    return -((-z)**(1/3))

def solve_cubic(a, b, c, d, tol=1e-14):
    """Solve a*x^3 + b*x^2 + c*x + d = 0
    Returns list of 1..3 roots (complex if needed).
    """
    # Handle degenerate cases
    if abs(a) < tol:
        # quadratic: b*x^2 + c*x + d = 0
        if abs(b) < tol:
            # linear: c*x + d = 0
            if abs(c) < tol:
                return []  # constant equation: no root (or infinite if d==0)
            return [ -d / c ]
        disc = c*c - 4*b*d
        if disc >= 0:
            r1 = (-c + math.sqrt(disc)) / (2*b)
            r2 = (-c - math.sqrt(disc)) / (2*b)
            return [r1, r2]
        else:
            r = math.sqrt(-disc) / (2*b)
            return [ complex(-c/(2*b), r), complex(-c/(2*b), -r) ]

    # ################################################################
    # Your solution goes here
    # ################################################################

    # Cleanup tiny imaginary/noise
    clean = []
    for z in roots:
        if isinstance(z, complex):
            re, im = z.real, z.imag
            if abs(im) <= 1e-12:
                clean.append(re)
            else:
                clean.append(z)
        else:
            clean.append(z)
    return clean


def main():
    tests = [
        (1, 0, 0, -1),     # roots of x^3 - 1 = 0 (1 and two complex cube roots)
        (1, -6, 11, -6),   # roots [1.0, 2.0, 3.0]
    ]
    for a, b, c, d in tests:
        roots = solve_cubic(a, b, c, d)
        print(f"solve_cubic({a}, {b}, {c}, {d}) -> {roots}")

if __name__ == "__main__":
    main()
