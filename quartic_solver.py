

from __future__ import annotations
import cmath
from typing import List
from cubic_solver import solve_cubic, csqrt, _cleanup, _is_close

def _solve_biquadratic(P: float, R: float) -> List[complex]:
    # u^4 + P u^2 + R = 0 -> y^2 + P y + R = 0 with y = u^2
    Dy = P*P - 4.0*R
    sy = csqrt(Dy)
    y1 = (-P + sy)/2.0
    y2 = (-P - sy)/2.0
    roots = []
    for y in (y1, y2):
        syi = csqrt(y)
        roots.extend([syi, -syi])
    return roots

def solve_quartic(a: float, b: float, c: float, d: float, e: float) -> List[complex]:
    """
    Solve a x^4 + b x^3 + c x^2 + d x + e = 0 using Ferrari.
    Returns four complex roots (with multiplicities).
    """
    # Degenerate cases
    if _is_close(a, 0.0):
        return [_cleanup(z) for z in solve_cubic(b, c, d, e)]

    # Normalize (monic)
    A = b / a
    B = c / a
    C = d / a
    D = e / a

    # Depressed quartic: x = u - A/4 -> u^4 + P u^2 + Q u + R = 0
    P = B - 3.0*(A*A)/8.0
    Q = C - 0.5*A*B + (A*A*A)/8.0
    R = D - 0.25*A*C + (A*A*B)/16.0 - 3.0*(A**4)/256.0

    if abs(Q) <= 1e-14:
        u_roots = _solve_biquadratic(P, R)
    else:
        # Resolvent cubic: 2y^3 - P y^2 - 2 R y + (P R - Q^2/4) = 0
        y_roots = solve_cubic(2.0, -P, -2.0*R, (P*R - (Q*Q)/4.0))

        # Choose a good y (prefer real with 2y-P > 0 to keep S nonzero)
        chosen_y = None
        for y in y_roots:
            if abs(y.imag) < 1e-10 and (2.0*y.real - P) > 1e-14:
                chosen_y = y.real
                break
        if chosen_y is None:
            y_roots_sorted = sorted(y_roots, key=lambda z: z.real, reverse=True)
            chosen_y = y_roots_sorted[0].real

        S = csqrt(2.0*chosen_y - P)
        term1 = -2.0*chosen_y - P

        # Beware S ~ 0; chosen_y selection tries to avoid that.
        term_plus  = csqrt(term1 + 2.0*Q / S)
        term_minus = csqrt(term1 - 2.0*Q / S)

        u1 = 0.5 * (-S + term_plus)
        u2 = 0.5 * (-S - term_plus)
        u3 = 0.5 * ( S + term_minus)
        u4 = 0.5 * ( S - term_minus)

        u_roots = [u1, u2, u3, u4]

    # Undo the shift
    shift = A / 4.0
    x_roots = [_cleanup(u - shift) for u in u_roots]
    # Ensure four outputs
    while len(x_roots) < 4:
        x_roots.append(x_roots[-1])
    return x_roots


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
