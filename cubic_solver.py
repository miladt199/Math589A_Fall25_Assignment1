from __future__ import annotations
import math
import cmath
from typing import List

_TOL = 1e-12

def _is_close(a: float, b: float = 0.0, tol: float = _TOL) -> bool:
    return abs(a - b) <= tol

def _cleanup(z: complex, tol: float = 1e-10) -> complex:
    zr, zi = z.real, z.imag
    if abs(zi) < tol:
        zi = 0.0
    return complex(0.0 if abs(zr) < tol else zr, zi)

def csqrt(z) -> complex:
    """
    Principal square root without using sqrt() or **0.5.
    Uses exp(0.5*log(z)). Works for real/complex. csqrt(0)=0.
    """
    zc = complex(z)
    if zc == 0:
        return 0j
    return cmath.exp(0.5 * cmath.log(zc))

def rsqrt_pos(x: float) -> float:
    """
    Real sqrt for strictly positive x without sqrt().
    Uses exp(0.5*log(x)). Assumes x>0.
    """
    return math.exp(0.5 * math.log(x))

def _real_cuberoot(x: float) -> float:
    # exact real cube root with sign (ok to use 1/3 power)
    return math.copysign(abs(x) ** (1.0/3.0), x)

# Precompute constants without sqrt()
SQRT3 = (csqrt(3.0)).real  # == √3

def solve_cubic(a: float, b: float, c: float, d: float) -> List[complex]:
    """
    Solve a x^3 + b x^2 + c x + d = 0 using trig/hyperbolic substitutions.
    Returns three complex roots (with multiplicities).
    """
    if _is_close(a, 0.0):
        # Degenerates to quadratic / linear
        if _is_close(b, 0.0):
            if _is_close(c, 0.0):
                return []  # constant
            return [_cleanup(-d / c)]
        # Quadratic: b x^2 + c x + d = 0
        D = c*c - 4.0*b*d
        sqrtD = csqrt(D)
        r1 = (-c + sqrtD) / (2.0*b)
        r2 = (-c - sqrtD) / (2.0*b)
        return [_cleanup(r1), _cleanup(r2)]

    # Normalize to monic
    A = b / a
    B = c / a
    C = d / a

    # Depressed cubic: x = t - A/3 -> t^3 + p t + q = 0
    p = B - (A*A)/3.0
    q = (2.0*A*A*A)/27.0 - (A*B)/3.0 + C

    # Discriminant
    Delta = (q*q)/4.0 + (p*p*p)/27.0

    roots_t: List[complex] = []

    if _is_close(p, 0.0) and _is_close(q, 0.0):
        roots_t = [0.0, 0.0, 0.0]
    elif Delta < -_TOL:
        # Three real roots (Casus irreducibilis): trigonometric form
        # r = 2 * sqrt(-p/3)
        r = 2.0 * rsqrt_pos(-p/3.0)
        # arg = (3q/(2p))*sqrt(-3/p)
        arg = (3.0*q/(2.0*p)) * rsqrt_pos(-3.0/p)
        # clamp
        arg = max(-1.0, min(1.0, arg))
        phi = math.acos(arg)
        roots_t = [r * math.cos((phi - 2.0*math.pi*k)/3.0) for k in (0,1,2)]
    else:
        # One real root; use hyperbolic substitution (no Cardano radicals)
        if p > 0.0:
            # t = 2*sqrt(p/3) * sinh(u),  sinh(3u) = -(3√3 q)/(2 p^{3/2})
            denom = 2.0 * (p * rsqrt_pos(p))  # p^{3/2} without sqrt
            rhs = -(3.0*SQRT3*q) / denom
            u = (1.0/3.0) * math.asinh(rhs)
            t0 = 2.0 * rsqrt_pos(p/3.0) * math.sinh(u)
        elif p < 0.0:
            # t = sgn(-q)*2*sqrt(-p/3)*cosh(u),  cosh(3u) = (3√3|q|)/(2 |p|^{3/2})
            pm = -p
            denom = 2.0 * (pm * rsqrt_pos(pm))  # |p|^{3/2}
            rhs = (3.0*SQRT3*abs(q)) / denom
            u = (1.0/3.0) * math.acosh(rhs)
            t0 = math.copysign(1.0, -q) * 2.0 * rsqrt_pos(pm/3.0) * math.cosh(u)
        else:
            # p == 0 -> t^3 + q = 0
            t0 = -_real_cuberoot(q)

        # Quadratic factor for remaining two roots
        A2 = 1.0
        B2 = t0
        C2 = t0*t0 + p
        Dq = B2*B2 - 4.0*A2*C2
        sqrtDq = csqrt(Dq)
        roots_t = [t0, (-B2 + sqrtDq)/2.0, (-B2 - sqrtDq)/2.0]

    # Undo shift
    shift = A/3.0
    roots = [_cleanup(complex(t) - shift) for t in roots_t]
    # Always return 3 items
    if len(roots) == 2:
        roots.append(roots[-1])
    if len(roots) == 1:
        roots += [roots[0], roots[0]]
    return roots


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
