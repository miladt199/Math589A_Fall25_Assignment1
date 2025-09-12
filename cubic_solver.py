
from __future__ import annotations
import math
import cmath
from typing import List, Iterable

_TOL = 1e-12

def _is_close(a: float, b: float = 0.0, tol: float = _TOL) -> bool:
    return abs(a - b) <= tol

def csqrt(z) -> complex:
    """Principal square root without sqrt() or **0.5."""
    zc = complex(z)
    if zc == 0:
        return 0j
    return cmath.exp(0.5 * cmath.log(zc))

def rsqrt_pos(x: float) -> float:
    """Real sqrt for strictly positive x without sqrt()."""
    return math.exp(0.5 * math.log(x))

def _real_cuberoot(x: float) -> float:
    return math.copysign(abs(x) ** (1.0/3.0), x)

def _cleanup(z: complex, tol: float = 1e-10) -> complex:
    zr, zi = z.real, z.imag
    if abs(zi) < tol:
        zi = 0.0
    if abs(zr) < tol:
        zr = 0.0
    return complex(zr, zi)

def _as_complex_list(vals: Iterable) -> List[complex]:
    return [complex(v) for v in vals]

def _sorted_roots(roots: Iterable[complex]) -> List[complex]:
    # sort by (rounded real, rounded imag) for deterministic order
    return sorted(
        (_cleanup(complex(r)) for r in roots),
        key=lambda z: (round(z.real, 12), round(z.imag, 12))
    )

# Precomputed √3 without using sqrt()
SQRT3 = (csqrt(3.0)).real

def solve_cubic(a: float, b: float, c: float, d: float) -> List[complex]:
    """
    Solve a x^3 + b x^2 + c x + d = 0 using trig/hyperbolic substitutions.
    Returns exactly three complex roots (with multiplicities), sorted.
    """
    roots: List[complex] = []

    # Handle degenerate leading coefficient -> quadratic/linear
    if _is_close(a, 0.0):
        if _is_close(b, 0.0):
            if _is_close(c, 0.0):
                roots = []  # constant == d
            else:
                roots = [complex(-d / c)]  # linear
        else:
            # Quadratic: b x^2 + c x + d = 0
            D = c*c - 4.0*b*d
            sqrtD = csqrt(D)
            r1 = (-c + sqrtD) / (2.0*b)
            r2 = (-c - sqrtD) / (2.0*b)
            roots = [r1, r2]

        # Pad to 3 and return sorted complex
        while len(roots) < 3:
            roots.append(roots[-1] if roots else 0j)
        return _sorted_roots(roots)

    # Normalize to monic: x^3 + A x^2 + B x + C
    A = b / a
    B = c / a
    C = d / a

    # Depressed: x = t - A/3 -> t^3 + p t + q = 0
    p = B - (A*A)/3.0
    q = (2.0*A*A*A)/27.0 - (A*B)/3.0 + C

    # Discriminant Δ = (q/2)^2 + (p/3)^3
    Delta = (q*q)/4.0 + (p*p*p)/27.0

    if _is_close(p, 0.0) and _is_close(q, 0.0):
        roots_t = [0.0, 0.0, 0.0]
    elif Delta < -_TOL:
        # Three real roots (Casus irreducibilis): trig form
        r = 2.0 * rsqrt_pos(-p/3.0)
        arg = (3.0*q/(2.0*p)) * rsqrt_pos(-3.0/p)
        arg = max(-1.0, min(1.0, arg))
        phi = math.acos(arg)
        roots_t = [r * math.cos((phi - 2.0*math.pi*k)/3.0) for k in (0, 1, 2)]
    else:
        # One real root: hyperbolic substitutions (no Cardano radicals)
        if p > 0.0:
            denom = 2.0 * (p * rsqrt_pos(p))  # p^{3/2}
            rhs = -(3.0*SQRT3*q) / denom
            u = (1.0/3.0) * math.asinh(rhs)
            t0 = 2.0 * rsqrt_pos(p/3.0) * math.sinh(u)
        elif p < 0.0:
            pm = -p
            denom = 2.0 * (pm * rsqrt_pos(pm))  # |p|^{3/2}
            rhs = (3.0*SQRT3*abs(q)) / denom
            u = (1.0/3.0) * math.acosh(rhs)
            t0 = math.copysign(1.0, -q) * 2.0 * rsqrt_pos(pm/3.0) * math.cosh(u)
        else:
            t0 = -_real_cuberoot(q)

        # Remaining two from quadratic factor t^2 + t0 t + (t0^2 + p) = 0
        Dq = t0*t0 - 4.0*(t0*t0 + p)
        sqrtDq = csqrt(Dq)
        roots_t = [t0, (-t0 + sqrtDq)/2.0, (-t0 - sqrtDq)/2.0]

    # Undo the shift: x = t - A/3
    shift = A/3.0
    roots = [complex(t) - shift for t in roots_t]

    return _sorted_roots(roots)


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
