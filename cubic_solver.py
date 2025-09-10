import math, cmath

SQRT3_OVER_2 = 0.8660254037844386467637231707529361834714026269051903140279034897  # sqrt(3)/2

def _half_power(z):
    """Principal square-root via exp(log)/2 (no sqrt or **0.5)."""
    if z == 0:
        return 0j
    return cmath.exp(0.5 * cmath.log(z))

def _third_power(z):
    """Principal cube-root via exp(log)/3 (no **(1/3))."""
    if z == 0:
        return 0j
    return cmath.exp(cmath.log(z) / 3.0)

def solve_cubic(a, b, c, d):
    """Solve a*x^3 + b*x^2 + c*x + d = 0
    Returns list of 1..3 roots (complex if needed).
    (No explicit sqrt/**0.5 or cube-roots like **(1/3)/pow(x,1/3).)
    """
    roots = []

    # Degenerate -> quadratic or linear
    if abs(a) < 1e-14:
        if abs(b) < 1e-14:
            if abs(c) > 1e-14:
                roots.append(-d / c)
            return roots
        disc = c*c - 4*b*d
        s = _half_power(disc)  # no sqrt()
        roots.append((-c + s) / (2*b))
        roots.append((-c - s) / (2*b))
        return roots

    # Normalize coefficients
    A = b/a
    B = c/a
    C = d/a

    # Depressed cubic: x = t - A/3  ->  t^3 + p t + q = 0
    shift = A / 3.0
    p = B - A*A/3.0
    q = 2*A*A*A/27.0 - A*B/3.0 + C

    # Discriminant
    half_q = 0.5 * q
    disc = half_q*half_q + (p/3.0)**3

    if abs(disc) < 1e-14:  # multiple roots
        if abs(half_q) < 1e-14:
            t1 = 0.0
            roots = [t1, t1, t1]
        else:
            u = _third_power(-half_q)  # no **(1/3)
            t1 = 2*u
            t2 = -u
            roots = [t1, t2, t2]

    elif disc > 0:  # one real + two complex
        sdisc = _half_power(disc)  # no sqrt()
        u = _third_power(-half_q + sdisc)
        v = _third_power(-half_q - sdisc)
        t1 = u + v
        omega = complex(-0.5,  SQRT3_OVER_2)
        omega2 = complex(-0.5, -SQRT3_OVER_2)
        t2 = u*omega  + v*omega2
        t3 = u*omega2 + v*omega
        roots = [t1, t2, t3]

    else:  # three real roots (casus irreducibilis): trigonometric form
        # r = sqrt(-p/3) but computed without sqrt()
        # r = sqrt(-p/3) > 0 ideally
        rp = -p / 3.0
        if rp <= 0:
            rp = max(rp, 1e-30)  # tiny nudge for round-off
        r = math.exp(0.5 * math.log(rp))
        r3 = r * r * r

        if r3 <= 1e-300:  # ultra-defensive; should not happen for Δ<0
            arg = 1.0 if q.real >= 0 else -1.0
        else:
            arg = q.real / (2.0 * r3)

        # Clamp to domain
        if arg > 1.0: arg = 1.0
        if arg < -1.0: arg = -1.0
        theta = math.acos(arg)
        t1 = 2*r*math.cos(theta/3.0)
        t2 = 2*r*math.cos((theta + 2*math.pi)/3.0)
        t3 = 2*r*math.cos((theta + 4*math.pi)/3.0)
        roots = [t1, t2, t3]

    # Shift back: x = t - A/3
    return [t - shift for t in roots]


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
