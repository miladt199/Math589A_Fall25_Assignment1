import math, cmath

SQRT3_OVER_2 = 0.8660254037844386467637231707529361834714026269051903140279034897  # sqrt(3)/2

def _is_close(a, b, tol=1e-12):
    return abs(a - b) <= tol

def _half_power(z):
    # principal "square root": exp( (1/2) log z )
    if z == 0:
        return 0j
    return cmath.exp(0.5 * cmath.log(z))

def _third_power(z):
    # principal "cube root": exp( (1/3) log z ), with zero guard
    if z == 0:
        return 0j
    return cmath.exp(cmath.log(z) / 3.0)

def _solve_linear(b, c):
    return [] if _is_close(b, 0.0) else [-c / b]

def _solve_quadratic(a, b, c):
    if _is_close(a, 0.0):
        return _solve_linear(b, c)
    disc = b*b - 4*a*c
    s = _half_power(disc)             # no sqrt()
    return [(-b + s)/(2*a), (-b - s)/(2*a)]

def solve_cubic(a, b, c, d):
    """Solve a*x^3 + b*x^2 + c*x + d = 0 (no explicit radicals)."""
    if _is_close(a, 0.0):
        return _solve_quadratic(b, c, d)

    # Normalize, depress: x = t - A/3  ->  t^3 + p t + q = 0
    A = b/a
    B = c/a
    C = d/a

    shift = A/3.0
    p = B - A*A/3.0
    q = 2*A*A*A/27.0 - A*B/3.0 + C

    half_q = 0.5*q
    disc = half_q*half_q + (p/3.0)**3

    if abs(disc) < 1e-14:
        if abs(half_q) < 1e-14:
            ts = [0.0, 0.0, 0.0]
        else:
            u = _third_power(-half_q)
            ts = [2*u, -u, -u]
    elif disc > 0:
        # One real, two complex: Cardano via exp(log)/3 (no **(1/3))
        sdisc = _half_power(disc)
        u = _third_power(-half_q + sdisc)
        v = _third_power(-half_q - sdisc)
        t1 = u + v
        w  = complex(-0.5,  SQRT3_OVER_2)
        w2 = complex(-0.5, -SQRT3_OVER_2)
        t2 = u*w  + v*w2
        t3 = u*w2 + v*w
        ts = [t1, t2, t3]
    else:
        # Three real roots (casus irreducibilis): trigonometric form
        # r = sqrt(-p/3)  -> compute without sqrt()
        rp = -p / 3.0
        if rp <= 0:
            rp = max(rp, 1e-30)  # guard tiny negatives
        r = math.exp(0.5 * math.log(rp))  # sqrt(-p/3) without sqrt()
        r3 = r * r * r

        # Correct identity: cos(3θ) = q / (2 r^3)
        arg = (q.real) / (2.0 * r3) if r3 != 0.0 else (1.0 if q.real >= 0 else -1.0)
        if arg > 1.0: arg = 1.0  # clamp to domain
        if arg < -1.0: arg = -1.0
        theta = math.acos(arg)

        t1 = 2 * r * math.cos(theta / 3.0)
        t2 = 2 * r * math.cos((theta + 2 * math.pi) / 3.0)
        t3 = 2 * r * math.cos((theta + 4 * math.pi) / 3.0)
        ts = [t1, t2, t3]

    return [t - shift for t in ts]


def solve_quartic(a, b, c, d, e):
    """Solve a*x^4 + b*x^3 + c*x^2 + d*x + e = 0.
    Returns a list of 1..4 roots (real or complex).
    No explicit sqrt or cube-root usage.
    """
    roots = []

    # Degenerate degree?
    if _is_close(a, 0.0):
        return solve_cubic(b, c, d, e)

    # Normalize to monic: x^4 + B x^3 + C x^2 + D x + E = 0
    B = b / a
    C = c / a
    D = d / a
    E = e / a

    # Depressed quartic: x = y - B/4  ->  y^4 + p y^2 + q y + r = 0
    alpha = B / 4.0
    B2 = B*B
    B3 = B2*B
    B4 = B2*B2

    p = C - 3.0*B2/8.0
    q = D - B*C/2.0 + B3/8.0
    r = E - B*D/4.0 + B2*C/16.0 - 3.0*B4/256.0

    # Biquadratic case
    if _is_close(q, 0.0):
        z_roots = _solve_quadratic(1.0, p, r)  # z = y^2
        for z in z_roots:
            s = _half_power(z)                 # no sqrt()
            roots += [s - alpha, -s - alpha]
        return roots

    # Ferrari resolvent: 8 m^3 - 4 p m^2 - 8 r m + (4 r p - q^2) = 0
    # Monic form: m^3 + (-p/2) m^2 + (-r) m + (r p - q^2/4)/2 = 0
    cb = -p/2.0
    cc = -r
    cd = (r*p - q*q/4.0) / 2.0
    m_roots = solve_cubic(1.0, cb, cc, cd)

    # Choose m (prefer real) with 2m - p >= 0
    m = None
    for mr in m_roots:
        if abs(mr.imag) <= 1e-10 and (2.0*mr.real - p) >= -1e-10:
            m = mr.real
            break
    if m is None:
        m = m_roots[0].real

    R = _half_power(2.0*m - p)  # no sqrt()

    # D1^2 = -(2m + p) + 2q/R,  D2^2 = -(2m + p) - 2q/R
    if abs(R) <= 1e-15:
        D1_sq = -(2.0*m + p)
        D2_sq = D1_sq
    else:
        D1_sq = -(2.0*m + p) + (2.0*q)/R
        D2_sq = -(2.0*m + p) - (2.0*q)/R

    D1 = _half_power(D1_sq)
    D2 = _half_power(D2_sq)

    y1 = 0.5*(+R + D1)
    y2 = 0.5*(+R - D1)
    y3 = 0.5*(-R + D2)
    y4 = 0.5*(-R - D2)

    roots = [y1 - alpha, y2 - alpha, y3 - alpha, y4 - alpha]
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
