import math, cmath

def _is_close(a, b, tol=1e-12):
    return abs(a - b) <= tol

def _solve_linear(b, c):
    # b*x + c = 0
    if _is_close(b, 0.0):
        return []
    return [-c / b]

def _solve_quadratic(a, b, c):
    if _is_close(a, 0.0):
        return _solve_linear(b, c)
    disc = b*b - 4*a*c
    s = cmath.sqrt(disc)
    return [(-b + s)/(2*a), (-b - s)/(2*a)]

def solve_cubic(a, b, c, d):
    """Solve a*x^3 + b*x^2 + c*x + d = 0
    Returns list of 1..3 roots (complex if needed).
    """
    roots = []

    # Degenerate to quadratic/linear
    if _is_close(a, 0.0):
        return _solve_quadratic(b, c, d)

    # Normalize to monic: x = t - A/3  =>  t^3 + p t + q = 0
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
            t1 = 0.0
            roots = [t1, t1, t1]
        else:
            u = (-half_q)**(1/3)
            t1 = 2*u
            t2 = -u
            roots = [t1, t2, t2]
    elif disc > 0:
        sqrt_disc = cmath.sqrt(disc)
        u = (-half_q + sqrt_disc)**(1/3)
        v = (-half_q - sqrt_disc)**(1/3)
        t1 = u + v
        omega = complex(-0.5, math.sqrt(3)/2)
        omega2 = complex(-0.5, -math.sqrt(3)/2)
        t2 = u*omega + v*omega2
        t3 = u*omega2 + v*omega
        roots = [t1, t2, t3]
    else:
        r = math.sqrt(-p/3.0)
        arg = (-half_q)/(r**3)
        arg = max(-1.0, min(1.0, arg))
        theta = math.acos(arg)
        t1 = 2*r*math.cos(theta/3.0)
        t2 = 2*r*math.cos((theta+2*math.pi)/3.0)
        t3 = 2*r*math.cos((theta+4*math.pi)/3.0)
        roots = [t1, t2, t3]

    return [t - shift for t in roots]


def solve_quartic(a, b, c, d, e):
    """Solve a*x^4 + b*x^3 + c*x^2 + d*x + e = 0.
    Returns a list of 1..4 roots (real numbers or complex numbers).
    If the leading coefficients are zero the function will
    handle lower-degree polynomials automatically.
    """
    roots = []

    # Degree < 4?
    if _is_close(a, 0.0):
        # Fallback to cubic
        return solve_cubic(b, c, d, e)

    # Normalize to monic: x^4 + B x^3 + C x^2 + D x + E = 0
    B = b / a
    C = c / a
    D = d / a
    E = e / a

    # Depress: x = y - B/4  ->  y^4 + p y^2 + q y + r = 0
    alpha = B / 4.0
    B2 = B*B
    B3 = B2*B
    B4 = B2*B2

    p = C - 3.0*B2/8.0
    q = D - B*C/2.0 + B3/8.0
    r = E - B*D/4.0 + B2*C/16.0 - 3.0*B4/256.0

    # If nearly biquadratic (q ~ 0), solve y^4 + p y^2 + r = 0
    if _is_close(q, 0.0):
        z_roots = _solve_quadratic(1.0, p, r)  # z = y^2
        for z in z_roots:
            s = cmath.sqrt(z)
            roots += [s - alpha, -s - alpha]
        return roots

    # Ferrari via resolvent cubic in m:
    # 8 m^3 - 4 p m^2 - 8 r m + (4 r p - q^2) = 0
    # Make it monic: m^3 + ( -p/2 ) m^2 + ( -r ) m + ( (r p - q^2/4) / 2 ) = 0
    cb = -p/2.0
    cc = -r
    cd = (r*p - q*q/4.0) / 2.0
    m_roots = solve_cubic(1.0, cb, cc, cd)

    # Choose an m that makes 2m - p >= 0 (prefer real)
    m = None
    for mr in m_roots:
        if abs(mr.imag) <= 1e-10 and (2.0*mr.real - p) >= -1e-10:
            m = mr.real
            break
    if m is None:
        # fallback: take real part of the first
        m = m_roots[0].real

    R = cmath.sqrt(2.0*m - p)

    # Build the two auxiliary terms:
    # D1^2 = -(2m + p) + 2q/R
    # D2^2 = -(2m + p) - 2q/R
    # guard for R ~ 0
    if abs(R) <= 1e-15:
        D1_sq = -(2.0*m + p)
        D2_sq = D1_sq
    else:
        D1_sq = -(2.0*m + p) + (2.0*q)/R
        D2_sq = -(2.0*m + p) - (2.0*q)/R

    D1 = cmath.sqrt(D1_sq)
    D2 = cmath.sqrt(D2_sq)

    # Four y-roots
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
