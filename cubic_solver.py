import math, cmath

def solve_cubic(a, b, c, d):
    """Solve a*x^3 + b*x^2 + c*x + d = 0
    Returns list of 1..3 roots (complex if needed).
    """
    roots = []

    if abs(a) < 1e-14:  # Degenerate -> quadratic
        if abs(b) < 1e-14:  # Linear
            if abs(c) > 1e-14:
                roots.append(-d/c)
            return roots
        disc = c*c - 4*b*d
        s = cmath.sqrt(disc)
        roots.append((-c + s)/(2*b))
        roots.append((-c - s)/(2*b))
        return roots

    # Normalize coefficients
    A = b/a
    B = c/a
    C = d/a

    # Depressed cubic: x = t - A/3  ->  t^3 + p t + q = 0
    shift = A/3.0
    p = B - A*A/3.0
    q = 2*A*A*A/27.0 - A*B/3.0 + C

    # Discriminant
    half_q = 0.5*q
    disc = half_q*half_q + (p/3.0)**3

    if abs(disc) < 1e-14:  # multiple roots
        if abs(half_q) < 1e-14:
            t1 = 0.0
            roots = [t1]*3
        else:
            u = (-half_q)**(1/3)
            t1 = 2*u
            t2 = -u
            roots = [t1, t2, t2]
    elif disc > 0:  # one real + two complex
        sqrt_disc = cmath.sqrt(disc)
        u = (-half_q + sqrt_disc)**(1/3)
        v = (-half_q - sqrt_disc)**(1/3)
        t1 = u+v
        # complex cube roots of unity
        omega = complex(-0.5, math.sqrt(3)/2)
        omega2 = complex(-0.5, -math.sqrt(3)/2)
        t2 = u*omega + v*omega2
        t3 = u*omega2 + v*omega
        roots = [t1, t2, t3]
    else:  # three real roots (casus irreducibilis)
        r = math.sqrt(-p/3.0)
        arg = (-half_q)/(r**3)
        arg = max(-1.0, min(1.0, arg))  # clamp
        theta = math.acos(arg)
        t1 = 2*r*math.cos(theta/3.0)
        t2 = 2*r*math.cos((theta+2*math.pi)/3.0)
        t3 = 2*r*math.cos((theta+4*math.pi)/3.0)
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
