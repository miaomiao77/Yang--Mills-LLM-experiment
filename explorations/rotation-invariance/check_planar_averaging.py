"""Exact algebra checks for planar-localization.zh.md (standard library only).

These check disk averages and Taylor coefficients, not infinite-dimensional
operator convergence, reflection positivity, or a field-theoretic theorem.
"""

from fractions import Fraction as Q
from functools import lru_cache
from math import comb, factorial


@lru_cache(None)
def angular_even(a, b):
    """Mean cos(theta)^(2a) sin(theta)^(2b), by integration by parts."""
    if a:
        return Q(2 * a - 1, 2 * (a + b)) * angular_even(a - 1, b)
    if b:
        return Q(2 * b - 1, 2 * b) * angular_even(0, b - 1)
    return Q(1)


def disk_moment(a, b, radius=Q(1)):
    """Normalized disk integral: angular recurrence times radial integral."""
    if a % 2 or b % 2:
        return Q(0)
    degree = a + b
    return radius**degree * angular_even(a // 2, b // 2) * Q(2, degree + 2)


def add_term(poly, key, value):
    poly[key] = poly.get(key, Q(0)) + value
    if not poly[key]:
        del poly[key]


def disk_average(poly, radius):
    """Integrate p(z + u) over u in a disk, using the binomial theorem."""
    out = {}
    for (a, b), coefficient in poly.items():
        for i in range(a + 1):
            for j in range(b + 1):
                add_term(out, (a - i, b - j), coefficient * comb(a, i)
                         * comb(b, j) * disk_moment(i, j, radius))
    return out


def laplacian(poly):
    out = {}
    for (a, b), coefficient in poly.items():
        if a >= 2:
            add_term(out, (a - 2, b), coefficient * a * (a - 1))
        if b >= 2:
            add_term(out, (a, b - 2), coefficient * b * (b - 1))
    return out


def differential_average(poly, radius):
    """Apply the claimed Laplacian series to a polynomial."""
    out, derivative, m = {}, poly, 0
    while derivative:
        factor = radius**(2 * m) / (4**m * factorial(m) * factorial(m + 1))
        for key, coefficient in derivative.items():
            add_term(out, key, factor * coefficient)
        derivative = laplacian(derivative)
        m += 1
    return out


def main():
    assert disk_moment(0, 0) == 1
    assert disk_moment(2, 0) == disk_moment(0, 2) == Q(1, 4)
    assert disk_moment(4, 0) == Q(1, 8)
    assert disk_moment(2, 2) == Q(1, 24)
    assert disk_moment(1, 1) == 0
    assert disk_moment(2, 0) / factorial(2) == Q(1, 8)
    assert disk_moment(4, 0) / factorial(4) == Q(1, 192)
    print("Disk moments and Taylor coefficients 1/8, 1/192: passed")

    radii = (Q(1, 7), Q(2, 9))
    count = 0
    for degree in range(11):
        for a in range(degree + 1):
            poly = {(a, degree - a): Q(1)}
            for radius in radii:
                assert disk_average(poly, radius) == differential_average(poly, radius)
            left = disk_average(disk_average(poly, radii[0]), radii[1])
            right = disk_average(disk_average(poly, radii[1]), radii[0])
            assert left == right
            count += 1
    print(f"Laplacian formula and commuting averages: {count} monomials passed")

    # A finite-dimensional local analytic illustration, not a field law.
    # h(z) = exp(-z_t H + i z_x P) A exp(z_t H - i z_x P).
    energies, momenta = (2, 3, 6), (1, -2, 0)
    assert all(h >= abs(p) for h, p in zip(energies, momenta))
    for i in range(3):
        for j in range(3):
            a, k = energies[j] - energies[i], momenta[i] - momenta[j]
            lam = a * a - k * k
            for m in range(6):
                integrated = sum(
                    Q(a**(2 * u) * (-1)**(m - u) * k**(2 * (m - u)),
                      factorial(2 * u) * factorial(2 * (m - u)))
                    * disk_moment(2 * u, 2 * (m - u))
                    for u in range(m + 1)
                )
                claimed = Q(lam**m, 4**m * factorial(m) * factorial(m + 1))
                assert integrated == claimed
                extrapolated = (4 * Q(1, 2)**(2 * m) - 1) * integrated / 3
                if m == 1:
                    assert extrapolated == 0
                if m == 2:
                    assert extrapolated == -Q(lam**2, 768)
    print("Finite matrix illustration and quadratic-bias cancellation: passed")
    assert Q(256, 3) * (4 * Q(1, 2)**4 + 1) / 3 == Q(320, 9)
    print("Fourth-order remainder constant 320/9: passed")
    print("Exact algebra checks passed; the analytic proof is in the note.")


if __name__ == "__main__":
    main()
