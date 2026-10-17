"""Exact finite checks for the spin-four Gaussian boundary example.

These checks validate a few formulas in local-note.zh.md. The all-order
positivity argument is the Gram representation in the note, not this script.
"""

from fractions import Fraction
from math import comb, factorial, prod


def gaussian_moment(order):
    """Return E[X**order] for X ~ N(0, 1), exactly."""
    return 0 if order % 2 else prod(range(1, order, 2))


def fourier_coefficients(order):
    """Coefficients of E[X**order] * cos(4*theta)**order."""
    moment = gaussian_moment(order)
    if not moment:
        return {}
    return {
        4 * (order - 2 * k): Fraction(moment * comb(order, k), 2**order)
        for k in range(order + 1)
    }


def moment_matrix_pivots(size):
    """Exact symmetric elimination of the Gaussian moment matrix."""
    matrix = [
        [Fraction(gaussian_moment(i + j)) for j in range(size)]
        for i in range(size)
    ]
    pivots = []
    for k in range(size):
        pivot = matrix[k][k]
        assert pivot > 0
        pivots.append(pivot)
        for i in range(k + 1, size):
            for j in range(k + 1, size):
                matrix[i][j] -= matrix[i][k] * matrix[k][j] / pivot
    return pivots


def main():
    for order in (2, 4, 6):
        coefficients = fourier_coefficients(order)
        assert max(map(abs, coefficients)) == 4 * order
        assert all(frequency % 8 == 0 for frequency in coefficients)
        assert sum(coefficients.values()) == gaussian_moment(order)
        print(f"N={order}: frequencies={sorted(coefficients)}")
    assert all(not fourier_coefficients(n) for n in (1, 3, 5))
    pivots = moment_matrix_pivots(7)
    assert pivots == [Fraction(factorial(k)) for k in range(7)]
    print("Gaussian moment-matrix LDL pivots:", [str(x) for x in pivots])
    print("Exact finite checks passed; see the note for the all-order arguments.")


if __name__ == "__main__":
    main()
