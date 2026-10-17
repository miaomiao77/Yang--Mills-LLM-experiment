"""Finite consistency checks for the abstract frequency zero-window lemma.

Standard library only. This is *not* a proof checker for the Yang--Mills paper.
"""

from itertools import combinations_with_replacement


def check_bounded_convex_sequences(max_k=8):
    """Exhaust nonnegative nondecreasing integer increments with fixed sum."""
    cases = 0
    for k in range(1, max_k + 1):
        for b in range(k + 1):
            j0 = k - b
            seen = 0
            for increments in combinations_with_replacement(range(b + 1), k):
                if sum(increments) != b:
                    continue
                seen += 1
                a = [0]
                for delta in increments:
                    a.append(a[-1] + delta)
                assert a[k] == b
                assert all(a[j] == 0 for j in range(j0 + 1)), (k, b, a)
            assert seen >= 1
            cases += seen
    return cases


def hinge(n, q, cutoff):
    return q * max(0, n - cutoff)


def check_hinge_examples(max_k=8, q=4):
    for k in range(1, max_k + 1):
        for b in range(k + 1):
            j0 = k - b
            cutoff = 2 * j0
            assert hinge(2 * k, q, cutoff) == 2 * q * b
            assert hinge(2 * j0 + 1, q, cutoff) == q
            for j in range(2 * k + 8):
                assert hinge(2 * j, q, cutoff) % (2 * q) == 0
            for a in range(2 * k + 8):
                for c in range(2 * k + 8):
                    assert (
                        2 * hinge(a + c, q, cutoff)
                        <= hinge(2 * a, q, cutoff) + hinge(2 * c, q, cutoff)
                    )

    assert hinge(20, 4, 10) == 40
    assert all(hinge(n, 4, 10) == 0 for n in range(11))
    assert hinge(11, 4, 10) == 4
    assert hinge(20, 4, 8) == 48
    assert hinge(10, 4, 8) == 8


if __name__ == "__main__":
    cases = check_bounded_convex_sequences()
    check_hinge_examples()
    print(f"Exact discrete checks passed ({cases} increment patterns).")
    print("Sharp hinge examples passed all tested midpoint inequalities.")
