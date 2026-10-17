# Local exploration: rotation invariance

AI-generated working notes on the upstream manuscript *Rotation invariance from lattice symmetry and reflection positivity*.

The underlying project is Ethan Sussman's collaboration with Levent Alpoge, building on the literature cited by that project. These supplementary notes retain that attribution and do not represent the upstream authors' endorsement.

Base revision: `785a68fdeef78a1da854c4807fe1f17b597d05e7`.

## Contents

- [Chinese mathematical note](local-note.zh.md): a conditional frequency-elimination lemma with a relaxed growth threshold, a finite-dimensional example at the threshold, and a shorter operator proof of the bilinear no-contact lemma.
- [Exact finite checks](check_boundary_example.py): Fourier coefficients and Gaussian moment-matrix pivots, using only the Python standard library.
- [Planar localization note (Chinese)](planar-localization.zh.md), 2026-10-18: a direct operator-norm construction of point insertions using commuting disk averages, with a quadratic regularization error and an optional fourth-order extrapolation.
- [Exact disk-average checks](check_planar_averaging.py): polar-integral moments, the Laplacian expansion, commuting polynomial averages, and cancellation of the quadratic bias.

Run the checks from the repository root:

```sh
python3 explorations/rotation-invariance/check_boundary_example.py
python3 explorations/rotation-invariance/check_planar_averaging.py
```

## Status

The note contains proposed complete local arguments, but has not received independent mathematical review. The finite checks are consistency checks, not a proof checker. The boundary example concerns the abstract frequency/positivity mechanism, not the full field-theoretic theorem.

No claim of novelty, a complete audit of the rotation theorem, completion of the YM4 construction, or a mass-gap proof is made.

## Next questions

1. Review the new norm-convergence argument independently. It starts from the bounded slab insertion and covariance identities and supplies a replacement for the compactness step in Appendix A.5; it does not complete an audit of A.3–A.7.
2. Audit the global angle continuation and kernel identification in Sections 4.2–4.4, especially how locally continued kernels are identified with the original moments.
3. Track precisely which field and normalization in a concrete model satisfy Gaussian domination, reflection positivity, and ordinary closedness.
4. Record proposed changes separately from the upstream text, with their exact hypotheses and verification status.
