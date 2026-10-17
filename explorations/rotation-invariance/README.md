# Local exploration: rotation invariance

AI-generated working notes on the upstream manuscript *Rotation invariance from lattice symmetry and reflection positivity*.

The underlying project is Ethan Sussman's collaboration with Levent Alpoge, building on the literature cited by that project. These supplementary notes retain that attribution and do not represent the upstream authors' endorsement.

Base revision: `785a68fdeef78a1da854c4807fe1f17b597d05e7`.

## Contents

- [Chinese mathematical note](local-note.zh.md): a conditional frequency-elimination lemma with a relaxed growth threshold, a finite-dimensional example at the threshold, and a shorter operator proof of the bilinear no-contact lemma.
- [Exact finite checks](check_boundary_example.py): Fourier coefficients and Gaussian moment-matrix pivots, using only the Python standard library.

Run the checks from the repository root:

```sh
python3 explorations/rotation-invariance/check_boundary_example.py
```

## Status

The note contains proposed complete local arguments, but has not received independent mathematical review. The finite checks are consistency checks, not a proof checker. The boundary example concerns the abstract frequency/positivity mechanism, not the full field-theoretic theorem.

No claim of novelty, a complete audit of the rotation theorem, completion of the YM4 construction, or a mass-gap proof is made.

## Next questions

1. Independently audit the insertion operators and planar localization in Appendix A.3–A.7 of the upstream short manuscript.
2. Track precisely which field and normalization in a concrete model satisfy Gaussian domination, reflection positivity, and ordinary closedness.
3. Record any proposed repair separately from the upstream text, with its exact hypotheses and verification status.
