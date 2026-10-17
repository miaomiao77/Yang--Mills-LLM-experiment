# A finite-order zero window for quantized midpoint-convex growth

*Working note, 2026-10-17. Independent companion to the local rotation-invariance exploration; **not** a result about the full Yang–Mills theory.*

## Scope and origin

The [rotation-invariance exploration on the separate branch](https://github.com/miaomiao77/Yang--Mills-LLM-experiment/blob/astris/rotation-invariance-notes/explorations/rotation-invariance/local-note.zh.md) studies an abstract frequency-elimination argument from *Rotation invariance from lattice symmetry and reflection positivity* in the Sussman–Alpoge project. Its asymptotic lemma asks when all Fourier-frequency bounds vanish.

Here we ask a different, finite-order question: **how many low orders must already vanish when one higher even order has a bounded frequency?** This note is a deduction under the *abstract sequence assumptions*; it neither audits their derivation from field theory nor establishes novelty relative to the mathematical literature.

## Assumptions

Fix \(q>0\). Let \((M_n)_{n\geq0}\) be a nonnegative real sequence such that

\[
M_0=0,\qquad 2M_{a+b}\le M_{2a}+M_{2b}
\quad(a,b\in\mathbb Z_{\geq0}), \tag{1}
\]

and with quantized even entries

\[
M_{2j}\in 2q\mathbb Z_{\geq0}\qquad(j\geq0). \tag{2}
\]

In the cited rotation-invariance abstract step, \(q=4\).

## Theorem: finite-order zero window

Fix an integer \(K\geq1\), write

\[
M_{2K}=2qB,\qquad B\in\mathbb Z_{\geq0},
\quad J=\max\{K-B,0\}. \tag{3}
\]

Then

\[
\boxed{M_n=0\quad\text{for every integer }0\le n\le2J.} \tag{4}
\]

When \(0\le B\le K\), this **guaranteed zero window is sharp**: there exists a sequence satisfying (1)–(2) and (3) with \(M_{2J+1}>0\).

Equivalently, for every \(1\le j\le K\),

\[
\boxed{M_{2K}<2q(K-j+1)
\quad\Longrightarrow\quad
M_n=0\quad(0\le n\le 2j).} \tag{5}
\]

The strict threshold in (5) cannot be relaxed to a weak inequality on the abstract assumptions alone.

### Proof

Write \(A_j=M_{2j}\). Taking \(a=j-1,\ b=j+1\) in (1) gives

\[
2A_j\le A_{j-1}+A_{j+1}\quad(j\ge1).
\]

Thus \(\Delta_j=A_j-A_{j-1}\) is nondecreasing. Since
\(\Delta_1=A_1\ge0\), and each \(\Delta_j\) is an integer multiple of \(2q\), every increment is either zero or at least \(2q\).

For \(1\le j\le K\), if \(A_j>0\), at least one increment among
\(\Delta_1,\dots,\Delta_j\) is positive. Monotonicity then forces
\(\Delta_r\ge2q\) for every \(r\ge j\), so

\[
A_K\ge 2q(K-j+1). \tag{6}
\]

Consequently, \(A_K<2q(K-j+1)\) implies \(A_j=0\). Since all increments are nonnegative, this also implies \(A_0=\cdots=A_j=0\).

For any \(n\le2j\), choose \(a=\lfloor n/2\rfloor\) and \(b=\lceil n/2\rceil\), which are at most \(j\). Using (1) once more,

\[
0\le2M_n\le A_a+A_b=0.
\]

This proves (5). If \(B<K\), take \(j=K-B\): then
\(A_K=2qB<2q(B+1)\), giving (4). If \(B\ge K\), then \(J=0\), and (4) is simply the assumption \(M_0=0\). This proves the theorem.

### Sharpness and two different boundaries

For any real cutoff \(c\ge0\), define the hinge sequence

\[
H_n^{(c)}=q(n-c)_+,\qquad(x)_+=\max\{x,0\}. \tag{7}
\]

The map \(x\mapsto q(x-c)_+\) is convex, so it obeys (1) by the midpoint convexity inequality. If \(c\) is an **even nonnegative integer**, (2) also holds.

1. **Sharp zero-window length at fixed \(B\le K\):** choose \(c=2(K-B)=2J\). Then
   \(H_{2K}^{(c)}=2qB\), all entries through order \(2J\) are zero, but
   \(H_{2J+1}^{(c)}=q>0\).
2. **Sharp threshold for a target order \(2j\):** choose \(c=2(j-1)\). Then
   \(H_{2K}^{(c)}=2q(K-j+1)\) while
   \(H_{2j-1}^{(c)}=q>0\) and \(H_{2j}^{(c)}=2q>0\).

These are countermodels for improving the **abstract inference**; they are **not** random closed two-form fields or counterexamples to the upstream rotation theorem.

## Concrete example (the \(q=4\) lattice-frequency setting)

Suppose \(K=10\) and the abstract assumptions hold, including
\[
M_{20}\le40.
\]
Quantization makes \(M_{20}=8B\) with \(0\le B\le5\), so the theorem yields
\[
M_0=M_1=\cdots=M_{10}=0.
\]

This is best possible from \(M_{20}\le40\) alone:
\[
M_n=4(n-10)_+
\]
satisfies the assumptions, has \(M_{20}=40\), and has \(M_{11}=4>0\).

The bound \(40\) is the largest admissible quantized value ensuring the first ten orders vanish. At \(M_{20}=48\), the model
\(M_n=4(n-8)_+\) has \(M_{10}=8>0\).

## What the statement does and does not offer

- It is a **finite-order certificate**: a bound at order \(2K\) can certify exact disappearance of a whole interval of lower orders, without an asymptotic limit.
- It still **presupposes** the nontrivial midpoint inequality and even-frequency quantization established in the surrounding argument. A raw Fourier degree estimate without them is insufficient.
- The hinge sequences demonstrate tightness only for the abstraction, not for realizable reflection-positive field theories.
- No claim of original mathematical novelty, completion of rotation invariance, construction of four-dimensional Yang–Mills, or a positive mass gap is made.

## Checks

Run the independent finite checks with Python standard library:

\`\`\`sh
python3 explorations/finite-order-zero-window/check_finite_window.py
\`\`\`

They exhaustively examine bounded quantized convex difference sequences for small \(K\) and verify midpoint inequalities for the hinge examples. The proof above, not the finite checks, establishes the all-order result.

## Possible next exploration

Can one derive *model-level* order-by-order bounds on \(M_{2K}\) with the required quantization and positivity, rather than merely assume them? That is where the abstract certificate would connect to the upstream analytic and field-theoretic work.
