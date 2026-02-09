#import "@preview/cetz:0.2.2": canvas, draw

#set page(width: 900pt)
#set text(size: 11pt)

= Lattice Geometry for ML-DSA (Intuition First)

This note explains lattice geometry without heavy formulas first, then connects it to the short-vector problems used in lattice cryptography.

== 1. What is a lattice?

Think of a lattice as an infinite grid made from two or more step vectors.

- Start at the origin.
- You can move by integer multiples of basis steps `b_1, b_2, ..., b_m`.
- Every reachable point is a lattice point.

In symbols (just one line):

$L(B) = { Bz : z in ZZ^m }$

where `B` stores basis vectors as columns and `z` is an integer vector.

#figure(
  caption: [A lattice as integer combinations of two basis vectors.],
  canvas(length: 8cm, {
    import draw: *

    for x in range(-3, 4) {
      for y in range(-3, 4) {
        circle((x * 0.9, y * 0.7), radius: 0.04, fill: rgb("#7A8594"))
      }
    }

    line((0, 0), (1.8, 0.0), stroke: (paint: rgb("#2563EB"), thickness: 1.2pt))
    line((0, 0), (0.9, 1.4), stroke: (paint: rgb("#DC2626"), thickness: 1.2pt))

    content((1.9, -0.1), $b_1$)
    content((0.95, 1.45), $b_2$)
    content((0.1, -0.2), [0])
  })
)

== 2. Same lattice, different bases

One lattice can have many valid bases.

- A good basis: vectors are short and not too skewed.
- A bad basis: vectors are long and very skewed.

This matters because algorithms are much easier with good bases.

== 3. What is a short vector?

A short vector is a nonzero lattice point close to the origin.

- Length is usually measured with Euclidean norm `||v||_2`.
- The shortest nonzero length is `lambda_1(L)`.

The Shortest Vector Problem (SVP):

- Input: a basis of lattice `L`.
- Goal: find nonzero `v in L` with minimal length.

#figure(
  caption: [Shortest vector intuition: nonzero lattice point nearest to origin.],
  canvas(length: 8cm, {
    import draw: *

    for x in range(-3, 4) {
      for y in range(-3, 4) {
        circle((x * 0.9, y * 0.9), radius: 0.035, fill: rgb("#9CA3AF"))
      }
    }

    circle((0, 0), radius: 0.06, fill: rgb("#111827"))
    line((0, 0), (0.9, 0.0), stroke: (paint: rgb("#059669"), thickness: 1.8pt))
    content((0.95, 0.07), [short])
  })
)

Why SVP is hard (high level):

- Dimension is large.
- Search space over integer combinations is huge.
- Best known algorithms are still expensive at cryptographic sizes.

== 4. Short relation form (SIS style)

Many lattice schemes use a linear equation plus a shortness condition.

Given `A in ZZ_q^(n x m)`, find `x != 0` such that:

$A x = 0 mod q$

and `x` is short.

This is hard because linear constraints alone are easy, and shortness alone is easy, but satisfying both together is hard.

Toy example over `q = 17`:

- `A = [3 5 7]`
- Need `3x_1 + 5x_2 + 7x_3 = 0 mod 17`
- `x = (1, 2, 3)` works since `3 + 10 + 21 = 34 = 0 mod 17`

== 5. Module lattices (ML-DSA view)

ML-DSA uses module lattices instead of plain integer lattices.

- Ring: `R_q = F_q[X] / (X^n + 1)` with `n = 256`.
- Vectors are in `R_q^k`.
- Constraints look like matrix equations over `R_q` with short vectors.

Each polynomial gives `n` coefficients, so module vectors still correspond to large high-dimensional lattices. The ring structure gives fast arithmetic (NTT) and compact keys/signatures.

#figure(
  caption: [Module-lattice idea: polynomial blocks flatten into many coefficients.],
  [
    $R_q^k ->$ polynomial blocks $->$ coefficient vector in dimension $n k$
  ]
)

== 6. Practical takeaway

- Geometry gives the hardness intuition: find very short constrained vectors.
- Algebra gives efficiency: structured polynomial arithmetic.
- ML-DSA combines both.
