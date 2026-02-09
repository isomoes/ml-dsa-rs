#import "@preview/basic-document-props:0.1.0": simple-page
#import "@preview/cetz:0.3.4": canvas, draw

#show: simple-page.with(
  "isomo",
  "",
  middle-text: "ntt intro",
  date: true,
  numbering: true,
  supress-mail-link: false,
)

#set heading(numbering: "1.1")
#set math.equation(numbering: "(1)")
#set text(size: 11pt, lang: "en")
#show table: it => align(center, it)

= Number Theoretic Transform (NTT) for ML-DSA

The NTT converts polynomial multiplication from $O(n^2)$ to $O(n log n)$, making lattice-based cryptography practical.

== 1. Setting

ML-DSA works in the polynomial ring:

$ R_q = ZZ_q [X] \/ (X^256 + 1), quad q = 8380417 $

A polynomial $f in R_q$ has 256 coefficients: $f = sum_(i=0)^(255) f_i X^i$.

Multiplying two polynomials naively costs $O(n^2)$. The NTT reduces this to $O(n log n)$ by evaluating at special points, multiplying pointwise, then interpolating back.

== 2. Core idea

The NTT is the finite-field analogue of the FFT. It exploits a primitive root of unity $zeta = 1753$ in $ZZ_q$ satisfying $zeta^256 equiv -1 mod q$.

#figure(
  caption: [NTT converts ring multiplication to pointwise multiplication.],
  canvas(length: 1.2cm, {
    import draw: *

    // Ring domain
    rect((-2.2, -0.4), (0.0, 0.8), stroke: rgb("#2563EB"), radius: 0.08)
    content((-1.1, 0.6), text(size: 9pt, fill: rgb("#2563EB"), weight: "bold")[Ring $R_q$])
    content((-1.1, 0.2), text(size: 9pt)[$f dot g$])
    content((-1.1, -0.1), text(size: 9pt, fill: rgb("#6B7280"))[$O(n^2)$])

    // NTT domain
    rect((2.2, -0.4), (4.4, 0.8), stroke: rgb("#DC2626"), radius: 0.08)
    content((3.3, 0.6), text(size: 9pt, fill: rgb("#DC2626"), weight: "bold")[NTT $T_q$])
    content((3.3, 0.2), text(size: 9pt)[$hat(f) #sym.circle.small hat(g)$])
    content((3.3, -0.1), text(size: 9pt, fill: rgb("#6B7280"))[$O(n)$])

    // Arrows
    line((0.1, 0.55), (2.1, 0.55), stroke: (paint: rgb("#059669"), thickness: 1pt), mark: (end: ">"))
    content((1.1, 0.72), text(size: 8pt, fill: rgb("#059669"))[NTT])

    line((2.1, 0.05), (0.1, 0.05), stroke: (paint: rgb("#7C3AED"), thickness: 1pt), mark: (end: ">"))
    content((1.1, -0.12), text(size: 8pt, fill: rgb("#7C3AED"))[$"NTT"^(-1)$])
  }),
)

The key identity:

$ "NTT"^(-1)(hat(f) circle.small hat(g)) = f dot g in R_q $

where $circle.small$ denotes pointwise (coefficient-by-coefficient) multiplication.

== 3. Forward NTT (Algorithm 41)

The forward NTT uses the Cooley-Tukey butterfly. At each layer, pairs of elements are combined using a twiddle factor $zeta^("brv"(m))$:

$
  cases(
    w[j] & <- w[j] + zeta^("brv"(m)) dot w[j + ell],
    w[j + ell] & <- w[j] - zeta^("brv"(m)) dot w[j + ell]
  )
$

where $"brv"$ is the 8-bit bit-reversal permutation.

The transform proceeds through 8 layers ($log_2 256 = 8$), halving the stride $ell$ each time:

#figure(
  caption: [Butterfly structure: 8 layers, stride halves each step.],
  canvas(length: 1.2cm, {
    import draw: *

    let sp = 1.4  // horizontal spacing between layers

    // Layer labels
    for (i, l) in ((0, "128"), (1, "64"), (2, "32"), (3, "..."), (4, "4"), (5, "2"), (6, "1")).enumerate() {
      let x = i * sp
      content((x, 1.2), text(size: 8pt, fill: rgb("#6B7280"))[$ell=#l$])
    }

    // Signal lines
    for k in range(0, 8) {
      let y = -k * 0.3
      line((-0.4, y), (6 * sp + 0.4, y), stroke: (paint: rgb("#D1D5DB"), thickness: 0.5pt))
      content((-0.7, y), text(size: 7pt)[$w_#k$])
    }

    // Layer 1 butterflies (stride 128 → pairs 0-4)
    for pair in ((0, 4),) {
      let y0 = -pair.at(0) * 0.3
      let y1 = -pair.at(1) * 0.3
      line((0 * sp, y0), (0 * sp, y1), stroke: (paint: rgb("#2563EB"), thickness: 1pt))
      circle((0 * sp, y0), radius: 0.03, fill: rgb("#2563EB"))
      circle((0 * sp, y1), radius: 0.03, fill: rgb("#2563EB"))
    }

    // Layer 2 butterflies (stride 64 → pairs 0-2, 4-6)
    for pair in ((0, 2), (4, 6)) {
      let y0 = -pair.at(0) * 0.3
      let y1 = -pair.at(1) * 0.3
      line((1 * sp, y0), (1 * sp, y1), stroke: (paint: rgb("#2563EB"), thickness: 1pt))
      circle((1 * sp, y0), radius: 0.03, fill: rgb("#2563EB"))
      circle((1 * sp, y1), radius: 0.03, fill: rgb("#2563EB"))
    }

    // Layer 3 butterflies (stride 32 → adjacent pairs)
    for pair in ((0, 1), (2, 3), (4, 5), (6, 7)) {
      let y0 = -pair.at(0) * 0.3
      let y1 = -pair.at(1) * 0.3
      line((2 * sp, y0), (2 * sp, y1), stroke: (paint: rgb("#2563EB"), thickness: 1pt))
      circle((2 * sp, y0), radius: 0.03, fill: rgb("#2563EB"))
      circle((2 * sp, y1), radius: 0.03, fill: rgb("#2563EB"))
    }

    // Dots for middle layers (layers 4-6: strides 16, 8)
    content((3 * sp, -1.05), text(size: 10pt, fill: rgb("#9CA3AF"))[...])

    // Layer 7 butterflies (stride 4 → adjacent pairs)
    for pair in ((0, 1), (2, 3), (4, 5), (6, 7)) {
      let y0 = -pair.at(0) * 0.3
      let y1 = -pair.at(1) * 0.3
      line((4 * sp, y0), (4 * sp, y1), stroke: (paint: rgb("#7C3AED"), thickness: 1pt))
      circle((4 * sp, y0), radius: 0.03, fill: rgb("#7C3AED"))
      circle((4 * sp, y1), radius: 0.03, fill: rgb("#7C3AED"))
    }

    // Layer 8a butterflies (stride 2 → adjacent pairs)
    for pair in ((0, 1), (2, 3), (4, 5), (6, 7)) {
      let y0 = -pair.at(0) * 0.3
      let y1 = -pair.at(1) * 0.3
      line((5 * sp, y0), (5 * sp, y1), stroke: (paint: rgb("#DC2626"), thickness: 1pt))
      circle((5 * sp, y0), radius: 0.03, fill: rgb("#DC2626"))
      circle((5 * sp, y1), radius: 0.03, fill: rgb("#DC2626"))
    }

    // Layer 8b butterflies (stride 1 → adjacent pairs)
    for pair in ((0, 1), (2, 3), (4, 5), (6, 7)) {
      let y0 = -pair.at(0) * 0.3
      let y1 = -pair.at(1) * 0.3
      line((6 * sp, y0), (6 * sp, y1), stroke: (paint: rgb("#DC2626"), thickness: 1pt))
      circle((6 * sp, y0), radius: 0.03, fill: rgb("#DC2626"))
      circle((6 * sp, y1), radius: 0.03, fill: rgb("#DC2626"))
    }
  }),
)

Rust implementation (from `src/ntt.rs`):

```rust
fn ntt_layer<const LEN: usize, const ITERATIONS: usize>(
    w: &mut [Elem; 256], m: &mut usize,
) {
    for i in 0..ITERATIONS {
        let start = i * 2 * LEN;
        *m += 1;
        let z = ZETA_POW_BITREV[*m];  // twiddle factor
        for j in start..(start + LEN) {
            let t = z * w[j + LEN];
            w[j + LEN] = w[j] - t;    // butterfly subtract
            w[j]       = w[j] + t;    // butterfly add
        }
    }
}

// Full forward NTT: 8 layers, stride halves each time
fn ntt(f: &Polynomial) -> NttPolynomial {
    let mut w = f.coefficients();
    let mut m = 0;
    ntt_layer::<128, 1>(&mut w, &mut m);   // ℓ = 128
    ntt_layer::< 64, 2>(&mut w, &mut m);   // ℓ = 64
    ntt_layer::< 32, 4>(&mut w, &mut m);   // ℓ = 32
    ntt_layer::< 16, 8>(&mut w, &mut m);   // ℓ = 16
    ntt_layer::<  8,16>(&mut w, &mut m);   // ℓ = 8
    ntt_layer::<  4,32>(&mut w, &mut m);   // ℓ = 4
    ntt_layer::<  2,64>(&mut w, &mut m);   // ℓ = 2
    ntt_layer::<  1,128>(&mut w, &mut m);  // ℓ = 1
    NttPolynomial::new(w)
}
```

The const generics `LEN` (stride) and `ITERATIONS` (number of blocks) ensure all loop bounds are compile-time constants, avoiding runtime division.

== 4. Inverse NTT (Algorithm 42)

The inverse NTT uses the Gentleman-Sande butterfly, reversing the layer order and negating the twiddle factors:

$
  cases(
    w[j] & <- w[j] + w[j + ell],
    w[j + ell] & <- (-zeta^("brv"(m))) dot (w[j] - w[j + ell])
  )
$

After all 8 layers, every coefficient is scaled by $256^(-1) mod q = 8347681$:

$ f_i = 256^(-1) dot w_i mod q $

```rust
fn ntt_inverse(f_hat: &NttPolynomial) -> Polynomial {
    const INVERSE_256: Elem = Elem::new(8_347_681);
    let mut w = f_hat.coefficients();
    let mut m = 256;
    ntt_inverse_layer::<  1, 128>(&mut w, &mut m);  // ℓ = 1
    ntt_inverse_layer::<  2,  64>(&mut w, &mut m);  // ℓ = 2
    ntt_inverse_layer::<  4,  32>(&mut w, &mut m);  // ℓ = 4
    // ... layers 8, 16, 32, 64, 128
    ntt_inverse_layer::<128,   1>(&mut w, &mut m);  // ℓ = 128
    INVERSE_256 * Polynomial::new(w)
}
```

== 5. Pointwise multiplication (Algorithm 45)

In the NTT domain, polynomial multiplication reduces to pointwise multiplication of 256 coefficients:

$ hat(h)_i = hat(f)_i dot hat(g)_i mod q, quad i = 0, dots, 255 $

This works because the NTT decomposes $R_q$ into 256 copies of $ZZ_q$.

```rust
fn multiply_ntt(f_hat: &NttPolynomial, g_hat: &NttPolynomial)
    -> NttPolynomial
{
    // Simply multiply corresponding coefficients
    NttPolynomial::new(
        f_hat.iter().zip(g_hat.iter())
            .map(|(&x, &y)| x * y)
            .collect()
    )
}
```

== 6. Twiddle factors

The twiddle factors are precomputed in bit-reversed order:

$ "ZETA\_POW\_BITREV"[i] = zeta^("brv"_8 (i)) mod q $

where $"brv"_8$ reverses the 8 bits of $i$. This ordering matches the access pattern of the butterfly, so the NTT reads twiddle factors sequentially ($m = 1, 2, 3, dots$).

Key constants:
- $zeta = 1753$ (primitive 512th root of unity mod $q$)
- $zeta^(128) mod q = 4808194$ (first twiddle factor, FIPS 204 Appendix B)
- $256^(-1) mod q = 8347681$

== 7. Why NTT matters for ML-DSA

ML-DSA key generation, signing, and verification all require matrix-vector products over $R_q^(k times ell)$. Each entry involves polynomial multiplication.

Without NTT: each multiplication is $O(n^2) = O(65536)$ operations.

With NTT: store the matrix $bold(A)$ in NTT form permanently. Multiply via:

$ bold(A) dot bold(s) = "NTT"^(-1)(hat(bold(A)) circle.small "NTT"(bold(s))) $

Cost per multiplication: $O(n log n) = O(2048)$ --- a *32x speedup*.

#figure(
  caption: [NTT usage in ML-DSA: matrix $bold(A)$ stays in NTT domain.],
  canvas(length: 1.2cm, {
    import draw: *

    // Key gen box
    rect((-0.3, -0.3), (1.8, 0.9), stroke: rgb("#2563EB"), radius: 0.06)
    content((0.75, 0.7), text(size: 8pt, fill: rgb("#2563EB"), weight: "bold")[KeyGen])
    content((0.75, 0.35), text(size: 8pt)[$hat(bold(A)) <- "SampleNTT"$])
    content((0.75, 0.05), text(size: 8pt)[$bold(t) = "NTT"^(-1)(hat(bold(A)) dot hat(bold(s)))$])

    // Sign box
    rect((2.3, -0.3), (4.4, 0.9), stroke: rgb("#059669"), radius: 0.06)
    content((3.35, 0.7), text(size: 8pt, fill: rgb("#059669"), weight: "bold")[Sign])
    content((3.35, 0.35), text(size: 8pt)[$hat(bold(w)) = hat(bold(A)) dot hat(bold(y))$])
    content((3.35, 0.05), text(size: 8pt)[$bold(w)_1 = "HighBits"(...)$])

    // Verify box
    rect((4.9, -0.3), (7.2, 0.9), stroke: rgb("#DC2626"), radius: 0.06)
    content((6.05, 0.7), text(size: 8pt, fill: rgb("#DC2626"), weight: "bold")[Verify])
    content((6.05, 0.35), text(size: 8pt)[$hat(bold(A)) dot hat(bold(z)) - hat(bold(c)) dot hat(bold(t))$])
    content((6.05, 0.05), text(size: 8pt)[all in NTT domain])

    // Arrows
    line((1.9, 0.3), (2.2, 0.3), stroke: (paint: rgb("#6B7280"), thickness: 0.8pt), mark: (end: ">"))
    line((4.5, 0.3), (4.8, 0.3), stroke: (paint: rgb("#6B7280"), thickness: 0.8pt), mark: (end: ">"))
  }),
)

== 8. Summary

#table(
  columns: (auto, auto, auto),
  inset: 8pt,
  align: left,
  [*Operation*], [*Algorithm*], [*Complexity*],
  [Forward NTT], [Alg. 41 (Cooley-Tukey)], [$O(n log n)$],
  [Inverse NTT], [Alg. 42 (Gentleman-Sande)], [$O(n log n)$],
  [Pointwise multiply], [Alg. 45], [$O(n)$],
  [Naive poly multiply], [schoolbook], [$O(n^2)$],
)

The NTT is the computational backbone of ML-DSA: it makes polynomial arithmetic fast enough for real-world digital signatures.
