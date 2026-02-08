//! Number Theoretic Transform (NTT) operations for ML-DSA
//!
//! Implements the forward NTT (Algorithm 41), inverse NTT (Algorithm 42), and
//! pointwise NTT multiplication (Algorithm 45) from FIPS 204.
//!
//! The NTT converts polynomials from the ring `R_q = Z_q[X]/(X^256 + 1)` to the
//! NTT algebra `T_q = Z_q^256`, where polynomial multiplication becomes pointwise
//! multiplication. This is the key optimization that makes lattice-based cryptography
//! practical, reducing multiplication from O(n^2) to O(n log n).

use crate::algebra::{BaseField, Elem, NttPolynomial, NttVector, Polynomial, Vector};
use crate::module_lattice::algebra::{Field, MultiplyNtt};
use hybrid_array::ArraySize;

// Since the powers of zeta used in the NTT and MultiplyNTTs are fixed, we use pre-computed tables
// to avoid the need to compute the exponentiations at runtime.
//
//   ZETA_POW_BITREV[i] = zeta^{BitRev_8(i)}
//
// Note that the const environment here imposes some annoying conditions.  Because operator
// overloading can't be const, we have to do all the reductions here manually.  Because `for` loops
// are forbidden in `const` functions, we do them manually with `while` loops.
//
// The values computed here match those provided in Appendix B of FIPS 204.
#[allow(clippy::cast_possible_truncation)]
#[allow(clippy::as_conversions)]
#[allow(clippy::integer_division_remainder_used, reason = "constant")]
const ZETA_POW_BITREV: [Elem; 256] = {
    const ZETA: u64 = 1753;
    const fn bitrev8(x: usize) -> usize {
        (x as u8).reverse_bits() as usize
    }

    // Compute the powers of zeta
    let mut pow = [Elem::new(0); 256];
    let mut i = 0;
    let mut curr = 1u64;
    while i < 256 {
        pow[i] = Elem::new(curr as u32);
        i += 1;
        curr = (curr * ZETA) % BaseField::QL;
    }

    // Reorder the powers according to bitrev8
    // Note that entry 0 is left as zero, in order to match the `zetas` array in the
    // specification.
    let mut pow_bitrev = [Elem::new(0); 256];
    let mut i = 1;
    while i < 256 {
        pow_bitrev[i] = pow[bitrev8(i)];
        i += 1;
    }
    pow_bitrev
};

/// Forward NTT transform trait.
///
/// Converts from the polynomial ring `R_q` to the NTT algebra `T_q`.
pub(crate) trait Ntt {
    /// The output type after NTT transformation
    type Output;
    /// Perform the forward NTT transform (Algorithm 41)
    fn ntt(&self) -> Self::Output;
}

/// Constant-time NTT butterfly layer.
///
/// Uses const generics to ensure loop bounds are compile-time constants,
/// avoiding UDIV instructions from runtime `step_by` calculations.
#[allow(clippy::inline_always)] // Required for constant-time guarantees in crypto code
#[inline(always)]
fn ntt_layer<const LEN: usize, const ITERATIONS: usize>(w: &mut [Elem; 256], m: &mut usize) {
    for i in 0..ITERATIONS {
        let start = i * 2 * LEN;
        *m += 1;
        let z = ZETA_POW_BITREV[*m];
        for j in start..(start + LEN) {
            let t = z * w[j + LEN];
            w[j + LEN] = w[j] - t;
            w[j] = w[j] + t;
        }
    }
}

impl Ntt for Polynomial {
    type Output = NttPolynomial;

    // Algorithm 41 NTT
    //
    // This implementation uses const-generic helper functions to ensure all loop
    // bounds are compile-time constants, avoiding potential UDIV instructions.
    fn ntt(&self) -> Self::Output {
        let mut w: [Elem; 256] = self.0.clone().into();
        let mut m = 0;

        ntt_layer::<128, 1>(&mut w, &mut m);
        ntt_layer::<64, 2>(&mut w, &mut m);
        ntt_layer::<32, 4>(&mut w, &mut m);
        ntt_layer::<16, 8>(&mut w, &mut m);
        ntt_layer::<8, 16>(&mut w, &mut m);
        ntt_layer::<4, 32>(&mut w, &mut m);
        ntt_layer::<2, 64>(&mut w, &mut m);
        ntt_layer::<1, 128>(&mut w, &mut m);

        NttPolynomial::new(w.into())
    }
}

impl<K: ArraySize> Ntt for Vector<K> {
    type Output = NttVector<K>;

    fn ntt(&self) -> Self::Output {
        NttVector::new(self.0.iter().map(Polynomial::ntt).collect())
    }
}

/// Inverse NTT transform trait.
///
/// Converts from the NTT algebra `T_q` back to the polynomial ring `R_q`.
#[allow(clippy::module_name_repetitions)]
pub(crate) trait NttInverse {
    /// The output type after inverse NTT transformation
    type Output;
    /// Perform the inverse NTT transform (Algorithm 42)
    fn ntt_inverse(&self) -> Self::Output;
}

/// Constant-time inverse NTT butterfly layer.
///
/// Uses const generics to ensure loop bounds are compile-time constants,
/// avoiding UDIV instructions from runtime `step_by` calculations.
#[allow(clippy::inline_always)] // Required for constant-time guarantees in crypto code
#[inline(always)]
fn ntt_inverse_layer<const LEN: usize, const ITERATIONS: usize>(
    w: &mut [Elem; 256],
    m: &mut usize,
) {
    for i in 0..ITERATIONS {
        let start = i * 2 * LEN;
        *m -= 1;
        let z = -ZETA_POW_BITREV[*m];
        for j in start..(start + LEN) {
            let t = w[j];
            w[j] = t + w[j + LEN];
            w[j + LEN] = z * (t - w[j + LEN]);
        }
    }
}

impl NttInverse for NttPolynomial {
    type Output = Polynomial;

    // Algorithm 42 NTT^{-1}
    //
    // This implementation uses const-generic helper functions to ensure all loop
    // bounds are compile-time constants, avoiding potential UDIV instructions.
    fn ntt_inverse(&self) -> Self::Output {
        /// The multiplicative inverse of 256 modulo Q: 256^{-1} mod 8380417 = 8347681
        const INVERSE_256: Elem = Elem::new(8_347_681);

        let mut w: [Elem; 256] = self.0.clone().into();
        let mut m = 256;

        ntt_inverse_layer::<1, 128>(&mut w, &mut m);
        ntt_inverse_layer::<2, 64>(&mut w, &mut m);
        ntt_inverse_layer::<4, 32>(&mut w, &mut m);
        ntt_inverse_layer::<8, 16>(&mut w, &mut m);
        ntt_inverse_layer::<16, 8>(&mut w, &mut m);
        ntt_inverse_layer::<32, 4>(&mut w, &mut m);
        ntt_inverse_layer::<64, 2>(&mut w, &mut m);
        ntt_inverse_layer::<128, 1>(&mut w, &mut m);

        INVERSE_256 * &Polynomial::new(w.into())
    }
}

impl<K: ArraySize> NttInverse for NttVector<K> {
    type Output = Vector<K>;

    fn ntt_inverse(&self) -> Self::Output {
        Vector::new(self.0.iter().map(NttPolynomial::ntt_inverse).collect())
    }
}

// Algorithm 45 MultiplyNTT
//
// For ML-DSA, the NTT multiplication is simply pointwise multiplication of the
// 256 coefficients. This is because the NTT decomposes the ring R_q into 256
// copies of Z_q.
impl MultiplyNtt for BaseField {
    fn multiply_ntt(lhs: &NttPolynomial, rhs: &NttPolynomial) -> NttPolynomial {
        NttPolynomial::new(
            lhs.0
                .iter()
                .zip(rhs.0.iter())
                .map(|(&x, &y)| x * y)
                .collect(),
        )
    }
}

#[cfg(test)]
#[allow(clippy::as_conversions)]
#[allow(clippy::cast_possible_truncation)]
mod test {
    use super::*;
    use crate::algebra::*;
    use hybrid_array::{
        typenum::{U2, U3},
        Array,
    };

    // Multiplication in R_q, modulo X^256 + 1
    fn poly_mul(lhs: &Polynomial, rhs: &Polynomial) -> Polynomial {
        let mut out = Polynomial::default();
        for (i, x) in lhs.0.iter().enumerate() {
            for (j, y) in rhs.0.iter().enumerate() {
                let (sign, index) = if i + j < 256 {
                    (Elem::new(1), i + j)
                } else {
                    (Elem::new(BaseField::Q - 1), i + j - 256)
                };

                out.0[index] = out.0[index] + (sign * *x * *y);
            }
        }
        out
    }

    // A polynomial with only a scalar component, to make simple test cases
    fn const_ntt(x: Int) -> NttPolynomial {
        let mut p = Polynomial::default();
        p.0[0] = Elem::new(x);
        p.ntt()
    }

    #[test]
    fn zeta_pow_bitrev_appendix_b() {
        // Verify first few entries from FIPS 204 Appendix B Table 2
        // zetas[0] is unused (left as 0)
        assert_eq!(ZETA_POW_BITREV[0].0, 0);
        // zetas[1] = zeta^{bitrev8(1)} = zeta^{128} mod Q
        // From FIPS 204 Appendix B: zetas[1] = 4808194
        assert_eq!(ZETA_POW_BITREV[1].0, 4_808_194);
    }

    #[test]
    fn inverse_256_correct() {
        // Verify that 256 * INVERSE_256 ≡ 1 (mod Q)
        let n = Elem::new(256);
        let inv = Elem::new(8_347_681);
        let product = n * inv;
        assert_eq!(product.0, 1);
    }

    #[test]
    fn ntt_round_trip() {
        // Verify that NTT and NTT^{-1} are actually inverses
        let f = Polynomial::new(Array::from_fn(|i| Elem::new(i as Int)));
        let f_hat = f.ntt();
        let f_unhat = f_hat.ntt_inverse();
        assert_eq!(f, f_unhat);
    }

    #[test]
    fn ntt_addition_homomorphism() {
        // Verify that NTT is a homomorphism with regard to addition
        let f = Polynomial::new(Array::from_fn(|i| Elem::new(i as Int)));
        let g = Polynomial::new(Array::from_fn(|i| Elem::new((2 * i) as Int)));
        let f_hat = f.ntt();
        let g_hat = g.ntt();

        let fg = &f + &g;
        let f_hat_g_hat = &f_hat + &g_hat;
        let fg_unhat = f_hat_g_hat.ntt_inverse();
        assert_eq!(fg, fg_unhat);
    }

    #[test]
    fn ntt_multiplication_homomorphism() {
        // Verify that NTT is a homomorphism with regard to multiplication
        let f = Polynomial::new(Array::from_fn(|i| Elem::new(i as Int)));
        let g = Polynomial::new(Array::from_fn(|i| Elem::new((2 * i) as Int)));
        let f_hat = f.ntt();
        let g_hat = g.ntt();

        let fg = poly_mul(&f, &g);
        let f_hat_g_hat = &f_hat * &g_hat;
        let fg_unhat = f_hat_g_hat.ntt_inverse();
        assert_eq!(fg, fg_unhat);
    }

    #[test]
    fn ntt_vector() {
        // Verify vector addition
        let v1: NttVector<U3> = NttVector::new(Array([const_ntt(1), const_ntt(1), const_ntt(1)]));
        let v2: NttVector<U3> = NttVector::new(Array([const_ntt(2), const_ntt(2), const_ntt(2)]));
        let v3: NttVector<U3> = NttVector::new(Array([const_ntt(3), const_ntt(3), const_ntt(3)]));
        assert_eq!((&v1 + &v2), v3);

        // Verify dot product
        assert_eq!((&v1 * &v2), const_ntt(6));
        assert_eq!((&v1 * &v3), const_ntt(9));
        assert_eq!((&v2 * &v3), const_ntt(18));
    }

    #[test]
    fn ntt_matrix() {
        use crate::algebra::NttMatrix;

        // Verify matrix multiplication by a vector
        let a: NttMatrix<U3, U2> = NttMatrix::new(Array([
            NttVector::new(Array([const_ntt(1), const_ntt(2)])),
            NttVector::new(Array([const_ntt(3), const_ntt(4)])),
            NttVector::new(Array([const_ntt(5), const_ntt(6)])),
        ]));
        let v_in: NttVector<U2> = NttVector::new(Array([const_ntt(1), const_ntt(2)]));
        let v_out: NttVector<U3> =
            NttVector::new(Array([const_ntt(5), const_ntt(11), const_ntt(17)]));
        assert_eq!(&a * &v_in, v_out);
    }

    #[test]
    fn ntt_vector_round_trip() {
        // Verify NTT round-trip for vectors
        let v: Vector<U3> = Vector::new(Array::from_fn(|k| {
            Polynomial::new(Array::from_fn(|i| {
                Elem::new(((k * 256 + i) % BaseField::Q as usize) as Int)
            }))
        }));
        let v_hat = v.ntt();
        let v_back = v_hat.ntt_inverse();
        assert_eq!(v, v_back);
    }

    #[test]
    fn ntt_zero_polynomial() {
        // NTT of zero should be zero
        let zero = Polynomial::default();
        let zero_hat = zero.ntt();
        assert_eq!(zero_hat, NttPolynomial::default());

        // And inverse should also give zero
        let zero_back = zero_hat.ntt_inverse();
        assert_eq!(zero_back, zero);
    }

    #[test]
    fn ntt_constant_polynomial() {
        // NTT of a constant polynomial c should give (c, c, ..., c)
        // because a constant is invariant under the NTT
        let c = Elem::new(42);
        let mut p = Polynomial::default();
        p.0[0] = c;
        let p_hat = p.ntt();

        // All coefficients in NTT domain should be 42
        for coeff in p_hat.0.iter() {
            assert_eq!(*coeff, c);
        }
    }

    #[test]
    fn multiply_ntt_commutativity() {
        let f = Polynomial::new(Array::from_fn(|i| Elem::new(i as Int)));
        let g = Polynomial::new(Array::from_fn(|i| Elem::new((3 * i + 1) as Int)));
        let f_hat = f.ntt();
        let g_hat = g.ntt();

        // Pointwise multiplication should be commutative
        let fg = &f_hat * &g_hat;
        let gf = &g_hat * &f_hat;
        assert_eq!(fg, gf);
    }
}
