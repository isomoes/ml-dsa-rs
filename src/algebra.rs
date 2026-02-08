//! Algebraic operations for ML-DSA
//!
//! This module defines the base field for ML-DSA (Q = 8,380,417) and provides
//! type aliases for the algebraic structures used throughout the implementation.
//! It also implements the ML-DSA-specific operations from FIPS 204:
//!
//! - Barrett modular reduction for arbitrary compile-time moduli
//! - Constant-time division for compile-time divisors
//! - Decompose (Algorithm 36)
//! - Power2Round (Algorithm 35)
//! - HighBits (Algorithm 37)
//! - LowBits (Algorithm 38)

use crate::module_lattice::algebra::Field;
use crate::module_lattice::util::Truncate;
use hybrid_array::{
    typenum::{Shleft, Unsigned, U1, U13},
    ArraySize,
};

crate::define_field!(BaseField, u32, u64, u128, 8_380_417);

/// The primitive integer type for the ML-DSA field
pub(crate) type Int = <BaseField as Field>::Int;

/// A field element in Z_q where q = 8,380,417
pub(crate) type Elem = crate::module_lattice::algebra::Elem<BaseField>;
/// A polynomial in R_q = Z_q[X] / (X^256 + 1)
pub(crate) type Polynomial = crate::module_lattice::algebra::Polynomial<BaseField>;
/// A vector of K polynomials from R_q
pub(crate) type Vector<K> = crate::module_lattice::algebra::Vector<BaseField, K>;
/// An NTT-domain polynomial in T_q
pub(crate) type NttPolynomial = crate::module_lattice::algebra::NttPolynomial<BaseField>;
/// A vector of K NTT-domain polynomials
pub(crate) type NttVector<K> = crate::module_lattice::algebra::NttVector<BaseField, K>;
/// A K x L matrix of NTT-domain polynomials
pub(crate) type NttMatrix<K, L> = crate::module_lattice::algebra::NttMatrix<BaseField, K, L>;

// We require modular reduction for three moduli: q, 2^d, and 2 * gamma2.  All three of these are
// greater than sqrt(q), which means that a number reduced mod q will always be less than M^2,
// which means that Barrett reduction will work.

/// Barrett modular reduction with precomputed constants for a compile-time modulus M.
///
/// This trait provides efficient modular reduction without expensive division,
/// using the identity: x mod M = x - floor(x * MULTIPLIER / 2^SHIFT) * M
pub(crate) trait BarrettReduce: Unsigned {
    /// Shift amount for Barrett reduction
    const SHIFT: usize;
    /// Precomputed multiplier: floor(2^SHIFT / M)
    const MULTIPLIER: u64;

    /// Reduce x modulo Self, where x < M^2
    fn reduce(x: u32) -> u32 {
        let m = Self::U64;
        let x: u64 = x.into();
        let quotient = (x * Self::MULTIPLIER) >> Self::SHIFT;
        let remainder = x - quotient * m;

        if remainder < m {
            Truncate::truncate(remainder)
        } else {
            Truncate::truncate(remainder - m)
        }
    }
}

impl<M> BarrettReduce for M
where
    M: Unsigned,
{
    #[allow(clippy::as_conversions)]
    const SHIFT: usize = 2 * (M::U64.ilog2() + 1) as usize;
    #[allow(clippy::integer_division_remainder_used)]
    const MULTIPLIER: u64 = (1 << Self::SHIFT) / M::U64;
}

/// Constant-time division by a compile-time constant divisor.
///
/// This trait provides a constant-time alternative to the hardware division
/// instruction, which has variable timing based on operand values.
/// Uses Barrett reduction to compute `x / M` where M is a compile-time constant.
pub(crate) trait ConstantTimeDiv: Unsigned {
    /// Bit shift for Barrett reduction, chosen to provide sufficient precision
    const CT_DIV_SHIFT: usize;
    /// Precomputed multiplier: ceil(2^SHIFT / M)
    const CT_DIV_MULTIPLIER: u64;

    /// Perform constant-time division of x by `Self::U32`
    /// Requires: x < Q (the field modulus, ~2^23)
    #[allow(clippy::inline_always)] // Required for constant-time guarantees in crypto code
    #[inline(always)]
    fn ct_div(x: u32) -> u32 {
        // Barrett reduction: q = (x * MULTIPLIER) >> SHIFT
        // This gives us floor(x / M) for x < 2^SHIFT / MULTIPLIER * M
        let x64 = u64::from(x);
        let quotient = (x64 * Self::CT_DIV_MULTIPLIER) >> Self::CT_DIV_SHIFT;
        // SAFETY: quotient is guaranteed to fit in u32 because:
        // - x < Q (~2^23), so quotient = x / M < x < 2^23 < 2^32
        #[allow(clippy::cast_possible_truncation, clippy::as_conversions)]
        let result = quotient as u32;
        result
    }
}

impl<M> ConstantTimeDiv for M
where
    M: Unsigned,
{
    // Use a shift that provides enough precision for the ML-DSA field (Q ~ 2^23)
    // We need SHIFT > log2(Q) + log2(M) to ensure accuracy
    // With Q < 2^24 and M < 2^20, SHIFT = 48 is sufficient
    const CT_DIV_SHIFT: usize = 48;

    // Precompute the multiplier at compile time
    // We add (M-1) before dividing to get ceiling division, ensuring we never underestimate
    #[allow(clippy::integer_division_remainder_used)]
    const CT_DIV_MULTIPLIER: u64 = (1u64 << Self::CT_DIV_SHIFT).div_ceil(M::U64);
}

/// Decompose a field element into high-order and low-order parts.
///
/// Implements Algorithm 36 (Decompose) from FIPS 204.
pub(crate) trait Decompose {
    /// Decompose self into (high_bits, low_bits) relative to modulus 2*gamma2
    fn decompose<TwoGamma2: Unsigned>(self) -> (Elem, Elem);
}

impl Decompose for Elem {
    // Algorithm 36 Decompose
    //
    // This implementation uses constant-time division to avoid timing side-channels.
    // The original algorithm used hardware division which has variable timing based
    // on operand values, potentially leaking secret information during signing.
    fn decompose<TwoGamma2: Unsigned>(self) -> (Elem, Elem) {
        let r_plus = self;
        let r0 = r_plus.mod_plus_minus::<TwoGamma2>();

        if r_plus - r0 == Elem::new(BaseField::Q - 1) {
            (Elem::new(0), r0 - Elem::new(1))
        } else {
            let diff = r_plus - r0;
            // Use constant-time division instead of hardware division
            let r1 = Elem::new(TwoGamma2::ct_div(diff.0));
            (r1, r0)
        }
    }
}

/// Extension trait for ML-DSA-specific algebraic operations.
///
/// Provides centered modular reduction, infinity norm computation,
/// and the Power2Round/HighBits/LowBits algorithms from FIPS 204.
#[allow(clippy::module_name_repetitions)]
pub(crate) trait AlgebraExt: Sized {
    /// Centered modular reduction: maps to (-M/2, M/2]
    fn mod_plus_minus<M: Unsigned>(&self) -> Self;
    /// Infinity norm (max absolute value under centered representation)
    fn infinity_norm(&self) -> Int;
    /// Algorithm 35 (Power2Round): split into high and low 2^d parts
    fn power2round(&self) -> (Self, Self);
    /// Algorithm 37 (HighBits): extract high-order bits relative to 2*gamma2
    fn high_bits<TwoGamma2: Unsigned>(&self) -> Self;
    /// Algorithm 38 (LowBits): extract low-order bits relative to 2*gamma2
    fn low_bits<TwoGamma2: Unsigned>(&self) -> Self;
}

impl AlgebraExt for Elem {
    fn mod_plus_minus<M: Unsigned>(&self) -> Self {
        let raw_mod = Elem::new(M::reduce(self.0));
        if raw_mod.0 <= M::U32 >> 1 {
            raw_mod
        } else {
            raw_mod - Elem::new(M::U32)
        }
    }

    // FIPS 204 defines the infinity norm differently for signed vs. unsigned integers:
    //
    // * For w in Z, |w|_\infinity = |w|, the absolute value of w
    // * For w in Z_q, |W|_infinity = |w mod^\pm q|
    //
    // Note that these two definitions are equivalent if |w| < q/2.  This property holds for all of
    // the signed integers used in this crate, so we can safely use the unsigned version.  However,
    // since mod_plus_minus is also unsigned, we need to unwrap the "negative" values.
    fn infinity_norm(&self) -> u32 {
        if self.0 <= BaseField::Q >> 1 {
            self.0
        } else {
            BaseField::Q - self.0
        }
    }

    // Algorithm 35 Power2Round
    //
    // In the specification, this function maps to signed integers rather than modular integers.
    // To avoid the need for a whole separate type for signed integer polynomials, we represent
    // these values using integers mod Q.  This is safe because Q is much larger than 2^13, so
    // there's no risk of overlap between positive numbers (x) and negative numbers (Q-x).
    fn power2round(&self) -> (Self, Self) {
        type D = U13;
        type Pow2D = Shleft<U1, D>;

        let r_plus = *self;
        let r0 = r_plus.mod_plus_minus::<Pow2D>();
        let r1 = Elem::new((r_plus - r0).0 >> D::USIZE);

        (r1, r0)
    }

    // Algorithm 37 HighBits
    fn high_bits<TwoGamma2: Unsigned>(&self) -> Self {
        self.decompose::<TwoGamma2>().0
    }

    // Algorithm 38 LowBits
    fn low_bits<TwoGamma2: Unsigned>(&self) -> Self {
        self.decompose::<TwoGamma2>().1
    }
}

impl AlgebraExt for Polynomial {
    fn mod_plus_minus<M: Unsigned>(&self) -> Self {
        Self(self.0.iter().map(AlgebraExt::mod_plus_minus::<M>).collect())
    }

    fn infinity_norm(&self) -> u32 {
        self.0.iter().map(AlgebraExt::infinity_norm).max().unwrap()
    }

    fn power2round(&self) -> (Self, Self) {
        let mut r1 = Self::default();
        let mut r0 = Self::default();

        for (i, x) in self.0.iter().enumerate() {
            (r1.0[i], r0.0[i]) = x.power2round();
        }

        (r1, r0)
    }

    fn high_bits<TwoGamma2: Unsigned>(&self) -> Self {
        Self(
            self.0
                .iter()
                .map(AlgebraExt::high_bits::<TwoGamma2>)
                .collect(),
        )
    }

    fn low_bits<TwoGamma2: Unsigned>(&self) -> Self {
        Self(
            self.0
                .iter()
                .map(AlgebraExt::low_bits::<TwoGamma2>)
                .collect(),
        )
    }
}

impl<K: ArraySize> AlgebraExt for Vector<K> {
    fn mod_plus_minus<M: Unsigned>(&self) -> Self {
        Self(self.0.iter().map(AlgebraExt::mod_plus_minus::<M>).collect())
    }

    fn infinity_norm(&self) -> u32 {
        self.0.iter().map(AlgebraExt::infinity_norm).max().unwrap()
    }

    fn power2round(&self) -> (Self, Self) {
        let mut r1 = Self::default();
        let mut r0 = Self::default();

        for (i, x) in self.0.iter().enumerate() {
            (r1.0[i], r0.0[i]) = x.power2round();
        }

        (r1, r0)
    }

    fn high_bits<TwoGamma2: Unsigned>(&self) -> Self {
        Self(
            self.0
                .iter()
                .map(AlgebraExt::high_bits::<TwoGamma2>)
                .collect(),
        )
    }

    fn low_bits<TwoGamma2: Unsigned>(&self) -> Self {
        Self(
            self.0
                .iter()
                .map(AlgebraExt::low_bits::<TwoGamma2>)
                .collect(),
        )
    }
}

#[cfg(test)]
#[allow(clippy::integer_division_remainder_used, reason = "tests")]
mod test {
    use super::*;
    use hybrid_array::typenum::{Prod, Quot, U2, U32, U4};

    // ML-DSA-65 parameters for testing: Gamma2 = (Q-1)/32, TwoGamma2 = 2*Gamma2
    type QMinus1 = hybrid_array::typenum::Diff<
        hybrid_array::typenum::Diff<Shleft<U1, hybrid_array::typenum::U23>, Shleft<U1, U13>>,
        U1,
    >;

    // Use simpler approach: define TwoGamma2 directly as a typenum
    // For ML-DSA-65: Gamma2 = (Q-1)/32 = 261888, TwoGamma2 = 523776
    type Gamma2 = Quot<QMinus1, U32>;
    type TwoGamma2 = Prod<U2, Gamma2>;

    const MOD: u32 = TwoGamma2::U32;
    const MOD_ELEM: Elem = Elem::new(MOD);

    #[test]
    fn base_field_constants() {
        assert_eq!(BaseField::Q, 8_380_417);
        // Verify Q = 2^23 - 2^13 + 1
        assert_eq!(BaseField::Q, (1 << 23) - (1 << 13) + 1);
    }

    #[test]
    fn mod_plus_minus() {
        for x in 0..MOD {
            let x = Elem::new(x);
            let x0 = x.mod_plus_minus::<TwoGamma2>();

            // Outputs from mod+- should be in the half-open interval (-gamma2, gamma2]
            let positive_bound = x0.0 <= MOD / 2;
            let negative_bound = x0.0 > BaseField::Q - MOD / 2;
            assert!(positive_bound || negative_bound);

            // The output should be equivalent to the input, mod 2 * gamma2.  We add 2 * gamma2
            // before comparing so that both values are "positive", avoiding interactions between
            // the mod-Q and mod-M operations.
            let xn = x + MOD_ELEM;
            let x0n = x0 + MOD_ELEM;
            assert_eq!(xn.0 % MOD, x0n.0 % MOD);
        }
    }

    #[test]
    fn decompose() {
        for x in 0..MOD {
            let x = Elem::new(x);
            let (x1, x0) = x.decompose::<TwoGamma2>();

            // The low-order output from decompose() is a mod+- output, optionally minus one.  So
            // they should be in the closed interval [-gamma2, gamma2].
            let positive_bound = x0.0 <= MOD / 2;
            let negative_bound = x0.0 >= BaseField::Q - MOD / 2;
            assert!(positive_bound || negative_bound);

            // The low-order and high-order outputs should combine to form the input.
            let xx = (MOD * x1.0 + x0.0) % BaseField::Q;
            assert_eq!(xx, x.0);
        }
    }

    #[test]
    fn barrett_reduce_boundary() {
        let m_minus_1 = TwoGamma2::U32 - 1;
        assert_eq!(TwoGamma2::reduce(m_minus_1), m_minus_1);
        assert_eq!(TwoGamma2::reduce(TwoGamma2::U32), 0);
        assert_eq!(TwoGamma2::reduce(TwoGamma2::U32 + 1), 1);
        assert_eq!(TwoGamma2::reduce(2 * TwoGamma2::U32 - 1), m_minus_1);
        assert_eq!(TwoGamma2::reduce(2 * TwoGamma2::U32), 0);
    }

    #[test]
    fn constant_time_div_accuracy() {
        for x in 0..1000 {
            assert_eq!(TwoGamma2::ct_div(x), x / TwoGamma2::U32);
        }
        for x in (BaseField::Q - 1000)..BaseField::Q {
            assert_eq!(TwoGamma2::ct_div(x), x / TwoGamma2::U32);
        }
    }

    #[test]
    fn decompose_edge_case() {
        let q_minus_1 = Elem::new(BaseField::Q - 1);
        let (r1, r0) = q_minus_1.decompose::<TwoGamma2>();
        let reconstructed = (MOD * r1.0 + r0.0) % BaseField::Q;
        assert_eq!(reconstructed, q_minus_1.0);
    }

    #[test]
    fn high_low_bits_consistency() {
        for x in [0, 1, MOD / 2, MOD - 1, MOD, MOD + 1, BaseField::Q - 1] {
            let elem = Elem::new(x);
            let (decomp_high, decomp_low) = elem.decompose::<TwoGamma2>();
            assert_eq!(elem.high_bits::<TwoGamma2>(), decomp_high);
            assert_eq!(elem.low_bits::<TwoGamma2>(), decomp_low);
        }
    }

    #[test]
    fn power2round_reconstruction() {
        // Power2Round should satisfy: r = r1 * 2^d + r0
        for x in [
            0,
            1,
            4096,
            8191,
            8192,
            8193,
            BaseField::Q / 2,
            BaseField::Q - 1,
        ] {
            let elem = Elem::new(x);
            let (r1, r0) = elem.power2round();

            // r0 should be in (-2^(d-1), 2^(d-1)]
            let d = 13;
            let half = 1u32 << (d - 1); // 4096
            let positive_bound = r0.0 <= half;
            let negative_bound = r0.0 > BaseField::Q - half;
            assert!(
                positive_bound || negative_bound,
                "r0={} out of range for x={}",
                r0.0,
                x
            );

            // Reconstruction: r1 * 2^d + r0 ≡ r (mod Q)
            let reconstructed = (r1.0 * (1 << d) + r0.0) % BaseField::Q;
            assert_eq!(
                reconstructed, x,
                "Power2Round reconstruction failed for x={}",
                x
            );
        }
    }

    #[test]
    fn infinity_norm_elem() {
        // Small positive values
        assert_eq!(Elem::new(0).infinity_norm(), 0);
        assert_eq!(Elem::new(1).infinity_norm(), 1);
        assert_eq!(Elem::new(100).infinity_norm(), 100);

        // Values near Q/2 boundary
        let half_q = BaseField::Q / 2;
        assert_eq!(Elem::new(half_q).infinity_norm(), half_q);

        // "Negative" values (close to Q)
        assert_eq!(Elem::new(BaseField::Q - 1).infinity_norm(), 1);
        assert_eq!(Elem::new(BaseField::Q - 100).infinity_norm(), 100);
    }

    #[test]
    fn infinity_norm_polynomial() {
        let mut coeffs = hybrid_array::Array::default();
        coeffs[0] = Elem::new(5);
        coeffs[1] = Elem::new(BaseField::Q - 10); // represents -10
        coeffs[2] = Elem::new(3);
        let poly = Polynomial::new(coeffs);

        assert_eq!(poly.infinity_norm(), 10);
    }

    #[test]
    fn infinity_norm_vector() {
        let mut coeffs1 = hybrid_array::Array::default();
        coeffs1[0] = Elem::new(5);
        let poly1 = Polynomial::new(coeffs1);

        let mut coeffs2 = hybrid_array::Array::default();
        coeffs2[0] = Elem::new(BaseField::Q - 20); // represents -20
        let poly2 = Polynomial::new(coeffs2);

        let vec: Vector<U4> = Vector::new(hybrid_array::Array::from_fn(|i| {
            if i == 0 {
                poly1.clone()
            } else if i == 1 {
                poly2.clone()
            } else {
                Polynomial::default()
            }
        }));

        assert_eq!(vec.infinity_norm(), 20);
    }

    #[test]
    fn mod_plus_minus_with_pow2() {
        // Test with 2^13 = 8192 (used in Power2Round)
        type Pow2D = Shleft<U1, U13>;
        let half = Pow2D::U32 / 2; // 4096

        // Value in positive range
        let x = Elem::new(100);
        let result = x.mod_plus_minus::<Pow2D>();
        assert_eq!(result.0, 100);

        // Value that wraps to negative
        let x = Elem::new(5000); // 5000 mod 8192 = 5000 > 4096, so negative
        let result = x.mod_plus_minus::<Pow2D>();
        // 5000 mod 8192 = 5000, which is > 4096, so result = 5000 - 8192 = -3192 ≡ Q - 3192
        assert_eq!(result.0, BaseField::Q - 3192);

        // Boundary: exactly half
        let x = Elem::new(half);
        let result = x.mod_plus_minus::<Pow2D>();
        assert_eq!(result.0, half); // half is included in positive range
    }

    #[test]
    fn decompose_all_ml_dsa_44_gamma2() {
        // ML-DSA-44: Gamma2 = (Q-1)/88 = 95232, TwoGamma2 = 190464
        type Gamma2_44 = Quot<QMinus1, hybrid_array::typenum::U88>;
        type TwoGamma2_44 = Prod<U2, Gamma2_44>;
        let mod_44 = TwoGamma2_44::U32;

        // Test a range of values
        for x in (0..mod_44).step_by(1000) {
            let elem = Elem::new(x);
            let (r1, r0) = elem.decompose::<TwoGamma2_44>();
            let reconstructed = (mod_44 * r1.0 + r0.0) % BaseField::Q;
            assert_eq!(
                reconstructed, x,
                "Decompose failed for x={} with ML-DSA-44 params",
                x
            );
        }
    }
}
