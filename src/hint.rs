//! Hint generation and verification for ML-DSA

use hybrid_array::{
    typenum::{Unsigned, U256},
    Array, ArraySize,
};

use crate::algebra::{AlgebraExt, BaseField, Decompose, Elem, Vector};
use crate::module_lattice::algebra::Field;

/// Algorithm 39 (MakeHint): determine whether low-order bits overflow into high-order bits.
pub fn make_hint<TwoGamma2: Unsigned>(z: Elem, r: Elem) -> bool {
    let r1 = r.high_bits::<TwoGamma2>();
    let v1 = (r + z).high_bits::<TwoGamma2>();
    r1 != v1
}

/// Algorithm 40 (UseHint): recover corrected high-order bits using a hint bit.
pub fn use_hint<TwoGamma2: Unsigned>(h: bool, r: Elem) -> Elem {
    let (r1, r0) = r.decompose::<TwoGamma2>();

    if !h {
        return r1;
    }

    let m = (BaseField::Q - 1) / TwoGamma2::U32;
    let gamma2 = TwoGamma2::U32 / 2;
    let r0_positive = r0.0 > 0 && r0.0 <= gamma2;

    if r0_positive {
        Elem::new((r1.0 + 1) % m)
    } else {
        Elem::new((r1.0 + m - 1) % m)
    }
}

/// Minimal parameter requirements for hint operations.
pub trait HintParams {
    /// Number of polynomials in vectors using this parameter set.
    type K: ArraySize;
    /// The value `2 * gamma2` as a type-level integer.
    type TwoGamma2: Unsigned;
    /// Maximum number of nonzero hint coefficients.
    const OMEGA: usize;
}

/// Hint bitfield over a vector of 256-coefficient polynomials.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Hint<P: HintParams> {
    /// Hint bits, indexed as `[poly_index][coeff_index]`.
    pub bits: Array<Array<bool, U256>, P::K>,
}

impl<P: HintParams> Hint<P> {
    /// Build a hint vector from `z` and `r`.
    pub fn new(z: &Vector<P::K>, r: &Vector<P::K>) -> Self {
        let bits = Array::from_fn(|i| {
            Array::from_fn(|j| make_hint::<P::TwoGamma2>(z.0[i].0[j], r.0[i].0[j]))
        });
        Self { bits }
    }

    /// Count the number of set hint bits.
    pub fn hamming_weight(&self) -> usize {
        self.bits
            .iter()
            .map(|poly| poly.iter().filter(|&&bit| bit).count())
            .sum()
    }

    /// Apply hints to a vector and recover corrected high-order bits.
    pub fn use_hint(&self, r: &Vector<P::K>) -> Vector<P::K> {
        Vector::new(Array::from_fn(|i| {
            crate::algebra::Polynomial::new(Array::from_fn(|j| {
                use_hint::<P::TwoGamma2>(self.bits[i][j], r.0[i].0[j])
            }))
        }))
    }

    /// Encode hint bits into the ML-DSA compact hint byte representation.
    ///
    /// Returns `None` if the hint exceeds `OMEGA` set bits.
    pub fn bit_pack(&self) -> Option<Vec<u8>> {
        let k = P::K::USIZE;
        let omega = P::OMEGA;

        if omega > usize::from(u8::MAX) {
            return None;
        }

        if self.hamming_weight() > omega {
            return None;
        }

        let mut out = vec![0u8; omega + k];
        let mut index = 0usize;

        for i in 0..k {
            for j in 0..U256::USIZE {
                if self.bits[i][j] {
                    #[allow(clippy::cast_possible_truncation, clippy::as_conversions)]
                    {
                        out[index] = j as u8;
                    }
                    index += 1;
                }
            }

            #[allow(clippy::cast_possible_truncation, clippy::as_conversions)]
            {
                out[omega + i] = index as u8;
            }
        }

        Some(out)
    }

    /// Decode a compact hint byte representation.
    ///
    /// Returns `None` for malformed encodings.
    pub fn bit_unpack(bytes: &[u8]) -> Option<Self> {
        let k = P::K::USIZE;
        let omega = P::OMEGA;

        if bytes.len() != omega + k {
            return None;
        }

        let mut bits = Array::from_fn(|_| Array::from_fn(|_| false));
        let mut index = 0usize;

        for i in 0..k {
            let next = usize::from(bytes[omega + i]);
            if next < index || next > omega {
                return None;
            }

            for j in index..next {
                let coeff = usize::from(bytes[j]);
                if coeff >= U256::USIZE {
                    return None;
                }

                if j > index && coeff <= usize::from(bytes[j - 1]) {
                    return None;
                }

                bits[i][coeff] = true;
            }

            index = next;
        }

        if bytes[index..omega].iter().any(|&x| x != 0) {
            return None;
        }

        Some(Self { bits })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::algebra::{AlgebraExt, Decompose, Polynomial};
    use crate::module_lattice::algebra::Field;
    use hybrid_array::typenum::{Diff, Prod, Quot, Shleft, U1, U13, U2, U23, U32, U4, U88};

    type QMinus1 = Diff<Diff<Shleft<U1, U23>, Shleft<U1, U13>>, U1>;
    type Gamma2_65 = Quot<QMinus1, U32>;
    type TwoGamma2 = Prod<U2, Gamma2_65>;
    type Gamma2_44 = Quot<QMinus1, U88>;
    type TwoGamma2_44 = Prod<U2, Gamma2_44>;

    #[derive(Debug, PartialEq, Eq)]
    struct TestParams;

    impl HintParams for TestParams {
        type K = U4;
        type TwoGamma2 = TwoGamma2;
        const OMEGA: usize = 8;
    }

    #[test]
    fn make_hint_matches_high_bits_change() {
        for &r in &[0, 1, 12345, 500000, BaseField::Q - 2, BaseField::Q - 1] {
            for &z in &[0, 1, 7, 1000, 50000, BaseField::Q - 10] {
                let rr = Elem::new(r);
                let zz = Elem::new(z);
                let expected = rr.high_bits::<TwoGamma2>() != (rr + zz).high_bits::<TwoGamma2>();
                assert_eq!(make_hint::<TwoGamma2>(zz, rr), expected);
            }
        }
    }

    #[test]
    fn use_hint_matches_reference_formula() {
        for &r in &[0, 1, 12345, 500000, BaseField::Q - 2, BaseField::Q - 1] {
            let rr = Elem::new(r);
            let (r1, r0) = rr.decompose::<TwoGamma2_44>();
            let m = (BaseField::Q - 1) / TwoGamma2_44::U32;
            let gamma2 = TwoGamma2_44::U32 / 2;

            let expected_h0 = r1;
            let expected_h1 = if r0.0 > 0 && r0.0 <= gamma2 {
                Elem::new((r1.0 + 1) % m)
            } else {
                Elem::new((r1.0 + m - 1) % m)
            };

            assert_eq!(use_hint::<TwoGamma2_44>(false, rr), expected_h0);
            assert_eq!(use_hint::<TwoGamma2_44>(true, rr), expected_h1);
        }
    }

    #[test]
    fn hint_new_hamming_weight_and_use_hint() {
        let mut z: Vector<U4> = Vector::default();
        let mut r: Vector<U4> = Vector::default();

        z.0[0].0[3] = Elem::new(1);
        z.0[1].0[5] = Elem::new(200_000);
        z.0[2].0[7] = Elem::new(100_000);
        r.0[2].0[7] = Elem::new(BaseField::Q - 10);

        let hint = Hint::<TestParams>::new(&z, &r);

        let expected_bits = Array::from_fn(|i| {
            Array::from_fn(|j| make_hint::<TwoGamma2>(z.0[i].0[j], r.0[i].0[j]))
        });
        assert_eq!(hint.bits, expected_bits);

        let expected_weight: usize = expected_bits
            .iter()
            .map(|poly| poly.iter().filter(|&&bit| bit).count())
            .sum();
        assert_eq!(hint.hamming_weight(), expected_weight);

        let used = hint.use_hint(&r);
        for i in 0..U4::USIZE {
            for j in 0..U256::USIZE {
                let expected = use_hint::<TwoGamma2>(hint.bits[i][j], r.0[i].0[j]);
                assert_eq!(used.0[i].0[j], expected);
            }
        }
    }

    #[test]
    fn hint_pack_unpack_round_trip() {
        let mut bits = Array::from_fn(|_| Array::from_fn(|_| false));
        bits[0][2] = true;
        bits[0][9] = true;
        bits[1][5] = true;
        bits[2][100] = true;
        bits[3][255] = true;

        let hint = Hint::<TestParams> { bits };
        let packed = hint.bit_pack().expect("valid hint should pack");
        let unpacked = Hint::<TestParams>::bit_unpack(&packed).expect("packed hint should decode");

        assert_eq!(unpacked, hint);
    }

    #[test]
    fn hint_unpack_rejects_unsorted_indices() {
        let mut bytes = vec![0u8; TestParams::OMEGA + U4::USIZE];
        bytes[0] = 7;
        bytes[1] = 3;
        bytes[TestParams::OMEGA] = 2;
        bytes[TestParams::OMEGA + 1] = 2;
        bytes[TestParams::OMEGA + 2] = 2;
        bytes[TestParams::OMEGA + 3] = 2;

        assert!(Hint::<TestParams>::bit_unpack(&bytes).is_none());
    }

    #[test]
    fn hint_unpack_rejects_nonzero_padding() {
        let mut bytes = vec![0u8; TestParams::OMEGA + U4::USIZE];
        bytes[0] = 10;
        bytes[TestParams::OMEGA] = 1;
        bytes[TestParams::OMEGA + 1] = 1;
        bytes[TestParams::OMEGA + 2] = 1;
        bytes[TestParams::OMEGA + 3] = 1;
        bytes[3] = 1;

        assert!(Hint::<TestParams>::bit_unpack(&bytes).is_none());
    }

    #[test]
    fn hint_pack_rejects_weight_above_omega() {
        let mut bits = Array::from_fn(|_| Array::from_fn(|_| false));
        for i in 0..(TestParams::OMEGA + 1) {
            bits[0][i] = true;
        }

        let hint = Hint::<TestParams> { bits };
        assert!(hint.bit_pack().is_none());
    }

    #[test]
    fn hint_use_hint_vector_shape() {
        let hint = Hint::<TestParams> {
            bits: Array::from_fn(|_| Array::from_fn(|_| false)),
        };
        let r: Vector<U4> = Vector::new(Array::from_fn(|_| Polynomial::default()));
        let out = hint.use_hint(&r);
        assert_eq!(out, r.high_bits::<TwoGamma2>());
    }
}
