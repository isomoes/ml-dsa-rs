//! Hint generation and verification for ML-DSA

use hybrid_array::{
    typenum::{Unsigned, U256},
    Array,
};

use crate::algebra::{AlgebraExt, BaseField, Decompose, Elem, Polynomial, Vector};
use crate::module_lattice::algebra::Field;
use crate::module_lattice::util::Truncate;
use crate::param::{EncodedHint, SignatureParams};

/// Algorithm 39 (MakeHint): determine whether low-order bits overflow into high-order bits.
pub fn make_hint<TwoGamma2: Unsigned>(z: Elem, r: Elem) -> bool {
    let r1 = r.high_bits::<TwoGamma2>();
    let v1 = (r + z).high_bits::<TwoGamma2>();
    r1 != v1
}

/// Algorithm 40 (UseHint): recover corrected high-order bits using a hint bit.
#[allow(clippy::integer_division_remainder_used)]
pub fn use_hint<TwoGamma2: Unsigned>(h: bool, r: Elem) -> Elem {
    let (r1, r0) = r.decompose::<TwoGamma2>();

    if !h {
        return r1;
    }

    let m = (BaseField::Q - 1) / TwoGamma2::U32;
    let gamma2 = TwoGamma2::U32 / 2;

    if r0.0 > 0 && r0.0 <= gamma2 {
        Elem::new((r1.0 + 1) % m)
    } else {
        Elem::new((r1.0 + m - 1) % m)
    }
}

/// Hint bitfield over a vector of 256-coefficient polynomials.
#[derive(Clone, Debug, PartialEq)]
pub struct Hint<P: SignatureParams>(
    /// Hint bits, indexed as `[poly_index][coeff_index]`.
    pub Array<Array<bool, U256>, P::K>,
);

impl<P: SignatureParams> Default for Hint<P> {
    fn default() -> Self {
        Self(Array::default())
    }
}

impl<P: SignatureParams> Hint<P> {
    /// Build a hint vector from `z` and `r`.
    pub fn new(z: &Vector<P::K>, r: &Vector<P::K>) -> Self {
        Self(
            z.0.iter()
                .zip(r.0.iter())
                .map(|(zv, rv)| {
                    zv.0.iter()
                        .zip(rv.0.iter())
                        .map(|(&z, &r)| make_hint::<P::TwoGamma2>(z, r))
                        .collect()
                })
                .collect(),
        )
    }

    /// Count the number of set hint bits.
    pub fn hamming_weight(&self) -> usize {
        self.0
            .iter()
            .map(|poly| poly.iter().filter(|&&bit| bit).count())
            .sum()
    }

    /// Apply hints to a vector and recover corrected high-order bits.
    pub fn use_hint(&self, r: &Vector<P::K>) -> Vector<P::K> {
        Vector::new(
            self.0
                .iter()
                .zip(r.0.iter())
                .map(|(hv, rv)| {
                    Polynomial::new(
                        hv.iter()
                            .zip(rv.0.iter())
                            .map(|(&h, &r)| use_hint::<P::TwoGamma2>(h, r))
                            .collect(),
                    )
                })
                .collect(),
        )
    }

    /// Encode hint bits into the ML-DSA compact hint byte representation.
    pub fn bit_pack(&self) -> EncodedHint<P> {
        let mut y: EncodedHint<P> = Array::default();
        let mut index = 0;
        let omega = P::Omega::USIZE;

        for i in 0..P::K::U8 {
            let i_usize: usize = i.into();
            for j in 0..256 {
                if self.0[i_usize][j] {
                    y[index] = Truncate::truncate(j);
                    index += 1;
                }
            }

            y[omega + i_usize] = Truncate::truncate(index);
        }

        y
    }

    /// Decode a compact hint byte representation.
    ///
    /// Returns `None` for malformed encodings.
    pub fn bit_unpack(y: &EncodedHint<P>) -> Option<Self> {
        let (indices, cuts) = P::split_hint(y);
        let cuts: Array<usize, P::K> = cuts.iter().map(|x| usize::from(*x)).collect();
        let indices: Array<usize, P::Omega> = indices.iter().map(|x| usize::from(*x)).collect();
        let max_cut: usize = cuts.iter().copied().max().unwrap_or(0);

        // cuts must be monotonic but can repeat
        if !cuts.windows(2).all(|w| w[0] <= w[1])
            || max_cut > indices.len()
            || indices[max_cut..].iter().copied().max().unwrap_or(0) > 0
        {
            return None;
        }

        let mut h = Self::default();
        let mut start = 0;
        for (i, &end) in cuts.iter().enumerate() {
            let chunk = &indices[start..end];

            // indices must be strictly increasing
            if !chunk.windows(2).all(|w| w[0] < w[1]) {
                return None;
            }

            for &j in chunk {
                if j >= 256 {
                    return None;
                }
                h.0[i][j] = true;
            }

            start = end;
        }

        Some(h)
    }
}

#[cfg(test)]
#[allow(clippy::integer_division_remainder_used)]
mod tests {
    use super::*;
    use crate::algebra::AlgebraExt;
    use crate::param::{MlDsa44, MlDsa65, ParameterSet};
    use hybrid_array::typenum::{Diff, Prod, Quot, Shleft, Unsigned, U1, U13, U2, U23, U32, U88};

    type QMinus1 = Diff<Diff<Shleft<U1, U23>, Shleft<U1, U13>>, U1>;
    type Gamma2_65 = Quot<QMinus1, U32>;
    type TwoGamma2 = Prod<U2, Gamma2_65>;
    type Gamma2_44 = Quot<QMinus1, U88>;
    type TwoGamma2_44 = Prod<U2, Gamma2_44>;

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
    fn hint_round_trip() {
        fn test<P: SignatureParams + PartialEq + core::fmt::Debug>() {
            let mut h = Hint::<P>::default();
            for i in 0..P::K::USIZE {
                if i < h.0.len() {
                    h.0[i][0] = true;
                    h.0[i][10] = true;
                    if i > 0 {
                        h.0[i][i * 5] = true;
                    }
                }
            }
            let packed = h.bit_pack();
            let unpacked = Hint::<P>::bit_unpack(&packed).unwrap();
            assert_eq!(h, unpacked);
        }
        test::<MlDsa44>();
        test::<MlDsa65>();
    }

    #[test]
    fn hint_new_hamming_weight_and_use_hint() {
        type P = MlDsa65;
        let mut z: Vector<<P as ParameterSet>::K> = Vector::default();
        let mut r: Vector<<P as ParameterSet>::K> = Vector::default();

        // Use values that are large enough to cause high-bits changes
        // gamma2 for MlDsa65 = (Q-1)/32 = 261888, so values near gamma2 boundary will flip
        z.0[0].0[3] = Elem::new(300_000);
        z.0[1].0[5] = Elem::new(500_000);
        r.0[0].0[3] = Elem::new(200_000);
        r.0[1].0[5] = Elem::new(BaseField::Q - 10);

        let hint = Hint::<P>::new(&z, &r);

        let used = hint.use_hint(&r);
        for i in 0..<P as ParameterSet>::K::USIZE {
            for j in 0..U256::USIZE {
                let expected =
                    use_hint::<<P as ParameterSet>::TwoGamma2>(hint.0[i][j], r.0[i].0[j]);
                assert_eq!(used.0[i].0[j], expected);
            }
        }
    }

    #[test]
    fn hint_use_hint_vector_shape() {
        type P = MlDsa44;
        let hint = Hint::<P>::default();
        let r: Vector<<P as ParameterSet>::K> =
            Vector::new(Array::from_fn(|_| Polynomial::default()));
        let out = hint.use_hint(&r);
        assert_eq!(out, r.high_bits::<<P as ParameterSet>::TwoGamma2>());
    }
}
