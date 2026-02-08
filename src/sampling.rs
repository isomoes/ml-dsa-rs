//! Sampling operations for ML-DSA.

use crate::{
    algebra::{BaseField, Elem, Int, NttMatrix, NttPolynomial, NttVector, Polynomial, Vector},
    crypto::{G, H},
};
use hybrid_array::{Array, ArraySize};

use crate::module_lattice::{algebra::Field, util::Truncate};

/// Eta parameter used by bounded sampling.
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub(crate) enum Eta {
    /// ETA = 2.
    Two,
    /// ETA = 4.
    Four,
}

/// Sampling-size adapter for Algorithm 34 (ExpandMask).
pub(crate) trait MaskSamplingSize {
    /// Number of bytes to squeeze for one polynomial.
    const SAMPLE_SIZE: usize;

    /// Unpack a squeezed byte slice into one polynomial.
    fn unpack(v: &[u8]) -> Polynomial;
}

/// Gamma1 = 2^17 mask unpacker.
pub(crate) struct Gamma1Power17;

/// Gamma1 = 2^19 mask unpacker.
pub(crate) struct Gamma1Power19;

impl MaskSamplingSize for Gamma1Power17 {
    const SAMPLE_SIZE: usize = 576;

    fn unpack(v: &[u8]) -> Polynomial {
        unpack_mask_poly(v, 18, 1 << 17)
    }
}

impl MaskSamplingSize for Gamma1Power19 {
    const SAMPLE_SIZE: usize = 640;

    fn unpack(v: &[u8]) -> Polynomial {
        unpack_mask_poly(v, 20, 1 << 19)
    }
}

// Algorithm 13 BytesToBits.
fn bit_set(z: &[u8], i: usize) -> bool {
    let bit_index = i & 0x07;
    let byte_index = i >> 3;

    z[byte_index] & (1 << bit_index) != 0
}

// Algorithm 14 CoeffFromThreeBytes.
fn coeff_from_three_bytes(b: [u8; 3]) -> Option<Elem> {
    let b0: Int = b[0].into();
    let b1: Int = b[1].into();
    let b2: Int = b[2].into();

    let b2p = if b2 > 127 { b2 - 128 } else { b2 };

    let z = (b2p << 16) + (b1 << 8) + b0;
    (z < BaseField::Q).then_some(Elem::new(z))
}

// Algorithm 15 CoeffFromHalfByte.
fn coeff_from_half_byte(b: u8, eta: Eta) -> Option<Elem> {
    match eta {
        Eta::Two if b < 15 => {
            let b = Int::from(match b {
                x if x < 5 => x,
                x if x < 10 => x - 5,
                _ => b - 10,
            });

            if b <= 2 {
                Some(Elem::new(2 - b))
            } else {
                Some(-Elem::new(b - 2))
            }
        }
        Eta::Four if b < 9 => {
            let b = Int::from(b);
            if b <= 4 {
                Some(Elem::new(4 - b))
            } else {
                Some(-Elem::new(b - 4))
            }
        }
        _ => None,
    }
}

fn coeffs_from_byte(z: u8, eta: Eta) -> (Option<Elem>, Option<Elem>) {
    (
        coeff_from_half_byte(z & 0x0F, eta),
        coeff_from_half_byte(z >> 4, eta),
    )
}

fn read_bits_le(v: &[u8], bit_offset: usize, width: usize) -> u32 {
    let mut out = 0u32;
    for i in 0..width {
        if bit_set(v, bit_offset + i) {
            out |= 1 << i;
        }
    }

    out
}

fn signed_to_elem(x: i32) -> Elem {
    if x >= 0 {
        Elem::new(u32::try_from(x).expect("non-negative i32 fits in u32"))
    } else {
        let abs = x.unsigned_abs();
        Elem::new(BaseField::Q - abs)
    }
}

fn unpack_mask_poly(v: &[u8], bits_per_coeff: usize, gamma1: u32) -> Polynomial {
    let mut out = Polynomial::default();
    for i in 0..256 {
        let t = read_bits_le(v, i * bits_per_coeff, bits_per_coeff);
        let c = i32::try_from(gamma1).expect("gamma1 fits i32")
            - i32::try_from(t).expect("packed coeff fits i32");
        out.0[i] = signed_to_elem(c);
    }

    out
}

// Algorithm 29 SampleInBall.
pub(crate) fn sample_in_ball(rho: &[u8], tau: usize) -> Polynomial {
    assert!(tau <= 256, "tau must be <= 256");

    const ONE: Elem = Elem::new(1);
    const MINUS_ONE: Elem = Elem::new(BaseField::Q - 1);

    let mut c = Polynomial::default();
    let mut ctx = H::default();
    ctx.absorb(rho);

    let mut s = [0u8; 8];
    ctx.squeeze(&mut s);

    let mut j = [0u8; 1];
    for i in (256 - tau)..256 {
        ctx.squeeze(&mut j);
        while usize::from(j[0]) > i {
            ctx.squeeze(&mut j);
        }

        let j = usize::from(j[0]);
        c.0[i] = c.0[j];
        c.0[j] = if bit_set(&s, i + tau - 256) {
            MINUS_ONE
        } else {
            ONE
        };
    }

    c
}

// Algorithm 30 RejNTTPoly.
fn rej_ntt_poly(rho: &[u8], r: u8, s: u8) -> NttPolynomial {
    let mut j = 0usize;

    let mut ctx = G::default();
    ctx.absorb(rho);
    ctx.absorb(&[s]);
    ctx.absorb(&[r]);

    let mut a = NttPolynomial::default();
    let mut sample = [0u8; 3];
    while j < 256 {
        ctx.squeeze(&mut sample);
        if let Some(x) = coeff_from_three_bytes(sample) {
            a.0[j] = x;
            j += 1;
        }
    }

    a
}

// Algorithm 31 RejBoundedPoly.
fn rej_bounded_poly(rho: &[u8], eta: Eta, r: u16) -> Polynomial {
    let mut j = 0usize;

    let mut ctx = H::default();
    ctx.absorb(rho);
    ctx.absorb(&r.to_le_bytes());

    let mut a = Polynomial::default();
    let mut z = [0u8; 1];
    while j < 256 {
        ctx.squeeze(&mut z);
        let (z0, z1) = coeffs_from_byte(z[0], eta);

        if let Some(coeff) = z0 {
            a.0[j] = coeff;
            j += 1;
        }

        if j == 256 {
            break;
        }

        if let Some(coeff) = z1 {
            a.0[j] = coeff;
            j += 1;
        }
    }

    a
}

// Algorithm 32 ExpandA.
pub(crate) fn expand_a<K: ArraySize, L: ArraySize>(rho: &[u8]) -> NttMatrix<K, L> {
    NttMatrix::new(Array::from_fn(|r| {
        NttVector::new(Array::from_fn(|s| {
            rej_ntt_poly(rho, u8::truncate(r), u8::truncate(s))
        }))
    }))
}

// Algorithm 33 ExpandS.
pub(crate) fn expand_s<K: ArraySize>(rho: &[u8], eta: Eta, base: usize) -> Vector<K> {
    Vector::new(Array::from_fn(|r| {
        let index = u16::truncate(r + base);
        rej_bounded_poly(rho, eta, index)
    }))
}

// Algorithm 34 ExpandMask.
pub(crate) fn expand_mask<K, Gamma1>(rho: &[u8], mu: u16) -> Vector<K>
where
    K: ArraySize,
    Gamma1: MaskSamplingSize,
{
    Vector::new(Array::from_fn(|r| {
        let mut ctx = H::default();
        ctx.absorb(rho);

        let r: u16 = u16::truncate(r);
        ctx.absorb(&(mu + r).to_le_bytes());

        let mut v = vec![0u8; Gamma1::SAMPLE_SIZE];
        ctx.squeeze(&mut v);

        Gamma1::unpack(&v)
    }))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::algebra::AlgebraExt;
    use hybrid_array::typenum::{U3, U4};

    fn hamming_weight(p: &Polynomial) -> usize {
        p.0.iter().filter(|x| x.0 != 0).count()
    }

    fn all_abs_le_eta(p: &Polynomial, eta: u32) -> bool {
        p.0.iter().all(|x| x.infinity_norm() <= eta)
    }

    fn all_mask_in_range(p: &Polynomial, gamma1: u32) -> bool {
        p.0.iter().all(|x| x.infinity_norm() <= gamma1)
    }

    #[test]
    fn bytes_to_bits_examples() {
        let data = [0b1010_0101u8, 0b0000_0001u8];

        assert!(bit_set(&data, 0));
        assert!(!bit_set(&data, 1));
        assert!(bit_set(&data, 2));
        assert!(bit_set(&data, 5));
        assert!(bit_set(&data, 8));
        assert!(!bit_set(&data, 9));
    }

    #[test]
    fn coeff_from_three_bytes_range() {
        assert_eq!(coeff_from_three_bytes([0, 0, 0]), Some(Elem::new(0)));
        assert_eq!(coeff_from_three_bytes([0, 0, 128]), Some(Elem::new(0)));

        let invalid = coeff_from_three_bytes([0xFF, 0xFF, 0x7F]);
        assert!(invalid.is_none());
    }

    #[test]
    fn coeff_from_half_byte_eta_two_and_four() {
        assert_eq!(coeff_from_half_byte(0, Eta::Two), Some(Elem::new(2)));
        assert_eq!(coeff_from_half_byte(2, Eta::Two), Some(Elem::new(0)));
        assert_eq!(
            coeff_from_half_byte(3, Eta::Two),
            Some(Elem::new(BaseField::Q - 1))
        );
        assert!(coeff_from_half_byte(15, Eta::Two).is_none());

        assert_eq!(coeff_from_half_byte(0, Eta::Four), Some(Elem::new(4)));
        assert_eq!(coeff_from_half_byte(4, Eta::Four), Some(Elem::new(0)));
        assert_eq!(
            coeff_from_half_byte(5, Eta::Four),
            Some(Elem::new(BaseField::Q - 1))
        );
        assert!(coeff_from_half_byte(9, Eta::Four).is_none());
    }

    #[test]
    fn sample_in_ball_has_expected_shape() {
        let rho = [7u8; 32];
        let tau = 39;
        let p = sample_in_ball(&rho, tau);

        assert_eq!(hamming_weight(&p), tau);
        let all_valid =
            p.0.iter()
                .all(|x| x.0 == 0 || x.0 == 1 || x.0 == BaseField::Q - 1);
        assert!(all_valid);
    }

    #[test]
    fn rej_ntt_poly_outputs_in_field_range() {
        let rho = [11u8; 32];
        let poly = rej_ntt_poly(&rho, 1, 2);
        assert!(poly.0.iter().all(|x| x.0 < BaseField::Q));
    }

    #[test]
    fn rej_bounded_poly_outputs_expected_ranges() {
        let rho = [13u8; 32];

        let p2 = rej_bounded_poly(&rho, Eta::Two, 0);
        assert!(all_abs_le_eta(&p2, 2));

        let p4 = rej_bounded_poly(&rho, Eta::Four, 0);
        assert!(all_abs_le_eta(&p4, 4));
    }

    #[test]
    fn expand_a_is_deterministic() {
        let rho = [17u8; 32];
        let a1 = expand_a::<U3, U4>(&rho);
        let a2 = expand_a::<U3, U4>(&rho);
        assert_eq!(a1, a2);
    }

    #[test]
    fn expand_s_is_deterministic_and_bounded() {
        let rho = [19u8; 32];
        let s1 = expand_s::<U4>(&rho, Eta::Two, 0);
        let s2 = expand_s::<U4>(&rho, Eta::Two, 0);

        assert_eq!(s1, s2);
        assert!(s1.0.iter().all(|p| all_abs_le_eta(p, 2)));
    }

    #[test]
    fn expand_mask_is_deterministic_for_both_gamma1_sizes() {
        let rho = [23u8; 64];

        let m17a = expand_mask::<U3, Gamma1Power17>(&rho, 5);
        let m17b = expand_mask::<U3, Gamma1Power17>(&rho, 5);
        assert_eq!(m17a, m17b);
        assert!(m17a.0.iter().all(|p| all_mask_in_range(p, 1 << 17)));

        let m19a = expand_mask::<U3, Gamma1Power19>(&rho, 7);
        let m19b = expand_mask::<U3, Gamma1Power19>(&rho, 7);
        assert_eq!(m19a, m19b);
        assert!(m19a.0.iter().all(|p| all_mask_in_range(p, 1 << 19)));
    }
}
