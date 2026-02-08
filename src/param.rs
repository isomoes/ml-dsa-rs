//! Parameter sets for ML-DSA
//!
//! This module encapsulates all compile-time logic related to parameter-set dependent sizes.
//! `ParameterSet` captures the parameters described by the FIPS 204 specification.
//! `SamplingSize` and `MaskSamplingSize` are upstream traits that provide basic logic about
//! sampling sizes. `SigningKeyParams`, `VerifyingKeyParams`, and `SignatureParams` are
//! blanket-implemented for all `ParameterSet` types, providing encoding/decoding and
//! concatenation/splitting logic for keys and signatures.

use crate::{
    algebra::{Polynomial, Vector},
    encode::{BitPack, RangeEncodedPolynomialSize, RangeEncodedVectorSize, RangeEncodingSize},
    module_lattice::encode::{
        ArraySize, Encode, EncodedPolynomialSize, EncodedVectorSize, EncodingSize,
    },
    util::{B32, B64},
};
use core::{
    fmt::Debug,
    ops::{Add, Div, Mul, Rem, Sub},
};
use hybrid_array::{
    typenum::{
        Diff, Len, Length, Prod, Shleft, Sum, Unsigned, U0, U1, U128, U13, U2, U23, U32, U320, U4,
        U416, U64,
    },
    Array,
};

// ============================================================================
// Compile-time constants from the FIPS 204 specification
// ============================================================================

/// Q = 2^23 - 2^13 + 1 = 8380417 as a type-level integer
pub(crate) type SpecQ = Sum<Diff<Shleft<U1, U23>, Shleft<U1, U13>>, U1>;

/// d = 13 (the number of dropped bits in Power2Round)
pub(crate) type SpecD = U13;

/// Q - 1 as a type-level integer
pub(crate) type QMinus1 = Diff<SpecQ, U1>;

/// bitlen(q) - d = 10 (bit width for t1 encoding)
pub(crate) type BitlenQMinusD = Diff<Length<SpecQ>, SpecD>;

/// 2^(d-1) as a type-level integer
pub(crate) type Pow2DMinus1 = Shleft<U1, Diff<SpecD, U1>>;

/// 2^(d-1) - 1 as a type-level integer
pub(crate) type Pow2DMinus1Minus1 = Diff<Pow2DMinus1, U1>;

// ============================================================================
// Eta enum and SamplingSize trait
// ============================================================================

/// Eta parameter used by bounded sampling (Algorithm 31 RejBoundedPoly).
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub(crate) enum Eta {
    /// ETA = 2 (used by ML-DSA-44 and ML-DSA-87).
    Two,
    /// ETA = 4 (used by ML-DSA-65).
    Four,
}

/// An integer that describes a bit length to be used in bounded sampling.
///
/// Only `U2` and `U4` implement this trait, corresponding to the two valid
/// values of the `eta` parameter in ML-DSA.
#[allow(private_interfaces)]
pub trait SamplingSize: ArraySize + Len {
    /// The runtime `Eta` variant for this sampling size.
    const ETA: Eta;
}

#[allow(private_interfaces)]
impl SamplingSize for U2 {
    const ETA: Eta = Eta::Two;
}

#[allow(private_interfaces)]
impl SamplingSize for U4 {
    const ETA: Eta = Eta::Four;
}

// ============================================================================
// MaskSamplingSize — mask sampling for ExpandMask (Algorithm 34)
// ============================================================================

/// An integer that describes a mask sampling size for ExpandMask.
///
/// This trait is blanket-implemented for any `Gamma1` value where the range
/// encoding `[-(Gamma1-1), Gamma1]` is valid.
pub trait MaskSamplingSize: Unsigned {
    /// The number of bytes in one encoded mask polynomial.
    type SampleSize: ArraySize;

    /// Unpack a byte array into a polynomial with coefficients in `[-(Gamma1-1), Gamma1]`.
    fn unpack(v: &Array<u8, Self::SampleSize>) -> Polynomial;
}

impl<G> MaskSamplingSize for G
where
    G: Unsigned + Sub<U1>,
    (Diff<G, U1>, G): RangeEncodingSize,
{
    type SampleSize = RangeEncodedPolynomialSize<Diff<G, U1>, G>;

    fn unpack(v: &Array<u8, Self::SampleSize>) -> Polynomial {
        BitPack::<Diff<G, U1>, G>::unpack(v)
    }
}

// ============================================================================
// ParameterSet — core trait capturing ML-DSA parameters
// ============================================================================

/// A `ParameterSet` captures the parameters that describe a particular instance of ML-DSA.
/// There are three variants (ML-DSA-44, ML-DSA-65, ML-DSA-87), corresponding to three
/// different security levels.
pub trait ParameterSet {
    /// Number of rows in the A matrix
    type K: ArraySize;

    /// Number of columns in the A matrix
    type L: ArraySize;

    /// Private key range parameter (eta)
    type Eta: SamplingSize;

    /// Error size bound for y (gamma1 = 2^17 or 2^19)
    type Gamma1: MaskSamplingSize;

    /// Low-order rounding range (gamma2)
    type Gamma2: Unsigned;

    /// 2 * gamma2 as a type-level integer
    type TwoGamma2: Unsigned;

    /// Encoding width of the w1 polynomial: bitlen((q-1)/(2*gamma2) - 1)
    type W1Bits: EncodingSize;

    /// Collision strength of c_tilde, in bytes (lambda/4 in the spec)
    type Lambda: ArraySize;

    /// Max number of true values in the hint
    type Omega: ArraySize;

    /// Number of nonzero values in the polynomial c
    const TAU: usize;

    /// Beta = Tau * Eta
    #[allow(clippy::as_conversions)]
    #[allow(clippy::cast_possible_truncation)]
    const BETA: u32 = (Self::TAU as u32) * Self::Eta::U32;
}

// ============================================================================
// SigningKeyParams — encoding/decoding for signing keys
// ============================================================================

/// Trait providing encoding/decoding operations for ML-DSA signing keys.
///
/// This is blanket-implemented for all types implementing `ParameterSet`.
pub trait SigningKeyParams: ParameterSet {
    /// Size of encoded s1 vector
    type S1Size: ArraySize;
    /// Size of encoded s2 vector
    type S2Size: ArraySize;
    /// Size of encoded t0 vector
    type T0Size: ArraySize;
    /// Total size of the expanded signing key
    type SigningKeySize: ArraySize;

    /// Encode s1 vector (Algorithm 17 BitPack with eta range)
    fn encode_s1(s1: &Vector<Self::L>) -> EncodedS1<Self>;
    /// Decode s1 vector
    fn decode_s1(enc: &EncodedS1<Self>) -> Vector<Self::L>;

    /// Encode s2 vector (Algorithm 17 BitPack with eta range)
    fn encode_s2(s2: &Vector<Self::K>) -> EncodedS2<Self>;
    /// Decode s2 vector
    fn decode_s2(enc: &EncodedS2<Self>) -> Vector<Self::K>;

    /// Encode t0 vector (Algorithm 17 BitPack with 2^(d-1) range)
    fn encode_t0(t0: &Vector<Self::K>) -> EncodedT0<Self>;
    /// Decode t0 vector
    fn decode_t0(enc: &EncodedT0<Self>) -> Vector<Self::K>;

    /// Concatenate signing key components into a single byte array (Algorithm 24 skEncode)
    fn concat_sk(
        rho: B32,
        k: B32,
        tr: B64,
        s1: EncodedS1<Self>,
        s2: EncodedS2<Self>,
        t0: EncodedT0<Self>,
    ) -> ExpandedSigningKey<Self>;

    /// Split an encoded signing key into its components (Algorithm 25 skDecode)
    fn split_sk(
        enc: &ExpandedSigningKey<Self>,
    ) -> (
        &B32,
        &B32,
        &B64,
        &EncodedS1<Self>,
        &EncodedS2<Self>,
        &EncodedT0<Self>,
    );
}

/// Encoded s1 vector
pub(crate) type EncodedS1<P> = Array<u8, <P as SigningKeyParams>::S1Size>;
/// Encoded s2 vector
pub(crate) type EncodedS2<P> = Array<u8, <P as SigningKeyParams>::S2Size>;
/// Encoded t0 vector
pub(crate) type EncodedT0<P> = Array<u8, <P as SigningKeyParams>::T0Size>;

/// Size of the expanded signing key
pub(crate) type SigningKeySize<P> = <P as SigningKeyParams>::SigningKeySize;

/// A signing key encoded as a byte array
pub type ExpandedSigningKey<P> = Array<u8, SigningKeySize<P>>;

impl<P> SigningKeyParams for P
where
    P: ParameterSet,
    // Eta encoding: range [-Eta, Eta]
    P::Eta: Add<P::Eta>,
    Sum<P::Eta, P::Eta>: Len,
    Length<Sum<P::Eta, P::Eta>>: EncodingSize,
    // S1 encoding (L polynomials)
    EncodedPolynomialSize<Length<Sum<P::Eta, P::Eta>>>: Mul<P::L>,
    Prod<EncodedPolynomialSize<Length<Sum<P::Eta, P::Eta>>>, P::L>: ArraySize
        + Div<P::L, Output = EncodedPolynomialSize<Length<Sum<P::Eta, P::Eta>>>>
        + Rem<P::L, Output = U0>,
    // S2 encoding (K polynomials)
    EncodedPolynomialSize<Length<Sum<P::Eta, P::Eta>>>: Mul<P::K>,
    Prod<EncodedPolynomialSize<Length<Sum<P::Eta, P::Eta>>>, P::K>: ArraySize
        + Div<P::K, Output = EncodedPolynomialSize<Length<Sum<P::Eta, P::Eta>>>>
        + Rem<P::K, Output = U0>,
    // T0 encoding: 13-bit range, 416 bytes per polynomial
    U416: Mul<P::K>,
    Prod<U416, P::K>: ArraySize + Div<P::K, Output = U416> + Rem<P::K, Output = U0>,
    // Signing key concatenation: rho(32) + K(32) + tr(64) + s1 + s2 + t0
    U128: Add<Prod<EncodedPolynomialSize<Length<Sum<P::Eta, P::Eta>>>, P::L>>,
    Sum<U128, Prod<EncodedPolynomialSize<Length<Sum<P::Eta, P::Eta>>>, P::L>>: ArraySize
        + Add<Prod<EncodedPolynomialSize<Length<Sum<P::Eta, P::Eta>>>, P::K>>
        + Sub<U128, Output = Prod<EncodedPolynomialSize<Length<Sum<P::Eta, P::Eta>>>, P::L>>,
    Sum<
        Sum<U128, Prod<EncodedPolynomialSize<Length<Sum<P::Eta, P::Eta>>>, P::L>>,
        Prod<EncodedPolynomialSize<Length<Sum<P::Eta, P::Eta>>>, P::K>,
    >: ArraySize
        + Add<Prod<U416, P::K>>
        + Sub<
            Sum<U128, Prod<EncodedPolynomialSize<Length<Sum<P::Eta, P::Eta>>>, P::L>>,
            Output = Prod<EncodedPolynomialSize<Length<Sum<P::Eta, P::Eta>>>, P::K>,
        >,
    Sum<
        Sum<
            Sum<U128, Prod<EncodedPolynomialSize<Length<Sum<P::Eta, P::Eta>>>, P::L>>,
            Prod<EncodedPolynomialSize<Length<Sum<P::Eta, P::Eta>>>, P::K>,
        >,
        Prod<U416, P::K>,
    >: ArraySize
        + Sub<
            Sum<
                Sum<U128, Prod<EncodedPolynomialSize<Length<Sum<P::Eta, P::Eta>>>, P::L>>,
                Prod<EncodedPolynomialSize<Length<Sum<P::Eta, P::Eta>>>, P::K>,
            >,
            Output = Prod<U416, P::K>,
        >,
{
    type S1Size = RangeEncodedVectorSize<P::Eta, P::Eta, P::L>;
    type S2Size = RangeEncodedVectorSize<P::Eta, P::Eta, P::K>;
    type T0Size = RangeEncodedVectorSize<Pow2DMinus1Minus1, Pow2DMinus1, P::K>;
    type SigningKeySize = Sum<
        Sum<
            Sum<U128, RangeEncodedVectorSize<P::Eta, P::Eta, P::L>>,
            RangeEncodedVectorSize<P::Eta, P::Eta, P::K>,
        >,
        RangeEncodedVectorSize<Pow2DMinus1Minus1, Pow2DMinus1, P::K>,
    >;

    fn encode_s1(s1: &Vector<Self::L>) -> EncodedS1<Self> {
        BitPack::<P::Eta, P::Eta>::pack(s1)
    }

    fn decode_s1(enc: &EncodedS1<Self>) -> Vector<Self::L> {
        BitPack::<P::Eta, P::Eta>::unpack(enc)
    }

    fn encode_s2(s2: &Vector<Self::K>) -> EncodedS2<Self> {
        BitPack::<P::Eta, P::Eta>::pack(s2)
    }

    fn decode_s2(enc: &EncodedS2<Self>) -> Vector<Self::K> {
        BitPack::<P::Eta, P::Eta>::unpack(enc)
    }

    fn encode_t0(t0: &Vector<Self::K>) -> EncodedT0<Self> {
        BitPack::<Pow2DMinus1Minus1, Pow2DMinus1>::pack(t0)
    }

    fn decode_t0(enc: &EncodedT0<Self>) -> Vector<Self::K> {
        BitPack::<Pow2DMinus1Minus1, Pow2DMinus1>::unpack(enc)
    }

    #[allow(non_snake_case)]
    fn concat_sk(
        rho: B32,
        k: B32,
        tr: B64,
        s1: EncodedS1<Self>,
        s2: EncodedS2<Self>,
        t0: EncodedT0<Self>,
    ) -> ExpandedSigningKey<Self> {
        rho.concat(k).concat(tr).concat(s1).concat(s2).concat(t0)
    }

    #[allow(non_snake_case)]
    fn split_sk(
        enc: &ExpandedSigningKey<Self>,
    ) -> (
        &B32,
        &B32,
        &B64,
        &EncodedS1<Self>,
        &EncodedS2<Self>,
        &EncodedT0<Self>,
    ) {
        let (enc, t0) = enc.split_ref();
        let (enc, s2) = enc.split_ref();
        let (enc, s1) = enc.split_ref();
        let (enc, tr) = enc.split_ref::<U64>();
        let (rho, k) = enc.split_ref();
        (rho, k, tr, s1, s2, t0)
    }
}

// ============================================================================
// VerifyingKeyParams — encoding/decoding for verifying keys
// ============================================================================

/// Trait providing encoding/decoding operations for ML-DSA verifying keys.
///
/// This is blanket-implemented for all types implementing `ParameterSet`.
pub trait VerifyingKeyParams: ParameterSet {
    /// Size of encoded t1 vector
    type T1Size: ArraySize;
    /// Total size of the encoded verifying key
    type VerifyingKeySize: ArraySize;

    /// Encode t1 vector (Algorithm 16 SimpleBitPack with bitlen(q)-d bits)
    fn encode_t1(t1: &Vector<Self::K>) -> EncodedT1<Self>;
    /// Decode t1 vector
    fn decode_t1(enc: &EncodedT1<Self>) -> Vector<Self::K>;

    /// Concatenate verifying key components (Algorithm 22 pkEncode)
    fn concat_vk(rho: B32, t1: EncodedT1<Self>) -> EncodedVerifyingKey<Self>;
    /// Split an encoded verifying key into its components (Algorithm 23 pkDecode)
    fn split_vk(enc: &EncodedVerifyingKey<Self>) -> (&B32, &EncodedT1<Self>);
}

/// Size of the encoded verifying key
pub(crate) type VerifyingKeySize<P> = <P as VerifyingKeyParams>::VerifyingKeySize;

/// Encoded t1 vector
pub(crate) type EncodedT1<P> = Array<u8, <P as VerifyingKeyParams>::T1Size>;

/// A verifying key encoded as a byte array
pub type EncodedVerifyingKey<P> = Array<u8, VerifyingKeySize<P>>;

impl<P> VerifyingKeyParams for P
where
    P: ParameterSet,
    // T1 encoding: 10-bit encoding (bitlen(q)-d = 10), 320 bytes per polynomial
    U320: Mul<P::K>,
    Prod<U320, P::K>: ArraySize + Div<P::K, Output = U320> + Rem<P::K, Output = U0>,
    // Verifying key: rho(32) + t1
    U32: Add<Prod<U320, P::K>>,
    Sum<U32, U32>: ArraySize,
    Sum<U32, Prod<U320, P::K>>: ArraySize + Sub<U32, Output = Prod<U320, P::K>>,
{
    type T1Size = EncodedVectorSize<BitlenQMinusD, P::K>;
    type VerifyingKeySize = Sum<U32, Self::T1Size>;

    fn encode_t1(t1: &Vector<P::K>) -> EncodedT1<Self> {
        Encode::<BitlenQMinusD>::encode(t1)
    }

    fn decode_t1(enc: &EncodedT1<Self>) -> Vector<Self::K> {
        Encode::<BitlenQMinusD>::decode(enc)
    }

    fn concat_vk(rho: B32, t1: EncodedT1<Self>) -> EncodedVerifyingKey<Self> {
        rho.concat(t1)
    }

    fn split_vk(enc: &EncodedVerifyingKey<Self>) -> (&B32, &EncodedT1<Self>) {
        enc.split_ref()
    }
}

// ============================================================================
// SignatureParams — encoding/decoding for signatures
// ============================================================================

/// Trait providing encoding/decoding operations for ML-DSA signatures.
///
/// This is blanket-implemented for all types implementing `ParameterSet`.
pub trait SignatureParams: ParameterSet {
    /// Size of encoded w1 vector
    type W1Size: ArraySize;
    /// Size of encoded z vector
    type ZSize: ArraySize;
    /// Size of encoded hint
    type HintSize: ArraySize;
    /// Total size of the encoded signature
    type SignatureSize: ArraySize;

    /// Gamma1 - Beta (used for infinity norm check on z)
    const GAMMA1_MINUS_BETA: u32;
    /// Gamma2 - Beta (used for infinity norm check on r0)
    const GAMMA2_MINUS_BETA: u32;

    /// Split hint bytes into indices and cuts
    fn split_hint(y: &EncodedHint<Self>) -> (&EncodedHintIndices<Self>, &EncodedHintCuts<Self>);

    /// Encode w1 vector (Algorithm 16 SimpleBitPack)
    fn encode_w1(w1: &Vector<Self::K>) -> EncodedW1<Self>;
    /// Decode w1 vector
    fn decode_w1(enc: &EncodedW1<Self>) -> Vector<Self::K>;

    /// Encode z vector (Algorithm 17 BitPack with gamma1 range)
    fn encode_z(z: &Vector<Self::L>) -> EncodedZ<Self>;
    /// Decode z vector
    fn decode_z(enc: &EncodedZ<Self>) -> Vector<Self::L>;

    /// Concatenate signature components (Algorithm 26 sigEncode)
    fn concat_sig(
        c_tilde: EncodedCTilde<Self>,
        z: EncodedZ<Self>,
        h: EncodedHint<Self>,
    ) -> EncodedSignature<Self>;

    /// Split an encoded signature into its components (Algorithm 27 sigDecode)
    fn split_sig(
        enc: &EncodedSignature<Self>,
    ) -> (&EncodedCTilde<Self>, &EncodedZ<Self>, &EncodedHint<Self>);
}

/// Size of the encoded signature
pub(crate) type SignatureSize<P> = <P as SignatureParams>::SignatureSize;

/// Encoded c_tilde (lambda bytes)
pub(crate) type EncodedCTilde<P> = Array<u8, <P as ParameterSet>::Lambda>;
/// Encoded w1 vector
pub(crate) type EncodedW1<P> = Array<u8, <P as SignatureParams>::W1Size>;
/// Encoded z vector
pub(crate) type EncodedZ<P> = Array<u8, <P as SignatureParams>::ZSize>;
/// Encoded hint indices (omega bytes)
pub(crate) type EncodedHintIndices<P> = Array<u8, <P as ParameterSet>::Omega>;
/// Encoded hint cuts (K bytes)
pub(crate) type EncodedHintCuts<P> = Array<u8, <P as ParameterSet>::K>;
/// Encoded hint (omega + K bytes)
pub(crate) type EncodedHint<P> = Array<u8, <P as SignatureParams>::HintSize>;

/// A signature encoded as a byte array
pub type EncodedSignature<P> = Array<u8, SignatureSize<P>>;

impl<P> SignatureParams for P
where
    P: ParameterSet,
    // W1 encoding
    U32: Mul<P::W1Bits>,
    EncodedPolynomialSize<P::W1Bits>: Mul<P::K>,
    Prod<EncodedPolynomialSize<P::W1Bits>, P::K>:
        ArraySize + Div<P::K, Output = EncodedPolynomialSize<P::W1Bits>> + Rem<P::K, Output = U0>,
    // Z encoding: range [-(Gamma1-1), Gamma1]
    P::Gamma1: Sub<U1>,
    (Diff<P::Gamma1, U1>, P::Gamma1): RangeEncodingSize,
    RangeEncodedPolynomialSize<Diff<P::Gamma1, U1>, P::Gamma1>: Mul<P::L>,
    Prod<RangeEncodedPolynomialSize<Diff<P::Gamma1, U1>, P::Gamma1>, P::L>: ArraySize
        + Div<P::L, Output = RangeEncodedPolynomialSize<Diff<P::Gamma1, U1>, P::Gamma1>>
        + Rem<P::L, Output = U0>,
    // Hint: omega + K bytes
    P::Omega: Add<P::K>,
    Sum<P::Omega, P::K>: ArraySize + Sub<P::Omega, Output = P::K>,
    // Signature: lambda + z_size + hint_size
    P::Lambda: Add<Prod<RangeEncodedPolynomialSize<Diff<P::Gamma1, U1>, P::Gamma1>, P::L>>,
    Sum<P::Lambda, Prod<RangeEncodedPolynomialSize<Diff<P::Gamma1, U1>, P::Gamma1>, P::L>>:
        ArraySize
            + Add<Sum<P::Omega, P::K>>
            + Sub<
                P::Lambda,
                Output = Prod<RangeEncodedPolynomialSize<Diff<P::Gamma1, U1>, P::Gamma1>, P::L>,
            >,
    Sum<
        Sum<P::Lambda, Prod<RangeEncodedPolynomialSize<Diff<P::Gamma1, U1>, P::Gamma1>, P::L>>,
        Sum<P::Omega, P::K>,
    >: ArraySize
        + Sub<
            Sum<P::Lambda, Prod<RangeEncodedPolynomialSize<Diff<P::Gamma1, U1>, P::Gamma1>, P::L>>,
            Output = Sum<P::Omega, P::K>,
        >,
{
    type W1Size = EncodedVectorSize<Self::W1Bits, P::K>;
    type ZSize = RangeEncodedVectorSize<Diff<P::Gamma1, U1>, P::Gamma1, P::L>;
    type HintSize = Sum<P::Omega, P::K>;
    type SignatureSize = Sum<Sum<P::Lambda, Self::ZSize>, Self::HintSize>;

    const GAMMA1_MINUS_BETA: u32 = P::Gamma1::U32 - P::BETA;
    const GAMMA2_MINUS_BETA: u32 = P::Gamma2::U32 - P::BETA;

    fn split_hint(y: &EncodedHint<Self>) -> (&EncodedHintIndices<Self>, &EncodedHintCuts<Self>) {
        y.split_ref()
    }

    fn encode_w1(w1: &Vector<Self::K>) -> EncodedW1<Self> {
        Encode::<Self::W1Bits>::encode(w1)
    }

    fn decode_w1(enc: &EncodedW1<Self>) -> Vector<Self::K> {
        Encode::<Self::W1Bits>::decode(enc)
    }

    fn encode_z(z: &Vector<Self::L>) -> EncodedZ<Self> {
        BitPack::<Diff<P::Gamma1, U1>, P::Gamma1>::pack(z)
    }

    fn decode_z(enc: &EncodedZ<Self>) -> Vector<Self::L> {
        BitPack::<Diff<P::Gamma1, U1>, P::Gamma1>::unpack(enc)
    }

    fn concat_sig(
        c_tilde: EncodedCTilde<P>,
        z: EncodedZ<P>,
        h: EncodedHint<P>,
    ) -> EncodedSignature<P> {
        c_tilde.concat(z).concat(h)
    }

    fn split_sig(enc: &EncodedSignature<P>) -> (&EncodedCTilde<P>, &EncodedZ<P>, &EncodedHint<P>) {
        let (enc, h) = enc.split_ref();
        let (c_tilde, z) = enc.split_ref();
        (c_tilde, z, h)
    }
}

// ============================================================================
// MlDsaParams — super-trait combining all parameter traits
// ============================================================================

/// An instance of `MlDsaParams` defines all of the parameters necessary for ML-DSA operations.
/// Typically this is done by implementing `ParameterSet` with values that will fit into the
/// blanket implementations of `SigningKeyParams`, `VerifyingKeyParams`, and `SignatureParams`.
pub trait MlDsaParams:
    SigningKeyParams + VerifyingKeyParams + SignatureParams + Debug + Default + PartialEq + Clone
{
}

impl<T> MlDsaParams for T where
    T: SigningKeyParams
        + VerifyingKeyParams
        + SignatureParams
        + Debug
        + Default
        + PartialEq
        + Clone
{
}

// ============================================================================
// Concrete parameter set implementations
// ============================================================================

/// ML-DSA parameter set 44 (Security Category 2)
#[derive(Default, Clone, Debug, PartialEq)]
pub struct MlDsa44;

impl ParameterSet for MlDsa44 {
    type K = hybrid_array::typenum::U4;
    type L = hybrid_array::typenum::U4;
    type Eta = U2;
    type Gamma1 = Shleft<U1, hybrid_array::typenum::U17>;
    type Gamma2 = hybrid_array::typenum::Quot<QMinus1, hybrid_array::typenum::U88>;
    type TwoGamma2 = Prod<U2, Self::Gamma2>;
    type W1Bits = Length<Diff<hybrid_array::typenum::Quot<hybrid_array::typenum::U88, U2>, U1>>;
    type Lambda = U32;
    type Omega = hybrid_array::typenum::U80;
    const TAU: usize = 39;
}

/// ML-DSA parameter set 65 (Security Category 3)
#[derive(Default, Clone, Debug, PartialEq)]
pub struct MlDsa65;

impl ParameterSet for MlDsa65 {
    type K = hybrid_array::typenum::U6;
    type L = hybrid_array::typenum::U5;
    type Eta = U4;
    type Gamma1 = Shleft<U1, hybrid_array::typenum::U19>;
    type Gamma2 = hybrid_array::typenum::Quot<QMinus1, U32>;
    type TwoGamma2 = Prod<U2, Self::Gamma2>;
    type W1Bits = Length<Diff<hybrid_array::typenum::Quot<U32, U2>, U1>>;
    type Lambda = hybrid_array::typenum::U48;
    type Omega = hybrid_array::typenum::U55;
    const TAU: usize = 49;
}

/// ML-DSA parameter set 87 (Security Category 5)
#[derive(Default, Clone, Debug, PartialEq)]
pub struct MlDsa87;

impl ParameterSet for MlDsa87 {
    type K = hybrid_array::typenum::U8;
    type L = hybrid_array::typenum::U7;
    type Eta = U2;
    type Gamma1 = Shleft<U1, hybrid_array::typenum::U19>;
    type Gamma2 = hybrid_array::typenum::Quot<QMinus1, U32>;
    type TwoGamma2 = Prod<U2, Self::Gamma2>;
    type W1Bits = Length<Diff<hybrid_array::typenum::Quot<U32, U2>, U1>>;
    type Lambda = U64;
    type Omega = hybrid_array::typenum::U75;
    const TAU: usize = 60;
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use hybrid_array::typenum::Unsigned;

    // Verify FIPS 204 Table 1 parameter values
    #[test]
    fn ml_dsa_44_parameters() {
        assert_eq!(<MlDsa44 as ParameterSet>::K::USIZE, 4);
        assert_eq!(<MlDsa44 as ParameterSet>::L::USIZE, 4);
        assert_eq!(<MlDsa44 as ParameterSet>::Eta::USIZE, 2);
        assert_eq!(<MlDsa44 as ParameterSet>::Gamma1::USIZE, 1 << 17);
        assert_eq!(
            <MlDsa44 as ParameterSet>::Gamma2::USIZE,
            (8_380_417 - 1) / 88
        );
        assert_eq!(
            <MlDsa44 as ParameterSet>::TwoGamma2::USIZE,
            2 * (8_380_417 - 1) / 88
        );
        assert_eq!(<MlDsa44 as ParameterSet>::Lambda::USIZE, 32);
        assert_eq!(<MlDsa44 as ParameterSet>::Omega::USIZE, 80);
        assert_eq!(MlDsa44::TAU, 39);
        assert_eq!(MlDsa44::BETA, 39 * 2);
    }

    #[test]
    fn ml_dsa_65_parameters() {
        assert_eq!(<MlDsa65 as ParameterSet>::K::USIZE, 6);
        assert_eq!(<MlDsa65 as ParameterSet>::L::USIZE, 5);
        assert_eq!(<MlDsa65 as ParameterSet>::Eta::USIZE, 4);
        assert_eq!(<MlDsa65 as ParameterSet>::Gamma1::USIZE, 1 << 19);
        assert_eq!(
            <MlDsa65 as ParameterSet>::Gamma2::USIZE,
            (8_380_417 - 1) / 32
        );
        assert_eq!(
            <MlDsa65 as ParameterSet>::TwoGamma2::USIZE,
            2 * (8_380_417 - 1) / 32
        );
        assert_eq!(<MlDsa65 as ParameterSet>::Lambda::USIZE, 48);
        assert_eq!(<MlDsa65 as ParameterSet>::Omega::USIZE, 55);
        assert_eq!(MlDsa65::TAU, 49);
        assert_eq!(MlDsa65::BETA, 49 * 4);
    }

    #[test]
    fn ml_dsa_87_parameters() {
        assert_eq!(<MlDsa87 as ParameterSet>::K::USIZE, 8);
        assert_eq!(<MlDsa87 as ParameterSet>::L::USIZE, 7);
        assert_eq!(<MlDsa87 as ParameterSet>::Eta::USIZE, 2);
        assert_eq!(<MlDsa87 as ParameterSet>::Gamma1::USIZE, 1 << 19);
        assert_eq!(
            <MlDsa87 as ParameterSet>::Gamma2::USIZE,
            (8_380_417 - 1) / 32
        );
        assert_eq!(
            <MlDsa87 as ParameterSet>::TwoGamma2::USIZE,
            2 * (8_380_417 - 1) / 32
        );
        assert_eq!(<MlDsa87 as ParameterSet>::Lambda::USIZE, 64);
        assert_eq!(<MlDsa87 as ParameterSet>::Omega::USIZE, 75);
        assert_eq!(MlDsa87::TAU, 60);
        assert_eq!(MlDsa87::BETA, 60 * 2);
    }

    // Verify FIPS 204 Table 2 key/signature sizes
    //           sk    pk    sig
    // ML-DSA-44 2560  1312  2420
    // ML-DSA-65 4032  1952  3309
    // ML-DSA-87 4896  2592  4627
    #[test]
    fn output_sizes() {
        assert_eq!(SigningKeySize::<MlDsa44>::USIZE, 2560);
        assert_eq!(VerifyingKeySize::<MlDsa44>::USIZE, 1312);
        assert_eq!(SignatureSize::<MlDsa44>::USIZE, 2420);

        assert_eq!(SigningKeySize::<MlDsa65>::USIZE, 4032);
        assert_eq!(VerifyingKeySize::<MlDsa65>::USIZE, 1952);
        assert_eq!(SignatureSize::<MlDsa65>::USIZE, 3309);

        assert_eq!(SigningKeySize::<MlDsa87>::USIZE, 4896);
        assert_eq!(VerifyingKeySize::<MlDsa87>::USIZE, 2592);
        assert_eq!(SignatureSize::<MlDsa87>::USIZE, 4627);
    }

    // Verify derived constants
    #[test]
    fn derived_constants() {
        // GAMMA1_MINUS_BETA
        assert_eq!(
            <MlDsa44 as SignatureParams>::GAMMA1_MINUS_BETA,
            (1u32 << 17) - 39 * 2
        );
        assert_eq!(
            <MlDsa65 as SignatureParams>::GAMMA1_MINUS_BETA,
            (1u32 << 19) - 49 * 4
        );
        assert_eq!(
            <MlDsa87 as SignatureParams>::GAMMA1_MINUS_BETA,
            (1u32 << 19) - 60 * 2
        );

        // GAMMA2_MINUS_BETA
        assert_eq!(
            <MlDsa44 as SignatureParams>::GAMMA2_MINUS_BETA,
            95232 - 39 * 2
        );
        assert_eq!(
            <MlDsa65 as SignatureParams>::GAMMA2_MINUS_BETA,
            261888 - 49 * 4
        );
        assert_eq!(
            <MlDsa87 as SignatureParams>::GAMMA2_MINUS_BETA,
            261888 - 60 * 2
        );
    }

    // Verify compile-time constants
    #[test]
    fn spec_constants() {
        assert_eq!(SpecQ::USIZE, 8_380_417);
        assert_eq!(QMinus1::USIZE, 8_380_416);
        assert_eq!(SpecD::USIZE, 13);
        assert_eq!(BitlenQMinusD::USIZE, 10);
        assert_eq!(Pow2DMinus1::USIZE, 4096);
        assert_eq!(Pow2DMinus1Minus1::USIZE, 4095);
    }

    // Verify that MlDsaParams is satisfied for all three parameter sets
    #[test]
    fn ml_dsa_params_trait_satisfied() {
        fn assert_ml_dsa_params<P: MlDsaParams>() {}
        assert_ml_dsa_params::<MlDsa44>();
        assert_ml_dsa_params::<MlDsa65>();
        assert_ml_dsa_params::<MlDsa87>();
    }

    // Verify W1Bits values
    #[test]
    fn w1_bits() {
        // ML-DSA-44: (q-1)/(2*gamma2) - 1 = 88/2 - 1 = 43, bitlen(43) = 6
        // But actually: W1Bits = Length<43> which is ceil(log2(44)) = 6
        assert_eq!(<MlDsa44 as ParameterSet>::W1Bits::USIZE, 6);

        // ML-DSA-65/87: (q-1)/(2*gamma2) - 1 = 32/2 - 1 = 15, bitlen(15) = 4
        assert_eq!(<MlDsa65 as ParameterSet>::W1Bits::USIZE, 4);
        assert_eq!(<MlDsa87 as ParameterSet>::W1Bits::USIZE, 4);
    }
}
