//! Pure Rust implementation of ML-DSA (formerly known as CRYSTALS-Dilithium) as
//! described in FIPS-204 (final)

#![warn(missing_docs, rust_2018_idioms)]
#![allow(non_snake_case)]
// Note: unsafe code is allowed only in module_lattice for performance-critical operations

// Foundation module that needs unsafe for array operations
pub mod module_lattice;

// All other modules forbid unsafe code
#[forbid(unsafe_code)]
pub mod algebra;
#[forbid(unsafe_code)]
pub mod crypto;
#[forbid(unsafe_code)]
pub mod encode;
#[forbid(unsafe_code)]
pub mod hint;
#[forbid(unsafe_code)]
pub mod ntt;
#[forbid(unsafe_code)]
pub mod param;
#[forbid(unsafe_code)]
pub mod sampling;
#[forbid(unsafe_code)]
pub mod util;

pub use crate::param::{
    EncodedSignature, EncodedVerifyingKey, ExpandedSigningKey, MlDsa44, MlDsa65, MlDsa87,
    MlDsaParams,
};
pub use signature::{self, Error};

use crate::algebra::{AlgebraExt, Elem, NttMatrix, NttVector, Vector};
use crate::crypto::H;
use crate::hint::Hint;
use crate::ntt::{Ntt, NttInverse};
use crate::param::{ParameterSet, SamplingSize, SpecQ};
use crate::sampling::{expand_a, expand_mask, expand_s, sample_in_ball};
use crate::util::{B32, B64};
use core::fmt;
use hybrid_array::{typenum::Unsigned, Array};
use signature::{MultipartSigner, MultipartVerifier, Signer};

#[cfg(feature = "rand_core")]
use rand_core::CryptoRng;

#[cfg(feature = "rand_core")]
use signature::{RandomizedMultipartSigner, RandomizedSigner};

/// ML-DSA seeds are signing (private) keys, which are consistently 32-bytes across all security
/// levels, and are the preferred serialization for representing such keys.
pub type Seed = B32;

// ============================================================================
// Signature
// ============================================================================

/// An ML-DSA signature
#[derive(Clone, PartialEq, Debug)]
pub struct Signature<P: MlDsaParams> {
    c_tilde: Array<u8, P::Lambda>,
    z: Vector<P::L>,
    h: Hint<P>,
}

impl<P: MlDsaParams> Signature<P> {
    /// Encode the signature in a fixed-size byte array.
    // Algorithm 26 sigEncode
    pub fn encode(&self) -> EncodedSignature<P> {
        let c_tilde = self.c_tilde.clone();
        let z = P::encode_z(&self.z);
        let h = self.h.bit_pack();
        P::concat_sig(c_tilde, z, h)
    }

    /// Decode the signature from an appropriately sized byte array.
    // Algorithm 27 sigDecode
    pub fn decode(enc: &EncodedSignature<P>) -> Option<Self> {
        let (c_tilde, z, h) = P::split_sig(enc);

        let c_tilde = c_tilde.clone();
        let z = P::decode_z(z);
        let h = Hint::bit_unpack(h)?;

        if z.infinity_norm() >= P::GAMMA1_MINUS_BETA {
            return None;
        }

        Some(Self { c_tilde, z, h })
    }
}

impl<'a, P: MlDsaParams> TryFrom<&'a [u8]> for Signature<P> {
    type Error = Error;

    fn try_from(value: &'a [u8]) -> Result<Self, Self::Error> {
        let enc = EncodedSignature::<P>::try_from(value).map_err(|_| Error::new())?;
        Self::decode(&enc).ok_or(Error::new())
    }
}

impl<P: MlDsaParams> TryFrom<Signature<P>> for EncodedSignature<P> {
    type Error = Error;

    fn try_from(sig: Signature<P>) -> Result<Self, Self::Error> {
        Ok(sig.encode())
    }
}

impl<P: MlDsaParams> signature::SignatureEncoding for Signature<P> {
    type Repr = EncodedSignature<P>;
}

// ============================================================================
// MuBuilder — domain separation helper
// ============================================================================

struct MuBuilder(H);

impl MuBuilder {
    fn new(tr: &[u8], ctx: &[u8]) -> Self {
        let mut h = H::default();
        h.absorb(tr);
        h.absorb(&[0]);
        #[allow(clippy::cast_possible_truncation, clippy::as_conversions)]
        h.absorb(&[ctx.len() as u8]);
        h.absorb(ctx);

        Self(h)
    }

    fn internal(tr: &[u8], Mp: &[&[u8]]) -> B64 {
        let mut h = H::default();
        h.absorb(tr);

        for m in Mp {
            h.absorb(m);
        }

        h.squeeze_new_array()
    }

    fn message(mut self, M: &[&[u8]]) -> B64 {
        for m in M {
            self.0.absorb(m);
        }

        self.0.squeeze_new_array()
    }
}

// ============================================================================
// KeyPair
// ============================================================================

/// An ML-DSA key pair
pub struct KeyPair<P: MlDsaParams> {
    signing_key: SigningKey<P>,
    verifying_key: VerifyingKey<P>,
    seed: B32,
}

impl<P: MlDsaParams> KeyPair<P> {
    /// The signing key of the key pair
    pub fn signing_key(&self) -> &SigningKey<P> {
        &self.signing_key
    }

    /// The verifying key of the key pair
    pub fn verifying_key(&self) -> &VerifyingKey<P> {
        &self.verifying_key
    }

    /// Serialize the [`Seed`] value: 32-bytes which can be used to reconstruct the
    /// [`KeyPair`].
    #[inline]
    pub fn to_seed(&self) -> Seed {
        self.seed.clone()
    }
}

impl<P: MlDsaParams> AsRef<VerifyingKey<P>> for KeyPair<P> {
    fn as_ref(&self) -> &VerifyingKey<P> {
        &self.verifying_key
    }
}

impl<P: MlDsaParams> fmt::Debug for KeyPair<P> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("KeyPair")
            .field("verifying_key", &self.verifying_key)
            .finish_non_exhaustive()
    }
}

impl<P: MlDsaParams> signature::KeypairRef for KeyPair<P> {
    type VerifyingKey = VerifyingKey<P>;
}

/// The `Signer` implementation for `KeyPair` uses the optional deterministic variant of ML-DSA,
/// and only supports signing with an empty context string.
impl<P: MlDsaParams> Signer<Signature<P>> for KeyPair<P> {
    fn try_sign(&self, msg: &[u8]) -> Result<Signature<P>, Error> {
        self.try_multipart_sign(&[msg])
    }
}

/// The `MultipartSigner` implementation for `KeyPair` uses the optional deterministic variant of
/// ML-DSA, and only supports signing with an empty context string.
impl<P: MlDsaParams> MultipartSigner<Signature<P>> for KeyPair<P> {
    fn try_multipart_sign(&self, msg: &[&[u8]]) -> Result<Signature<P>, Error> {
        self.signing_key.raw_sign_deterministic(msg, &[])
    }
}

// ============================================================================
// SigningKey
// ============================================================================

/// An ML-DSA signing key
#[derive(Clone, PartialEq)]
pub struct SigningKey<P: MlDsaParams> {
    rho: B32,
    K: B32,
    tr: B64,
    s1: Vector<P::L>,
    s2: Vector<P::K>,
    t0: Vector<P::K>,

    // Derived values
    s1_hat: NttVector<P::L>,
    s2_hat: NttVector<P::K>,
    t0_hat: NttVector<P::K>,
    A_hat: NttMatrix<P::K, P::L>,
}

impl<P: MlDsaParams> fmt::Debug for SigningKey<P> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("SigningKey").finish_non_exhaustive()
    }
}

impl<P: MlDsaParams> SigningKey<P> {
    fn new(
        rho: B32,
        K: B32,
        tr: B64,
        s1: Vector<P::L>,
        s2: Vector<P::K>,
        t0: Vector<P::K>,
        A_hat: Option<NttMatrix<P::K, P::L>>,
    ) -> Self {
        let A_hat = A_hat.unwrap_or_else(|| expand_a(&rho));
        let s1_hat = s1.ntt();
        let s2_hat = s2.ntt();
        let t0_hat = t0.ntt();

        Self {
            rho,
            K,
            tr,
            s1,
            s2,
            t0,
            s1_hat,
            s2_hat,
            t0_hat,
            A_hat,
        }
    }

    /// Deterministically generate a signing key from the specified seed.
    ///
    /// This method reflects the ML-DSA.KeyGen_internal algorithm from FIPS 204, but only returns a
    /// signing key.
    #[must_use]
    pub fn from_seed(seed: &Seed) -> Self {
        let kp = P::from_seed(seed);
        kp.signing_key
    }

    /// This method reflects the ML-DSA.Sign_internal algorithm from FIPS 204. It does not
    /// include the domain separator that distinguishes between the normal and pre-hashed cases,
    /// and it does not separate the context string from the rest of the message.
    // Algorithm 7 ML-DSA.Sign_internal
    pub fn sign_internal(&self, Mp: &[&[u8]], rnd: &B32) -> Signature<P> {
        let mu = MuBuilder::internal(&self.tr, Mp);
        self.raw_sign_mu(&mu, rnd)
    }

    fn raw_sign_mu(&self, mu: &B64, rnd: &B32) -> Signature<P> {
        // Compute the private random seed
        let mut rhopp_ctx = H::default();
        rhopp_ctx.absorb(&self.K);
        rhopp_ctx.absorb(rnd);
        rhopp_ctx.absorb(mu);
        let rhopp: B64 = rhopp_ctx.squeeze_new_array();

        // Rejection sampling loop
        for kappa in (0..u16::MAX).step_by(P::L::USIZE) {
            let y = expand_mask::<P::L, P::Gamma1>(&rhopp, kappa);
            let w = (&self.A_hat * &y.ntt()).ntt_inverse();
            let w1 = w.high_bits::<P::TwoGamma2>();

            let w1_tilde = P::encode_w1(&w1);
            let mut c_ctx = H::default();
            c_ctx.absorb(mu);
            c_ctx.absorb(&w1_tilde);
            let c_tilde: Array<u8, P::Lambda> = c_ctx.squeeze_new_array();
            let c = sample_in_ball(&c_tilde, P::TAU);
            let c_hat = c.ntt();

            let cs1 = (&c_hat * &self.s1_hat).ntt_inverse();
            let cs2 = (&c_hat * &self.s2_hat).ntt_inverse();

            let z = &y + &cs1;
            let r0 = (&w - &cs2).low_bits::<P::TwoGamma2>();

            if z.infinity_norm() >= P::GAMMA1_MINUS_BETA
                || r0.infinity_norm() >= P::GAMMA2_MINUS_BETA
            {
                continue;
            }

            let ct0 = (&c_hat * &self.t0_hat).ntt_inverse();
            let minus_ct0 = -&ct0;
            let w_cs2_ct0 = &(&w - &cs2) + &ct0;
            let h = Hint::<P>::new(&minus_ct0, &w_cs2_ct0);

            if ct0.infinity_norm() >= P::Gamma2::U32 || h.hamming_weight() > P::Omega::USIZE {
                continue;
            }

            let z = z.mod_plus_minus::<SpecQ>();
            return Signature { c_tilde, z, h };
        }

        unreachable!("Rejection sampling failed to find a valid signature");
    }

    /// This method reflects the optional deterministic variant of the ML-DSA.Sign algorithm.
    ///
    /// # Errors
    ///
    /// This method will return an opaque error if the context string is more than 255 bytes long.
    // Algorithm 2 ML-DSA.Sign (optional deterministic variant)
    pub fn sign_deterministic(&self, M: &[u8], ctx: &[u8]) -> Result<Signature<P>, Error> {
        self.raw_sign_deterministic(&[M], ctx)
    }

    /// This method reflects the optional deterministic variant of the ML-DSA.Sign algorithm with a
    /// pre-computed mu.
    // Algorithm 2 ML-DSA.Sign (optional deterministic and pre-computed mu variant)
    pub fn sign_mu_deterministic(&self, mu: &B64) -> Signature<P> {
        let rnd = B32::default();
        self.raw_sign_mu(mu, &rnd)
    }

    fn raw_sign_deterministic(&self, Mp: &[&[u8]], ctx: &[u8]) -> Result<Signature<P>, Error> {
        if ctx.len() > 255 {
            return Err(Error::new());
        }

        let mu = MuBuilder::new(&self.tr, ctx).message(Mp);
        Ok(self.sign_mu_deterministic(&mu))
    }

    /// This method reflects the randomized ML-DSA.Sign algorithm.
    ///
    /// # Errors
    ///
    /// This method will return an opaque error if the context string is more than 255 bytes long,
    /// or if it fails to get enough randomness.
    // Algorithm 2 ML-DSA.Sign
    #[cfg(feature = "rand_core")]
    pub fn sign_randomized<R: rand_core::TryCryptoRng + ?Sized>(
        &self,
        M: &[u8],
        ctx: &[u8],
        rng: &mut R,
    ) -> Result<Signature<P>, Error> {
        self.raw_sign_randomized(&[M], ctx, rng)
    }

    #[cfg(feature = "rand_core")]
    fn raw_sign_randomized<R: rand_core::TryCryptoRng + ?Sized>(
        &self,
        Mp: &[&[u8]],
        ctx: &[u8],
        rng: &mut R,
    ) -> Result<Signature<P>, Error> {
        if ctx.len() > 255 {
            return Err(Error::new());
        }

        let mut rnd = B32::default();
        rng.try_fill_bytes(&mut rnd).map_err(|_| Error::new())?;

        let mu = MuBuilder::new(&self.tr, ctx).message(Mp);
        Ok(self.raw_sign_mu(&mu, &rnd))
    }

    /// Derive a `VerifyingKey` from this `SigningKey`.
    ///
    /// This is a convenience method that delegates to the [`signature::Keypair`] trait
    /// implementation.
    pub fn verifying_key(&self) -> VerifyingKey<P> {
        let kp: &dyn signature::Keypair<VerifyingKey = VerifyingKey<P>> = self;
        kp.verifying_key()
    }

    /// Encode the key in a fixed-size byte array.
    // Algorithm 24 skEncode
    pub fn to_expanded(&self) -> ExpandedSigningKey<P> {
        let s1_enc = P::encode_s1(&self.s1);
        let s2_enc = P::encode_s2(&self.s2);
        let t0_enc = P::encode_t0(&self.t0);
        P::concat_sk(
            self.rho.clone(),
            self.K.clone(),
            self.tr.clone(),
            s1_enc,
            s2_enc,
            t0_enc,
        )
    }

    /// Decode the key from an appropriately sized byte array.
    // Algorithm 25 skDecode
    pub fn from_expanded(enc: &ExpandedSigningKey<P>) -> Self {
        let (rho, K, tr, s1_enc, s2_enc, t0_enc) = P::split_sk(enc);
        Self::new(
            rho.clone(),
            K.clone(),
            tr.clone(),
            P::decode_s1(s1_enc),
            P::decode_s2(s2_enc),
            P::decode_t0(t0_enc),
            None,
        )
    }
}

/// The `Signer` implementation for `SigningKey` uses the optional deterministic variant of ML-DSA,
/// and only supports signing with an empty context string. If you would like to include a context
/// string, use the [`SigningKey::sign_deterministic`] method.
impl<P: MlDsaParams> Signer<Signature<P>> for SigningKey<P> {
    fn try_sign(&self, msg: &[u8]) -> Result<Signature<P>, Error> {
        self.try_multipart_sign(&[msg])
    }
}

/// The `MultipartSigner` implementation for `SigningKey` uses the optional deterministic variant
/// of ML-DSA, and only supports signing with an empty context string. If you would like to
/// include a context string, use the [`SigningKey::sign_deterministic`] method.
impl<P: MlDsaParams> MultipartSigner<Signature<P>> for SigningKey<P> {
    fn try_multipart_sign(&self, msg: &[&[u8]]) -> Result<Signature<P>, Error> {
        self.raw_sign_deterministic(msg, &[])
    }
}

/// The `Keypair` implementation for `SigningKey` allows deriving a `VerifyingKey` from a bare
/// `SigningKey` (even in the absence of the original seed).
impl<P: MlDsaParams> signature::Keypair for SigningKey<P> {
    type VerifyingKey = VerifyingKey<P>;

    fn verifying_key(&self) -> Self::VerifyingKey {
        let As1 = &self.A_hat * &self.s1_hat;
        let t = &As1.ntt_inverse() + &self.s2;

        // Discard t0
        let (t1, _) = t.power2round();

        VerifyingKey::new(self.rho.clone(), t1, Some(self.A_hat.clone()), None)
    }
}

/// The `RandomizedSigner` implementation for `SigningKey` only supports signing with an empty
/// context string. If you would like to include a context string, use the
/// [`SigningKey::sign_randomized`] method.
#[cfg(feature = "rand_core")]
impl<P: MlDsaParams> RandomizedSigner<Signature<P>> for SigningKey<P> {
    fn try_sign_with_rng<R: rand_core::TryCryptoRng + ?Sized>(
        &self,
        rng: &mut R,
        msg: &[u8],
    ) -> Result<Signature<P>, Error> {
        self.try_multipart_sign_with_rng(rng, &[msg])
    }
}

/// The `RandomizedMultipartSigner` implementation for `SigningKey` only supports signing with an
/// empty context string. If you would like to include a context string, use the
/// [`SigningKey::sign_randomized`] method.
#[cfg(feature = "rand_core")]
impl<P: MlDsaParams> RandomizedMultipartSigner<Signature<P>> for SigningKey<P> {
    fn try_multipart_sign_with_rng<R: rand_core::TryCryptoRng + ?Sized>(
        &self,
        rng: &mut R,
        msg: &[&[u8]],
    ) -> Result<Signature<P>, Error> {
        self.raw_sign_randomized(msg, &[], rng)
    }
}

// ============================================================================
// VerifyingKey
// ============================================================================

/// An ML-DSA verification key
#[derive(Clone, Debug, PartialEq)]
pub struct VerifyingKey<P: ParameterSet> {
    rho: B32,
    t1: Vector<P::K>,

    // Derived values
    A_hat: NttMatrix<P::K, P::L>,
    t1_2d_hat: NttVector<P::K>,
    tr: B64,
}

impl<P: MlDsaParams> VerifyingKey<P> {
    fn new(
        rho: B32,
        t1: Vector<P::K>,
        A_hat: Option<NttMatrix<P::K, P::L>>,
        enc: Option<EncodedVerifyingKey<P>>,
    ) -> Self {
        let A_hat = A_hat.unwrap_or_else(|| expand_a(&rho));
        let enc = enc.unwrap_or_else(|| Self::encode_internal(&rho, &t1));

        let t1_2d_hat = (Elem::new(1 << 13) * &t1).ntt();
        let mut tr_ctx = H::default();
        tr_ctx.absorb(&enc);
        let tr: B64 = tr_ctx.squeeze_new_array();

        Self {
            rho,
            t1,
            A_hat,
            t1_2d_hat,
            tr,
        }
    }

    /// This algorithm reflects the ML-DSA.Verify_internal algorithm from FIPS 204. It does not
    /// include the domain separator that distinguishes between the normal and pre-hashed cases,
    /// and it does not separate the context string from the rest of the message.
    // Algorithm 8 ML-DSA.Verify_internal
    pub fn verify_internal(&self, M: &[u8], sigma: &Signature<P>) -> bool {
        let mu = MuBuilder::internal(&self.tr, &[M]);
        self.raw_verify_mu(&mu, sigma)
    }

    fn raw_verify_mu(&self, mu: &B64, sigma: &Signature<P>) -> bool {
        // Reconstruct w
        let c = sample_in_ball(&sigma.c_tilde, P::TAU);

        let z_hat = sigma.z.ntt();
        let c_hat = c.ntt();
        let Az_hat = &self.A_hat * &z_hat;
        let ct1_2d_hat = &c_hat * &self.t1_2d_hat;

        let wp_approx = (&Az_hat - &ct1_2d_hat).ntt_inverse();
        let w1p = sigma.h.use_hint(&wp_approx);

        let w1p_tilde = P::encode_w1(&w1p);
        let mut cp_ctx = H::default();
        cp_ctx.absorb(mu);
        cp_ctx.absorb(&w1p_tilde);
        let cp_tilde: Array<u8, P::Lambda> = cp_ctx.squeeze_new_array();

        sigma.c_tilde == cp_tilde
    }

    /// This algorithm reflects the ML-DSA.Verify algorithm from FIPS 204.
    // Algorithm 3 ML-DSA.Verify
    pub fn verify_with_context(&self, M: &[u8], ctx: &[u8], sigma: &Signature<P>) -> bool {
        self.raw_verify_with_context(&[M], ctx, sigma)
    }

    fn raw_verify_with_context(&self, M: &[&[u8]], ctx: &[u8], sigma: &Signature<P>) -> bool {
        if ctx.len() > 255 {
            return false;
        }

        let mu = MuBuilder::new(&self.tr, ctx).message(M);
        self.raw_verify_mu(&mu, sigma)
    }

    fn encode_internal(rho: &B32, t1: &Vector<P::K>) -> EncodedVerifyingKey<P> {
        let t1_enc = P::encode_t1(t1);
        P::concat_vk(rho.clone(), t1_enc)
    }

    /// Encode the key in a fixed-size byte array.
    // Algorithm 22 pkEncode
    pub fn encode(&self) -> EncodedVerifyingKey<P> {
        Self::encode_internal(&self.rho, &self.t1)
    }

    /// Decode the key from an appropriately sized byte array.
    // Algorithm 23 pkDecode
    pub fn decode(enc: &EncodedVerifyingKey<P>) -> Self {
        let (rho, t1_enc) = P::split_vk(enc);
        let t1 = P::decode_t1(t1_enc);
        Self::new(rho.clone(), t1, None, Some(enc.clone()))
    }
}

impl<P: MlDsaParams> signature::Verifier<Signature<P>> for VerifyingKey<P> {
    fn verify(&self, msg: &[u8], signature: &Signature<P>) -> Result<(), Error> {
        self.multipart_verify(&[msg], signature)
    }
}

impl<P: MlDsaParams> MultipartVerifier<Signature<P>> for VerifyingKey<P> {
    fn multipart_verify(&self, msg: &[&[u8]], signature: &Signature<P>) -> Result<(), Error> {
        self.raw_verify_with_context(msg, &[], signature)
            .then_some(())
            .ok_or(Error::new())
    }
}

// ============================================================================
// KeyGen trait
// ============================================================================

/// A parameter set that knows how to generate key pairs
pub trait KeyGen: MlDsaParams {
    /// Generate a signing key pair from the specified RNG
    #[cfg(feature = "rand_core")]
    fn key_gen<R: CryptoRng + ?Sized>(rng: &mut R) -> KeyPair<Self>;

    /// Deterministically generate a signing key pair from the specified seed
    ///
    /// This method reflects the ML-DSA.KeyGen_internal algorithm from FIPS 204.
    fn from_seed(xi: &B32) -> KeyPair<Self>;
}

impl<P> KeyGen for P
where
    P: MlDsaParams,
{
    /// Generate a signing key pair from the specified RNG
    // Algorithm 1 ML-DSA.KeyGen()
    #[cfg(feature = "rand_core")]
    fn key_gen<R: CryptoRng + ?Sized>(rng: &mut R) -> KeyPair<P> {
        let mut xi = B32::default();
        rng.fill_bytes(&mut xi);
        Self::from_seed(&xi)
    }

    /// Deterministically generate a signing key pair from the specified seed
    // Algorithm 6 ML-DSA.KeyGen_internal
    fn from_seed(xi: &Seed) -> KeyPair<P> {
        // Derive seeds
        let mut h = H::default();
        h.absorb(xi);
        h.absorb(&[P::K::U8]);
        h.absorb(&[P::L::U8]);

        let rho: B32 = h.squeeze_new_array();
        let rhop: B64 = h.squeeze_new_array();
        let K: B32 = h.squeeze_new_array();

        // Sample private key components
        let A_hat = expand_a::<P::K, P::L>(&rho);
        let s1 = expand_s::<P::L>(&rhop, P::Eta::ETA, 0);
        let s2 = expand_s::<P::K>(&rhop, P::Eta::ETA, P::L::USIZE);

        // Compute derived values
        let As1_hat = &A_hat * &s1.ntt();
        let t = &As1_hat.ntt_inverse() + &s2;

        // Compress and encode
        let (t1, t0) = t.power2round();

        let verifying_key = VerifyingKey::new(rho.clone(), t1, Some(A_hat.clone()), None);
        let signing_key =
            SigningKey::new(rho, K, verifying_key.tr.clone(), s1, s2, t0, Some(A_hat));

        KeyPair {
            signing_key,
            verifying_key,
            seed: xi.clone(),
        }
    }
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod test {
    use super::*;
    use crate::param::*;

    #[test]
    fn output_sizes() {
        //           priv pub  sig
        // ML-DSA-44 2560 1312 2420
        // ML-DSA-65 4032 1952 3309
        // ML-DSA-87 4896 2592 4627
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

    fn encode_decode_round_trip_test<P>()
    where
        P: MlDsaParams + PartialEq,
    {
        let seed = Array::default();
        let kp = P::from_seed(&seed);
        assert_eq!(kp.to_seed(), seed);

        let sk = &kp.signing_key;
        let vk = &kp.verifying_key;

        let vk_bytes = vk.encode();
        let vk2 = VerifyingKey::<P>::decode(&vk_bytes);
        assert!(*vk == vk2);

        let sk_bytes = sk.to_expanded();
        let sk2 = SigningKey::<P>::from_expanded(&sk_bytes);
        assert!(*sk == sk2);

        let M = b"Hello world";
        let rnd = Array([0u8; 32]);
        let sig = sk.sign_internal(&[M], &rnd);
        let sig_bytes = sig.encode();
        let sig2 = Signature::<P>::decode(&sig_bytes).unwrap();
        assert!(sig == sig2);
    }

    #[test]
    fn encode_decode_round_trip() {
        encode_decode_round_trip_test::<MlDsa44>();
        encode_decode_round_trip_test::<MlDsa65>();
        encode_decode_round_trip_test::<MlDsa87>();
    }

    fn public_from_private_test<P>()
    where
        P: MlDsaParams + PartialEq,
    {
        let kp = P::from_seed(&Array::default());
        let sk = &kp.signing_key;
        let vk = &kp.verifying_key;
        let vk_derived = sk.verifying_key();

        assert!(*vk == vk_derived);
    }

    #[test]
    fn public_from_private() {
        public_from_private_test::<MlDsa44>();
        public_from_private_test::<MlDsa65>();
        public_from_private_test::<MlDsa87>();
    }

    fn sign_verify_round_trip_test<P>()
    where
        P: MlDsaParams,
    {
        let kp = P::from_seed(&Array::default());
        let sk = &kp.signing_key;
        let vk = &kp.verifying_key;

        let M = b"Hello world";
        let rnd = Array([0u8; 32]);
        let sig = sk.sign_internal(&[M], &rnd);

        assert!(vk.verify_internal(M, &sig));
    }

    #[test]
    fn sign_verify_round_trip() {
        sign_verify_round_trip_test::<MlDsa44>();
        sign_verify_round_trip_test::<MlDsa65>();
        sign_verify_round_trip_test::<MlDsa87>();
    }

    #[test]
    fn from_seed_implementations_match() {
        fn assert_from_seed_equality<P>()
        where
            P: MlDsaParams,
        {
            let seed = Seed::default();
            let kp1 = P::from_seed(&seed);
            let sk1 = SigningKey::<P>::from_seed(&seed);
            let vk1 = sk1.verifying_key();
            assert_eq!(kp1.signing_key, sk1);
            assert_eq!(kp1.verifying_key, vk1);
        }
        assert_from_seed_equality::<MlDsa44>();
        assert_from_seed_equality::<MlDsa65>();
        assert_from_seed_equality::<MlDsa87>();
    }

    #[test]
    fn to_seed_returns_correct_seed() {
        fn test_to_seed<P: MlDsaParams>() {
            let seed = Array([
                1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23,
                24, 25, 26, 27, 28, 29, 30, 31, 32,
            ]);
            let kp = P::from_seed(&seed);
            assert_eq!(kp.to_seed(), seed);
        }
        test_to_seed::<MlDsa44>();
        test_to_seed::<MlDsa65>();
        test_to_seed::<MlDsa87>();
    }

    #[test]
    fn verification_rejects_invalid_signature() {
        fn test_invalid_sig<P: MlDsaParams>() {
            let kp = P::from_seed(&Array::default());
            let vk = kp.verifying_key();

            let msg = b"Hello world";
            let rnd = Array([0u8; 32]);
            let mut sig = kp.signing_key().sign_internal(&[msg], &rnd);
            sig.c_tilde[0] ^= 0xFF;

            assert!(!vk.verify_with_context(msg, &[], &sig));
        }
        test_invalid_sig::<MlDsa44>();
        test_invalid_sig::<MlDsa65>();
        test_invalid_sig::<MlDsa87>();
    }

    #[test]
    fn verification_rejects_wrong_message() {
        fn test_wrong_msg<P: MlDsaParams>() {
            let kp = P::from_seed(&Array::default());
            let vk = kp.verifying_key();

            let msg1 = b"Hello world";
            let msg2 = b"Wrong message";
            let rnd = Array([0u8; 32]);
            let sig = kp.signing_key().sign_internal(&[msg1], &rnd);

            assert!(!vk.verify_with_context(msg2, &[], &sig));
        }
        test_wrong_msg::<MlDsa44>();
        test_wrong_msg::<MlDsa65>();
        test_wrong_msg::<MlDsa87>();
    }

    #[test]
    fn context_length_validation() {
        fn test_ctx_length<P: MlDsaParams>() {
            let kp = P::from_seed(&Array::default());
            let sk = kp.signing_key();
            let vk = kp.verifying_key();

            let msg = b"Hello world";
            let long_ctx = [0u8; 256];
            let short_ctx = [0u8; 255];

            assert!(sk.sign_deterministic(msg, &long_ctx).is_err());

            let sig = sk.sign_deterministic(msg, &short_ctx).unwrap();
            assert!(!vk.verify_with_context(msg, &long_ctx, &sig));
            assert!(vk.verify_with_context(msg, &short_ctx, &sig));
        }
        test_ctx_length::<MlDsa44>();
        test_ctx_length::<MlDsa65>();
        test_ctx_length::<MlDsa87>();
    }

    #[test]
    fn derived_verifying_key_validates_signatures() {
        fn test_derived_vk<P: MlDsaParams>() {
            let seed = Array([42u8; 32]);
            let kp = P::from_seed(&seed);
            let sk = kp.signing_key();
            let derived_vk = sk.verifying_key();

            let msg = b"Test message for derived key";
            let rnd = Array([0u8; 32]);
            let sig = sk.sign_internal(&[msg], &rnd);

            assert!(derived_vk.verify_internal(msg, &sig));
            assert_eq!(derived_vk.encode(), kp.verifying_key().encode());
        }
        test_derived_vk::<MlDsa44>();
        test_derived_vk::<MlDsa65>();
        test_derived_vk::<MlDsa87>();
    }

    // ========================================================================
    // Signature trait integration tests
    // ========================================================================

    #[test]
    fn signer_verifier_trait_round_trip() {
        use signature::{Signer, Verifier};

        fn test_signer_verifier<P: MlDsaParams>() {
            let kp = P::from_seed(&Array::default());
            let sk = kp.signing_key();
            let vk = kp.verifying_key();

            let msg = b"Hello world";
            let sig: Signature<P> = sk.sign(msg);
            assert!(vk.verify(msg, &sig).is_ok());
            assert!(vk.verify(b"Wrong message", &sig).is_err());
        }
        test_signer_verifier::<MlDsa44>();
        test_signer_verifier::<MlDsa65>();
        test_signer_verifier::<MlDsa87>();
    }

    #[test]
    fn keypair_signer_trait() {
        use signature::{Signer, Verifier};

        fn test_keypair_signer<P: MlDsaParams>() {
            let kp = P::from_seed(&Array::default());

            let msg = b"Signed by keypair";
            let sig: Signature<P> = kp.sign(msg);
            assert!(kp.verifying_key().verify(msg, &sig).is_ok());
        }
        test_keypair_signer::<MlDsa44>();
        test_keypair_signer::<MlDsa65>();
        test_keypair_signer::<MlDsa87>();
    }

    #[test]
    fn keypair_ref_trait() {
        fn test_keypair_ref<P: MlDsaParams>() {
            let kp = P::from_seed(&Array::default());
            let vk_from_trait: VerifyingKey<P> = signature::Keypair::verifying_key(&kp);
            assert_eq!(vk_from_trait.encode(), kp.verifying_key().encode());
        }
        test_keypair_ref::<MlDsa44>();
        test_keypair_ref::<MlDsa65>();
        test_keypair_ref::<MlDsa87>();
    }

    #[test]
    fn signing_key_keypair_trait() {
        use signature::{Keypair, Signer, Verifier};

        fn test_sk_keypair<P: MlDsaParams>() {
            let kp = P::from_seed(&Array([42u8; 32]));
            let sk = kp.signing_key();

            // Derive verifying key via Keypair trait
            let vk: VerifyingKey<P> = Keypair::verifying_key(sk);
            assert_eq!(vk.encode(), kp.verifying_key().encode());

            // Sign and verify using derived key
            let msg = b"Keypair trait test";
            let sig: Signature<P> = sk.sign(msg);
            assert!(vk.verify(msg, &sig).is_ok());
        }
        test_sk_keypair::<MlDsa44>();
        test_sk_keypair::<MlDsa65>();
        test_sk_keypair::<MlDsa87>();
    }

    #[test]
    fn multipart_signer_verifier_trait() {
        use signature::{MultipartSigner, MultipartVerifier};

        fn test_multipart<P: MlDsaParams>() {
            let kp = P::from_seed(&Array::default());
            let sk = kp.signing_key();
            let vk = kp.verifying_key();

            let parts: &[&[u8]] = &[b"Hello", b" ", b"world"];
            let sig: Signature<P> = sk.multipart_sign(parts);
            assert!(vk.multipart_verify(parts, &sig).is_ok());

            // Different parts should fail
            let wrong_parts: &[&[u8]] = &[b"Hello", b" ", b"wrong"];
            assert!(vk.multipart_verify(wrong_parts, &sig).is_err());
        }
        test_multipart::<MlDsa44>();
        test_multipart::<MlDsa65>();
        test_multipart::<MlDsa87>();
    }

    #[test]
    fn signature_encoding_trait() {
        use signature::SignatureEncoding;

        fn test_encoding<P: MlDsaParams>() {
            let kp = P::from_seed(&Array::default());
            let sk = kp.signing_key();

            let msg = b"Encoding test";
            let sig: Signature<P> = sk.sign_deterministic(msg, &[]).unwrap();

            // to_bytes round-trip
            let bytes = sig.to_bytes();
            let sig2 = Signature::<P>::try_from(bytes.as_ref()).unwrap();
            assert_eq!(sig, sig2);
        }
        test_encoding::<MlDsa44>();
        test_encoding::<MlDsa65>();
        test_encoding::<MlDsa87>();
    }

    #[test]
    fn randomized_signer_trait() {
        use rand::rngs::StdRng;
        use rand::SeedableRng;
        use signature::{RandomizedSigner, Verifier};

        fn test_randomized<P: MlDsaParams>() {
            let kp = P::from_seed(&Array::default());
            let sk = kp.signing_key();
            let vk = kp.verifying_key();

            let mut rng = StdRng::seed_from_u64(42);
            let msg = b"Randomized signing";
            let sig: Signature<P> = sk.sign_with_rng(&mut rng, msg);
            assert!(vk.verify(msg, &sig).is_ok());
        }
        test_randomized::<MlDsa44>();
        test_randomized::<MlDsa65>();
        test_randomized::<MlDsa87>();
    }
}
