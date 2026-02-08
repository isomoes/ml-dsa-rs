//! Cryptographic functions for ML-DSA

use hybrid_array::{Array, ArraySize};
use sha3::{
    digest::{ExtendableOutput, Update, XofReader},
    Shake128, Shake256,
};

/// State machine wrapper around a SHAKE instance.
///
/// SHAKE supports two phases:
/// - absorbing input bytes
/// - squeezing an arbitrary amount of output bytes
///
/// Once squeezing starts, additional absorption is not allowed.
pub enum ShakeState<Shake>
where
    Shake: Default + Update + ExtendableOutput,
{
    /// The state is still absorbing input bytes.
    Absorbing(Shake),
    /// The state is squeezing output bytes.
    Squeezing(<Shake as ExtendableOutput>::Reader),
}

impl<Shake> Default for ShakeState<Shake>
where
    Shake: Default + Update + ExtendableOutput,
{
    fn default() -> Self {
        Self::Absorbing(Shake::default())
    }
}

impl<Shake> ShakeState<Shake>
where
    Shake: Default + Update + ExtendableOutput,
{
    /// Absorb input bytes into the SHAKE state.
    ///
    /// # Panics
    ///
    /// Panics if called after squeezing has started.
    pub fn absorb(&mut self, input: &[u8]) {
        match self {
            Self::Absorbing(shake) => shake.update(input),
            Self::Squeezing(_) => panic!("cannot absorb after squeezing has started"),
        }
    }

    /// Squeeze output bytes from the SHAKE state.
    ///
    /// The first call transitions from absorbing to squeezing state.
    pub fn squeeze(&mut self, output: &mut [u8]) {
        let reader = match self {
            Self::Squeezing(reader) => reader,
            Self::Absorbing(_) => {
                let old_state = core::mem::replace(self, Self::Absorbing(Shake::default()));
                let reader = match old_state {
                    Self::Absorbing(shake) => shake.finalize_xof(),
                    Self::Squeezing(_) => unreachable!("replacement guarantees absorbing state"),
                };
                *self = Self::Squeezing(reader);

                match self {
                    Self::Squeezing(reader) => reader,
                    Self::Absorbing(_) => unreachable!("state set to squeezing"),
                }
            }
        };

        reader.read(output);
    }

    /// Squeeze a newly allocated fixed-size output array.
    pub fn squeeze_new<const N: usize>(&mut self) -> [u8; N] {
        let mut output = [0u8; N];
        self.squeeze(&mut output);
        output
    }

    /// Squeeze a newly allocated `hybrid_array::Array` of the given typenum size.
    pub fn squeeze_new_array<N: ArraySize>(&mut self) -> Array<u8, N> {
        let mut output = Array::<u8, N>::default();
        self.squeeze(&mut output);
        output
    }
}

/// SHAKE128 state used in ML-DSA (`G` in FIPS 204).
pub type G = ShakeState<Shake128>;

/// SHAKE256 state used in ML-DSA (`H` in FIPS 204).
pub type H = ShakeState<Shake256>;

#[cfg(test)]
mod tests {
    use super::{G, H};
    use hex_literal::hex;

    #[test]
    fn shake128_empty_message_known_vector() {
        let mut g = G::default();
        g.absorb(b"");

        let out = g.squeeze_new::<16>();
        assert_eq!(out, hex!("7f9c2ba4e88f827d616045507605853e"));
    }

    #[test]
    fn shake256_empty_message_known_vector() {
        let mut h = H::default();
        h.absorb(b"");

        let out = h.squeeze_new::<64>();
        assert_eq!(
            out,
            hex!(
                "46b9dd2b0ba88d13233b3feb743eeb243fcd52ea62b81b82"
                "b50c27646ed5762f"
                "d75dc4ddd8c0f200cb05019d67b592f6"
                "fc821c49479ab48640292eacb3b7c4be"
            )
        );
    }

    #[test]
    fn squeeze_is_continuous_across_calls() {
        let mut g = G::default();
        g.absorb(b"abc");

        let full = g.squeeze_new::<64>();

        let mut g_split = G::default();
        g_split.absorb(b"abc");
        let mut first = [0u8; 32];
        let mut second = [0u8; 32];
        g_split.squeeze(&mut first);
        g_split.squeeze(&mut second);

        assert_eq!(&full[..32], &first);
        assert_eq!(&full[32..], &second);
    }
}
