//! Property-based tests for the `ml-dsa` crate.

macro_rules! signature_round_trip_encode {
    ($alg:ident, $sig:expr) => {{
        let sig_enc = $sig.encode();
        let sig_dec_result = Signature::<$alg>::decode(&sig_enc);
        prop_assert_eq!(&sig_dec_result, &Some($sig));
        sig_dec_result.unwrap()
    }};
}

macro_rules! mldsa_proptests {
    ($name:ident, $alg:ident) => {
        mod $name {
            use ml_dsa::{$alg, KeyGen, Signature};
            use proptest::{collection, prelude::*};

            proptest! {
                #[test]
                fn round_trip_test(
                    seed in any::<[u8; 32]>(),
                    msg in collection::vec(0u8..u8::MAX, 0..65536),
                    rnd in any::<[u8; 32]>()
                ) {
                    let kp = $alg::from_seed(&seed.into());
                    let sk = kp.signing_key();
                    let vk = kp.verifying_key();

                    let sig = sk.sign_internal(&[&msg], &rnd.into());
                    let sig_dec = signature_round_trip_encode!($alg, sig);
                    assert!(vk.verify_internal(&msg, &sig_dec));
                }

                #[test]
                fn deterministic_sign_verify(
                    seed in any::<[u8; 32]>(),
                    msg in collection::vec(0u8..u8::MAX, 0..1024),
                    ctx in collection::vec(0u8..u8::MAX, 0..255),
                ) {
                    let kp = $alg::from_seed(&seed.into());
                    let sk = kp.signing_key();
                    let vk = kp.verifying_key();

                    let sig = sk.sign_deterministic(&msg, &ctx).unwrap();
                    assert!(vk.verify_with_context(&msg, &ctx, &sig));
                }

                #[test]
                fn deterministic_signing_is_deterministic(
                    seed in any::<[u8; 32]>(),
                    msg in collection::vec(0u8..u8::MAX, 0..1024),
                ) {
                    let kp = $alg::from_seed(&seed.into());
                    let sk = kp.signing_key();

                    let sig1 = sk.sign_deterministic(&msg, &[]).unwrap();
                    let sig2 = sk.sign_deterministic(&msg, &[]).unwrap();
                    assert_eq!(sig1.encode(), sig2.encode());
                }

                #[test]
                fn key_encode_decode_round_trip(
                    seed in any::<[u8; 32]>(),
                ) {
                    let kp = $alg::from_seed(&seed.into());
                    let sk = kp.signing_key();
                    let vk = kp.verifying_key();

                    // Verifying key round-trip
                    let vk_enc = vk.encode();
                    let vk_dec = ml_dsa::VerifyingKey::<$alg>::decode(&vk_enc);
                    assert_eq!(*vk, vk_dec);

                    // Signing key round-trip
                    let sk_enc = sk.to_expanded();
                    let sk_dec = ml_dsa::SigningKey::<$alg>::from_expanded(&sk_enc);
                    assert_eq!(*sk, sk_dec);
                }

                #[test]
                fn seed_round_trip(
                    seed in any::<[u8; 32]>(),
                ) {
                    let kp = $alg::from_seed(&seed.into());
                    assert_eq!(kp.to_seed().as_slice(), &seed);
                }
            }
        }
    };
}

mldsa_proptests!(mldsa44, MlDsa44);
mldsa_proptests!(mldsa65, MlDsa65);
mldsa_proptests!(mldsa87, MlDsa87);
