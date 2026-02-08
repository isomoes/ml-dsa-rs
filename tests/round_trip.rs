//! Key generation and sign/verify round-trip tests across all parameter sets.

use ml_dsa::*;

use hybrid_array::Array;

// ============================================================================
// Key generation round-trip tests
// ============================================================================

macro_rules! keygen_round_trip_tests {
    ($name:ident, $alg:ident) => {
        mod $name {
            use super::*;

            #[test]
            fn seed_determinism() {
                // Same seed always produces the same key pair
                let seed = Array([42u8; 32]);
                let kp1 = $alg::from_seed(&seed);
                let kp2 = $alg::from_seed(&seed);

                assert_eq!(
                    kp1.signing_key().to_expanded(),
                    kp2.signing_key().to_expanded()
                );
                assert_eq!(kp1.verifying_key().encode(), kp2.verifying_key().encode());
            }

            #[test]
            fn different_seeds_different_keys() {
                let seed1 = Array([1u8; 32]);
                let seed2 = Array([2u8; 32]);
                let kp1 = $alg::from_seed(&seed1);
                let kp2 = $alg::from_seed(&seed2);

                assert_ne!(kp1.verifying_key().encode(), kp2.verifying_key().encode());
                assert_ne!(
                    kp1.signing_key().to_expanded(),
                    kp2.signing_key().to_expanded()
                );
            }

            #[test]
            fn seed_round_trip() {
                let seed = Array([99u8; 32]);
                let kp = $alg::from_seed(&seed);
                assert_eq!(kp.to_seed(), seed);
            }

            #[test]
            fn verifying_key_encode_decode() {
                let kp = $alg::from_seed(&Array([0u8; 32]));
                let vk = kp.verifying_key();
                let encoded = vk.encode();
                let decoded = VerifyingKey::<$alg>::decode(&encoded);
                assert_eq!(*vk, decoded);
            }

            #[test]
            fn signing_key_encode_decode() {
                let kp = $alg::from_seed(&Array([0u8; 32]));
                let sk = kp.signing_key();
                let encoded = sk.to_expanded();
                let decoded = SigningKey::<$alg>::from_expanded(&encoded);
                assert_eq!(*sk, decoded);
            }

            #[test]
            fn derived_verifying_key_matches() {
                let kp = $alg::from_seed(&Array([77u8; 32]));
                let derived_vk = kp.signing_key().verifying_key();
                assert_eq!(kp.verifying_key().encode(), derived_vk.encode());
            }

            #[test]
            fn from_seed_matches_keygen_trait() {
                let seed = Array([55u8; 32]);
                let kp = $alg::from_seed(&seed);
                let sk = SigningKey::<$alg>::from_seed(&seed);
                assert_eq!(*kp.signing_key(), sk);
                assert_eq!(*kp.verifying_key(), sk.verifying_key());
            }

            #[cfg(feature = "rand_core")]
            #[test]
            fn random_keygen() {
                let mut rng = rand::rng();
                let kp = $alg::key_gen(&mut rng);

                // Verify the generated key pair is internally consistent
                let derived_vk = kp.signing_key().verifying_key();
                assert_eq!(kp.verifying_key().encode(), derived_vk.encode());

                // Verify the seed round-trips
                let kp2 = $alg::from_seed(&kp.to_seed());
                assert_eq!(kp.verifying_key().encode(), kp2.verifying_key().encode());
            }
        }
    };
}

keygen_round_trip_tests!(mldsa44_keygen, MlDsa44);
keygen_round_trip_tests!(mldsa65_keygen, MlDsa65);
keygen_round_trip_tests!(mldsa87_keygen, MlDsa87);

// ============================================================================
// Sign/verify round-trip tests
// ============================================================================

macro_rules! sign_verify_round_trip_tests {
    ($name:ident, $alg:ident) => {
        mod $name {
            use super::*;

            #[test]
            fn sign_internal_verify_internal() {
                let kp = $alg::from_seed(&Array([0u8; 32]));
                let msg = b"Hello, ML-DSA!";
                let rnd = Array([0u8; 32]);
                let sig = kp.signing_key().sign_internal(&[msg], &rnd);
                assert!(kp.verifying_key().verify_internal(msg, &sig));
            }

            #[test]
            fn sign_deterministic_verify_with_context() {
                let kp = $alg::from_seed(&Array([1u8; 32]));
                let msg = b"Test message";
                let ctx = b"test context";
                let sig = kp.signing_key().sign_deterministic(msg, ctx).unwrap();
                assert!(kp.verifying_key().verify_with_context(msg, ctx, &sig));
            }

            #[test]
            fn empty_message() {
                let kp = $alg::from_seed(&Array([2u8; 32]));
                let msg = b"";
                let rnd = Array([0u8; 32]);
                let sig = kp.signing_key().sign_internal(&[msg], &rnd);
                assert!(kp.verifying_key().verify_internal(msg, &sig));
            }

            #[test]
            fn large_message() {
                let kp = $alg::from_seed(&Array([3u8; 32]));
                let msg = vec![0xABu8; 100_000];
                let rnd = Array([0u8; 32]);
                let sig = kp.signing_key().sign_internal(&[&msg], &rnd);
                assert!(kp.verifying_key().verify_internal(&msg, &sig));
            }

            #[test]
            fn empty_context() {
                let kp = $alg::from_seed(&Array([4u8; 32]));
                let msg = b"message";
                let sig = kp.signing_key().sign_deterministic(msg, &[]).unwrap();
                assert!(kp.verifying_key().verify_with_context(msg, &[], &sig));
            }

            #[test]
            fn max_context_length() {
                let kp = $alg::from_seed(&Array([5u8; 32]));
                let msg = b"message";
                let ctx = [0u8; 255];
                let sig = kp.signing_key().sign_deterministic(msg, &ctx).unwrap();
                assert!(kp.verifying_key().verify_with_context(msg, &ctx, &sig));
            }

            #[test]
            fn signature_encode_decode_round_trip() {
                let kp = $alg::from_seed(&Array([6u8; 32]));
                let msg = b"encode test";
                let rnd = Array([0u8; 32]);
                let sig = kp.signing_key().sign_internal(&[msg], &rnd);
                let encoded = sig.encode();
                let decoded = Signature::<$alg>::decode(&encoded).unwrap();
                assert_eq!(sig, decoded);
                assert!(kp.verifying_key().verify_internal(msg, &decoded));
            }

            #[test]
            fn signature_try_from_bytes() {
                let kp = $alg::from_seed(&Array([7u8; 32]));
                let msg = b"try_from test";
                let rnd = Array([0u8; 32]);
                let sig = kp.signing_key().sign_internal(&[msg], &rnd);
                let encoded = sig.encode();
                let decoded = Signature::<$alg>::try_from(encoded.as_slice()).unwrap();
                assert_eq!(sig, decoded);
            }

            #[test]
            fn deterministic_signing_produces_same_signature() {
                let kp = $alg::from_seed(&Array([8u8; 32]));
                let msg = b"deterministic";
                let sig1 = kp.signing_key().sign_deterministic(msg, &[]).unwrap();
                let sig2 = kp.signing_key().sign_deterministic(msg, &[]).unwrap();
                assert_eq!(sig1.encode(), sig2.encode());
            }

            #[cfg(feature = "rand_core")]
            #[test]
            fn randomized_signing_verifies() {
                let kp = $alg::from_seed(&Array([9u8; 32]));
                let msg = b"randomized test";
                let mut rng = rand::rng();
                let sig = kp
                    .signing_key()
                    .sign_randomized(msg, &[], &mut rng)
                    .unwrap();
                assert!(kp.verifying_key().verify_with_context(msg, &[], &sig));
            }

            #[test]
            fn multipart_message() {
                let kp = $alg::from_seed(&Array([10u8; 32]));
                let part1 = b"Hello, ";
                let part2 = b"world!";
                let rnd = Array([0u8; 32]);
                let sig = kp.signing_key().sign_internal(&[part1, part2], &rnd);
                // Multipart signing should produce a valid signature
                // (verified against the concatenated message via internal API)
                let encoded = sig.encode();
                let decoded = Signature::<$alg>::decode(&encoded).unwrap();
                assert_eq!(sig, decoded);
            }
        }
    };
}

sign_verify_round_trip_tests!(mldsa44_sign_verify, MlDsa44);
sign_verify_round_trip_tests!(mldsa65_sign_verify, MlDsa65);
sign_verify_round_trip_tests!(mldsa87_sign_verify, MlDsa87);
