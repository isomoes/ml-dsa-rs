//! Boundary condition tests: invalid signatures, wrong keys, edge cases.

use ml_dsa::*;

use hybrid_array::Array;

// ============================================================================
// Invalid signature tests
// ============================================================================

macro_rules! boundary_tests {
    ($name:ident, $alg:ident) => {
        mod $name {
            use super::*;

            #[test]
            fn tampered_c_tilde_rejects() {
                let kp = $alg::from_seed(&Array([0u8; 32]));
                let msg = b"test message";
                let rnd = Array([0u8; 32]);
                let sig = kp.signing_key().sign_internal(&[msg], &rnd);
                let mut sig_bytes = sig.encode();

                // Flip a bit in c_tilde (first bytes of the signature)
                sig_bytes[0] ^= 0xFF;

                if let Some(tampered_sig) = Signature::<$alg>::decode(&sig_bytes) {
                    assert!(
                        !kp.verifying_key().verify_internal(msg, &tampered_sig),
                        "Tampered c_tilde should not verify"
                    );
                }
                // If decode fails, that's also acceptable (invalid signature)
            }

            #[test]
            fn tampered_z_rejects() {
                let kp = $alg::from_seed(&Array([0u8; 32]));
                let msg = b"test message";
                let rnd = Array([0u8; 32]);
                let sig = kp.signing_key().sign_internal(&[msg], &rnd);
                let mut sig_bytes = sig.encode();

                // Tamper with z portion (after c_tilde, which is Lambda bytes)
                let lambda_bytes = sig_bytes.len() / 4; // approximate
                if lambda_bytes + 10 < sig_bytes.len() {
                    sig_bytes[lambda_bytes + 5] ^= 0xFF;
                }

                if let Some(tampered_sig) = Signature::<$alg>::decode(&sig_bytes) {
                    assert!(
                        !kp.verifying_key().verify_internal(msg, &tampered_sig),
                        "Tampered z should not verify"
                    );
                }
            }

            #[test]
            fn tampered_hint_rejects() {
                let kp = $alg::from_seed(&Array([0u8; 32]));
                let msg = b"test message";
                let rnd = Array([0u8; 32]);
                let sig = kp.signing_key().sign_internal(&[msg], &rnd);
                let mut sig_bytes = sig.encode();

                // Tamper with hint portion (last bytes of the signature)
                let last = sig_bytes.len() - 1;
                sig_bytes[last] ^= 0xFF;

                if let Some(tampered_sig) = Signature::<$alg>::decode(&sig_bytes) {
                    assert!(
                        !kp.verifying_key().verify_internal(msg, &tampered_sig),
                        "Tampered hint should not verify"
                    );
                }
            }

            #[test]
            fn wrong_verifying_key_rejects() {
                let kp1 = $alg::from_seed(&Array([1u8; 32]));
                let kp2 = $alg::from_seed(&Array([2u8; 32]));
                let msg = b"test message";
                let rnd = Array([0u8; 32]);
                let sig = kp1.signing_key().sign_internal(&[msg], &rnd);

                assert!(
                    !kp2.verifying_key().verify_internal(msg, &sig),
                    "Signature should not verify with wrong key"
                );
            }

            #[test]
            fn wrong_message_rejects() {
                let kp = $alg::from_seed(&Array([0u8; 32]));
                let msg1 = b"correct message";
                let msg2 = b"wrong message";
                let rnd = Array([0u8; 32]);
                let sig = kp.signing_key().sign_internal(&[msg1], &rnd);

                assert!(
                    !kp.verifying_key().verify_internal(msg2, &sig),
                    "Signature should not verify with wrong message"
                );
            }

            #[test]
            fn wrong_context_rejects() {
                let kp = $alg::from_seed(&Array([0u8; 32]));
                let msg = b"test message";
                let ctx1 = b"context1";
                let ctx2 = b"context2";
                let sig = kp.signing_key().sign_deterministic(msg, ctx1).unwrap();

                assert!(
                    !kp.verifying_key().verify_with_context(msg, ctx2, &sig),
                    "Signature should not verify with wrong context"
                );
            }

            #[test]
            fn context_too_long_errors() {
                let kp = $alg::from_seed(&Array([0u8; 32]));
                let msg = b"test message";
                let long_ctx = [0u8; 256];

                assert!(
                    kp.signing_key().sign_deterministic(msg, &long_ctx).is_err(),
                    "Context > 255 bytes should error"
                );
            }

            #[test]
            fn context_too_long_verify_rejects() {
                let kp = $alg::from_seed(&Array([0u8; 32]));
                let msg = b"test message";
                let valid_ctx = [0u8; 255];
                let long_ctx = [0u8; 256];

                let sig = kp
                    .signing_key()
                    .sign_deterministic(msg, &valid_ctx)
                    .unwrap();

                assert!(
                    !kp.verifying_key().verify_with_context(msg, &long_ctx, &sig),
                    "Verification with context > 255 bytes should fail"
                );
            }

            #[test]
            fn all_zero_signature_rejects() {
                let kp = $alg::from_seed(&Array([0u8; 32]));
                let msg = b"test message";
                let zero_sig = EncodedSignature::<$alg>::default();

                if let Some(sig) = Signature::<$alg>::decode(&zero_sig) {
                    assert!(
                        !kp.verifying_key().verify_internal(msg, &sig),
                        "All-zero signature should not verify"
                    );
                }
                // If decode fails, that's also acceptable
            }

            #[test]
            fn signature_wrong_length_rejects() {
                // Too short
                let short_bytes = vec![0u8; 10];
                assert!(
                    Signature::<$alg>::try_from(short_bytes.as_slice()).is_err(),
                    "Too-short bytes should fail to decode"
                );

                // Too long
                let sig_size = std::mem::size_of::<EncodedSignature<$alg>>();
                let long_bytes = vec![0u8; sig_size + 1];
                assert!(
                    Signature::<$alg>::try_from(long_bytes.as_slice()).is_err(),
                    "Too-long bytes should fail to decode"
                );
            }

            #[test]
            fn cross_parameter_set_rejects() {
                // Sign with one parameter set, try to verify with another
                // This is a compile-time check in Rust (different types), but we can
                // verify that the encoded bytes don't accidentally work across types
                let kp44 = MlDsa44::from_seed(&Array([0u8; 32]));
                let msg = b"cross-param test";
                let rnd = Array([0u8; 32]);
                let sig44 = kp44.signing_key().sign_internal(&[msg], &rnd);
                let sig44_bytes = sig44.encode();

                // The signature sizes differ between parameter sets, so trying to
                // decode a MlDsa44 signature as MlDsa65 should fail
                let sig44_slice = sig44_bytes.as_slice();
                // This should fail because the sizes don't match
                if sig44_slice.len() != std::mem::size_of::<EncodedSignature<MlDsa65>>() {
                    assert!(
                        Signature::<MlDsa65>::try_from(sig44_slice).is_err(),
                        "Cross-parameter-set signature should fail to decode"
                    );
                }
            }

            #[test]
            fn single_bit_flip_in_message_rejects() {
                let kp = $alg::from_seed(&Array([0u8; 32]));
                let msg = b"bit flip test message";
                let rnd = Array([0u8; 32]);
                let sig = kp.signing_key().sign_internal(&[msg], &rnd);

                let mut modified_msg = msg.to_vec();
                modified_msg[0] ^= 1; // Flip single bit

                assert!(
                    !kp.verifying_key().verify_internal(&modified_msg, &sig),
                    "Single bit flip in message should cause verification failure"
                );
            }

            #[test]
            fn verify_with_derived_key() {
                let kp = $alg::from_seed(&Array([0u8; 32]));
                let msg = b"derived key test";
                let rnd = Array([0u8; 32]);
                let sig = kp.signing_key().sign_internal(&[msg], &rnd);

                // Derive verifying key from signing key
                let derived_vk = kp.signing_key().verifying_key();
                assert!(derived_vk.verify_internal(msg, &sig));
            }

            #[test]
            fn verify_with_decoded_key() {
                let kp = $alg::from_seed(&Array([0u8; 32]));
                let msg = b"decoded key test";
                let rnd = Array([0u8; 32]);
                let sig = kp.signing_key().sign_internal(&[msg], &rnd);

                // Encode and decode the verifying key
                let vk_bytes = kp.verifying_key().encode();
                let decoded_vk = VerifyingKey::<$alg>::decode(&vk_bytes);
                assert!(decoded_vk.verify_internal(msg, &sig));
            }

            #[test]
            fn sign_with_decoded_key() {
                let kp = $alg::from_seed(&Array([0u8; 32]));
                let msg = b"decoded sk test";
                let rnd = Array([0u8; 32]);

                // Encode and decode the signing key
                let sk_bytes = kp.signing_key().to_expanded();
                let decoded_sk = SigningKey::<$alg>::from_expanded(&sk_bytes);
                let sig = decoded_sk.sign_internal(&[msg], &rnd);

                assert!(kp.verifying_key().verify_internal(msg, &sig));
            }
        }
    };
}

boundary_tests!(mldsa44_boundary, MlDsa44);
boundary_tests!(mldsa65_boundary, MlDsa65);
boundary_tests!(mldsa87_boundary, MlDsa87);
