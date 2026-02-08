# ML-DSA Implementation Roadmap

This document outlines the implementation plan for ML-DSA (FIPS 204) based on the [RustCrypto ML-DSA](https://github.com/RustCrypto/signatures/tree/master/ml-dsa) reference implementation.

## Current Status

The `module_lattice` foundation layer is complete. Track A (core algebra) is complete. Remaining ML-DSA-specific modules are stubs.

| Module                   | Status      | Notes                                                                          |
| ------------------------ | ----------- | ------------------------------------------------------------------------------ |
| `module_lattice/util`    | **Done**    | Truncate, Flatten, Unflatten                                                   |
| `module_lattice/algebra` | **Done**    | Field, Elem, Polynomial, Vector, NttPolynomial, NttVector, NttMatrix           |
| `module_lattice/encode`  | **Done**    | byte_encode/decode, Encode trait for Polynomial/Vector                         |
| `param`                  | **Done**    | Constants for ML-DSA-44/65/87 (raw constants, no traits yet)                   |
| `lib`                    | **Partial** | Signature/SigningKey/VerifyingKey type shells only                             |
| `algebra`                | **Done**    | BaseField, type aliases, BarrettReduce, ConstantTimeDiv, Decompose, AlgebraExt |
| `crypto`                 | **Done**    | ShakeState (`G`/`H`) + SHAKE known vectors                                     |
| `encode`                 | **Done**    | RangeEncodingSize, BitPack for Polynomial/Vector, Algorithm 17                 |
| `hint`                   | Stub        |                                                                                |
| `ntt`                    | **Done**    | NTT/NTT^-1, MultiplyNtt, ZETA_POW_BITREV table, round-trip tests              |
| `sampling`               | Stub        |                                                                                |
| `util`                   | Stub        |                                                                                |

## Parallel Task Plan

The dependency graph allows significant parallelism. Below, tasks at the same level can be worked on **concurrently**.

```
                    ┌─────────────────────────────────────────────┐
                    │  Phase 1: Foundation (DONE)                 │
                    │  module_lattice/{util, algebra, encode}     │
                    └──────────────────┬──────────────────────────┘
                                       │
              ┌────────────────────────┼────────────────────────┐
              ▼                        ▼                        ▼
   ┌──────────────────┐   ┌──────────────────┐   ┌──────────────────┐
   │ Track A: Algebra  │   │ Track B: Crypto  │   │ Track C: Util    │
   │ algebra.rs        │   │ crypto.rs        │   │ util.rs          │
   └────────┬─────────┘   └────────┬─────────┘   └──────────────────┘
            │                      │
     ┌──────┼──────┐               │
     ▼      ▼      ▼              ▼
   ┌────┐ ┌────┐ ┌────┐   ┌──────────────┐
   │NTT │ │Enc │ │Hint│   │  Sampling    │
   │ntt │ │enc │ │hint│   │  sampling.rs │
   └──┬─┘ └──┬─┘ └──┬─┘   └──────┬───────┘
      │       │      │            │
      └───────┴──────┴────────────┘
                     │
              ┌──────┴──────┐
              ▼             ▼
   ┌──────────────┐  ┌──────────────┐
   │ Param Traits │  │ Core Logic   │
   │ param.rs     │  │ lib.rs       │
   └──────────────┘  └──────┬───────┘
                            │
                     ┌──────┴──────┐
                     ▼             ▼
              ┌──────────┐  ┌──────────┐
              │ Trait     │  │ Testing  │
              │ Impls     │  │ & KATs   │
              └──────────┘  └──────────┘
```

---

## Level 0 — Foundation (DONE)

- [x] `module_lattice/util.rs` — Truncate, Flatten, Unflatten traits
- [x] `module_lattice/algebra.rs` — Field trait, define_field! macro, Elem/Polynomial/Vector/NttPolynomial/NttVector/NttMatrix with arithmetic
- [x] `module_lattice/encode.rs` — byte_encode/byte_decode, EncodingSize, Encode trait
- [x] `module_lattice/mod.rs` — Module re-exports

---

## Level 1 — Parallel Tracks A/B/C

These three tracks have **no dependencies on each other** and can be implemented concurrently.

### Track A: Core Algebra (`src/algebra.rs`) — DONE

**Depends on**: module_lattice (done)

- [x] Define `BaseField` using `define_field!(BaseField, u32, u64, u128, 8_380_417)`
- [x] Create type aliases: `Int`, `Elem`, `Polynomial`, `Vector<K>`, `NttPolynomial`, `NttVector<K>`, `NttMatrix<K, L>`
- [x] Implement `BarrettReduce` trait — Barrett modular reduction with precomputed constants
- [x] Implement `ConstantTimeDiv` trait — constant-time division via fixed-point multiplication
- [x] Implement `Decompose` trait for `Elem` — Algorithm 36 (Decompose)
- [x] Implement `AlgebraExt` trait:
  - [x] `mod_plus_minus<M>` — centered modular reduction
  - [x] `infinity_norm` — infinity norm for Elem, Polynomial, Vector
  - [x] `power2round` — Algorithm 35 (Power2Round)
  - [x] `high_bits<TwoGamma2>` — Algorithm 37 (HighBits)
  - [x] `low_bits<TwoGamma2>` — Algorithm 38 (LowBits)
- [x] Add unit tests for all algebra operations

### Track B: Cryptographic Primitives (`src/crypto.rs`)

**Depends on**: nothing (uses external sha3 crate)

- [x] Implement `ShakeState<Shake>` enum with `Absorbing`/`Squeezing` variants
- [x] Implement `absorb(&[u8])` method
- [x] Implement `squeeze(&mut [u8])` and `squeeze_new<N>()` methods
- [x] Define type aliases `G = ShakeState<Shake128>` and `H = ShakeState<Shake256>`
- [x] Add tests with known SHAKE test vectors

### Track C: Utility Types (`src/util.rs`)

**Depends on**: nothing

- [x] Define `B32 = Array<u8, U32>` byte array type
- [x] Define `B64 = Array<u8, U64>` byte array type
- [x] Add any additional utility functions as needed

---

## Level 2 — Parallel Tracks D/E/F/G

These tracks depend on Level 1 completions as noted, but are **independent of each other**.

### Track D: NTT (`src/ntt.rs`)

**Depends on**: Track A (algebra)

- [x] Define `ZETA_POW_BITREV` precomputed table (powers of zeta=1753, bit-reversed order)
- [x] Define `INVERSE_256` constant (8_347_681 mod Q)
- [x] Implement `ntt_layer` and `ntt_inverse_layer` helper functions
- [x] Implement `Ntt` trait — Algorithm 41 (NTT forward transform)
  - [x] For `Polynomial` → `NttPolynomial`
  - [x] For `Vector<K>` → `NttVector<K>`
- [x] Implement `NttInverse` trait — Algorithm 42 (NTT inverse transform)
  - [x] For `NttPolynomial` → `Polynomial`
  - [x] For `NttVector<K>` → `Vector<K>`
- [x] Implement `MultiplyNtt` for `BaseField` — Algorithm 45 (pointwise NTT multiplication)
- [x] Add NTT round-trip tests and verify against FIPS 204 Appendix B

### Track E: ML-DSA Encoding (`src/encode.rs`)

**Depends on**: Track A (algebra)

- [x] Implement `RangeEncodingSize` trait for `(A, B)` pairs
- [x] Define range encoding type aliases: `RangeMin`, `RangeMax`, `RangeEncodingBits`, `RangeEncodedPolynomial`, `RangeEncodedVector`
- [x] Implement `BitPack` trait for `Polynomial` — Algorithm 17 (BitPack/BitUnPack)
- [x] Implement `BitPack` trait for `Vector<K>`
- [x] Add encoding round-trip tests for various bit widths

### Track F: Hint Operations (`src/hint.rs`)

**Depends on**: Track A (algebra)

- [ ] Implement `make_hint<TwoGamma2>(z, r) -> bool` — Algorithm 39 (MakeHint)
- [ ] Implement `use_hint<TwoGamma2>(h, r) -> Elem` — Algorithm 40 (UseHint)
- [ ] Implement `Hint<P>` struct with `Array<Array<bool, U256>, P::K>` field
- [ ] Implement `Hint::new()` — create hint from two vectors
- [ ] Implement `Hint::hamming_weight()` — count set bits
- [ ] Implement `Hint::use_hint()` — apply hint to vector
- [ ] Implement `Hint::bit_pack()` / `Hint::bit_unpack()` — hint encoding/decoding
- [ ] Add hint tests

### Track G: Sampling Functions (`src/sampling.rs`)

**Depends on**: Track A (algebra) + Track B (crypto)

- [ ] Implement `bit_set` — Algorithm 13 (BytesToBits)
- [ ] Implement `coeff_from_three_bytes` — Algorithm 14
- [ ] Implement `coeff_from_half_byte` — Algorithm 15
- [ ] Implement `sample_in_ball` — Algorithm 29 (SampleInBall)
- [ ] Implement `rej_ntt_poly` — Algorithm 30 (RejNTTPoly) using `G` (SHAKE-128)
- [ ] Implement `rej_bounded_poly` — Algorithm 31 (RejBoundedPoly) using `H` (SHAKE-256)
- [ ] Implement `expand_a` — Algorithm 32 (ExpandA) → `NttMatrix<K, L>`
- [ ] Implement `expand_s` — Algorithm 33 (ExpandS) → `Vector<K>`
- [ ] Implement `expand_mask` — Algorithm 34 (ExpandMask) → `Vector<K>`
- [ ] Add sampling tests

---

## Level 3 — Parameter Traits (`src/param.rs`)

**Depends on**: Track A (algebra) + Track E (encode)

- [ ] Refactor raw constants into trait-based parameterization
- [ ] Define `ParameterSet` trait with associated types: `K`, `L`, `Eta`, `Gamma1`, `Gamma2`, `TwoGamma2`, `W1Bits`, `Lambda`, `Omega`, `TAU`, `BETA`
- [ ] Define `SamplingSize` trait with `ETA` enum (`Two`, `Four`)
- [ ] Define `MaskSamplingSize` trait with `SampleSize` and `unpack()`
- [ ] Implement `SigningKeyParams` trait (blanket impl):
  - [ ] `encode_s1/decode_s1`, `encode_s2/decode_s2`, `encode_t0/decode_t0`
  - [ ] `concat_sk/split_sk`
- [ ] Implement `VerifyingKeyParams` trait (blanket impl):
  - [ ] `encode_t1/decode_t1`, `concat_vk/split_vk`
- [ ] Implement `SignatureParams` trait (blanket impl):
  - [ ] `encode_w1/decode_w1`, `encode_z/decode_z`
  - [ ] `concat_sig/split_sig`, `split_hint`
- [ ] Define `MlDsaParams` super-trait combining all param traits
- [ ] Implement all traits for `MlDsa44`, `MlDsa65`, `MlDsa87`
- [ ] Add parameter validation tests

---

## Level 4 — Core ML-DSA Logic (`src/lib.rs`)

**Depends on**: All Level 2 tracks + Level 3

### Key Types

- [ ] Redefine `Signature<P>` with fields: `c_tilde`, `z`, `h`
- [ ] Redefine `SigningKey<P>` with fields: `rho`, `K`, `tr`, `s1`, `s2`, `t0` + derived NTT values (`s1_hat`, `s2_hat`, `t0_hat`, `A_hat`)
- [ ] Redefine `VerifyingKey<P>` with fields: `rho`, `t1` + derived `A_hat`, `t1_2d_hat`, `tr`
- [ ] Define `KeyPair<P>` struct
- [ ] Define `Seed = B32` type alias

### Key Generation

- [ ] Implement `KeyGen` trait — Algorithm 1 (ML-DSA.KeyGen)
- [ ] Implement `key_gen_internal` — Algorithm 6 (ML-DSA.KeyGen_internal)
- [ ] Implement `SigningKey::from_seed()`

### Signing

- [ ] Implement `sign_internal` — Algorithm 7 (ML-DSA.Sign_internal)
- [ ] Implement `sign_deterministic` — Algorithm 2 (deterministic variant)
- [ ] Implement `sign_randomized` — Algorithm 2 (randomized variant, feature-gated)
- [ ] Implement `MuBuilder` for domain separation

### Verification

- [ ] Implement `verify_internal` — Algorithm 8 (ML-DSA.Verify_internal)
- [ ] Implement `verify_with_context` — Algorithm 3 (ML-DSA.Verify)

### Key/Signature Encoding

- [ ] Implement `Signature::encode` / `Signature::decode` — Algorithms 26/27 (sigEncode/sigDecode)
- [ ] Implement `VerifyingKey::encode` / `VerifyingKey::decode` — Algorithms 22/23 (pkEncode/pkDecode)
- [ ] Implement `SigningKey::to_expanded` / `SigningKey::from_expanded` — Algorithms 24/25 (skEncode/skDecode)

---

## Level 5 — Parallel: Trait Impls + Testing

These two tracks can proceed **concurrently**.

### Track H: Signature Crate Integration

**Depends on**: Level 4

- [ ] Implement `signature::Signer` for `SigningKey<P>` and `KeyPair<P>`
- [ ] Implement `signature::Verifier` for `VerifyingKey<P>`
- [ ] Implement `signature::Keypair` for `KeyPair<P>`
- [ ] Implement `signature::SignatureEncoding` for `Signature<P>`
- [ ] Implement `MultipartSigner` / `MultipartVerifier`
- [ ] Implement `DigestSigner<Shake256>` / `DigestVerifier<Shake256>`
- [ ] Implement `RandomizedSigner` (feature-gated on `rand_core`)

### Track I: Testing & Validation

**Depends on**: Level 4

- [ ] Add FIPS 204 Known Answer Tests (key-gen.json, sig-gen.json, sig-ver.json)
- [ ] Add key generation round-trip tests
- [ ] Add sign/verify round-trip tests
- [ ] Add boundary condition tests (invalid signatures, wrong keys)
- [ ] Add interoperability tests against reference implementation
- [ ] Add property-based tests (proptest)
- [ ] Add Wycheproof test vectors

---

## Level 6 — Optional Features & Polish

### Track J: Optional Features

**Depends on**: Level 5

- [ ] Implement PKCS#8 support (`pkcs8` feature):
  - [ ] `AssociatedAlgorithmIdentifier` for parameter sets
  - [ ] `TryFrom<PrivateKeyInfoRef>` for `KeyPair`/`SigningKey`
  - [ ] `TryFrom<SubjectPublicKeyInfoRef>` for `VerifyingKey`
  - [ ] `EncodePrivateKey` / `EncodePublicKey` (with `alloc`)
  - [ ] `SignatureBitStringEncoding` for `Signature` (with `alloc`)
- [ ] Implement `Zeroize` / `ZeroizeOnDrop` support (`zeroize` feature)
- [ ] Add proper feature gates in `Cargo.toml`

### Track K: Documentation & Benchmarks

**Depends on**: Level 5

- [ ] Complete API documentation (resolve `missing_docs` warnings)
- [ ] Add usage examples in doc comments
- [ ] Add benchmarks for key generation, signing, verification
- [ ] Profile and optimize critical paths (NTT, sampling)
- [ ] Run clippy and fix all warnings
- [ ] Update README with features and usage

---

## Summary: Parallel Execution Plan

| Level | Tracks                                      | Can Run In Parallel                      |
| ----- | ------------------------------------------- | ---------------------------------------- |
| 0     | Foundation                                  | **DONE**                                 |
| 1     | A (algebra), B (crypto), C (util)           | **A DONE** ∥ B ∥ C                       |
| 2     | D (ntt), E (encode), F (hint), G (sampling) | D ∥ E ∥ F ∥ G (after respective L1 deps) |
| 3     | Param traits                                | Sequential (needs A + E)                 |
| 4     | Core logic (lib.rs)                         | Sequential (needs all above)             |
| 5     | H (trait impls), I (testing)                | H ∥ I                                    |
| 6     | J (optional features), K (docs/bench)       | J ∥ K                                    |

**Critical path**: A → D → Level 3 → Level 4 → Level 5

---

## Key Design Decisions

1. **Type-level Programming**: Use `typenum` for compile-time parameter verification
2. **No Unsafe Code**: Maintain `#![forbid(unsafe_code)]`
3. **Modular Design**: Keep components loosely coupled for testability
4. **Performance**: Use precomputed NTT tables and Barrett reduction
5. **Security**: Constant-time operations where applicable (ConstantTimeDiv, etc.)
6. **Compatibility**: Match upstream RustCrypto API surface for interoperability

## Success Criteria

- [ ] All ML-DSA-44/65/87 operations implemented correctly
- [ ] Passes all FIPS 204 test vectors (key-gen, sig-gen, sig-ver)
- [ ] Passes Wycheproof test vectors
- [ ] Compatible with RustCrypto `signature` traits
- [ ] No unsafe code
- [ ] Comprehensive test coverage (>90%)
- [ ] Performance comparable to reference implementation
