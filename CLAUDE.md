# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is a pure Rust implementation of ML-DSA (Module-Lattice-Based Digital Signature Standard), formerly known as CRYSTALS-Dilithium, following FIPS 204 specification. The project structure is based on the RustCrypto ML-DSA implementation.

## Development Commands

### Building and Checking

```bash
cargo check        # Quick syntax and type checking
cargo build        # Build the library
cargo build --release  # Release build
```

### Testing

```bash
cargo test         # Run all tests (currently none implemented)
cargo test --doc   # Run documentation tests
```

### Benchmarking

```bash
cargo bench        # Run performance benchmarks
```

### Feature Testing

```bash
cargo check --no-default-features  # Check minimal build
cargo check --all-features         # Check with all features
```

## Code Quality

The project enforces:

- `#![forbid(unsafe_code)]` - No unsafe code allowed
- `#![warn(missing_docs, rust_2018_idioms)]` - Documentation and style warnings
- Current warnings exist for missing documentation on Q constants in param.rs

## Reference Implementation

This project follows the structure of the official RustCrypto ML-DSA implementation at https://github.com/RustCrypto/signatures/tree/master/ml-dsa
