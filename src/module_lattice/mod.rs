//! This module contains functions that should be common across ML-KEM and ML-DSA:
//!
//! * Linear algebra with degree-256 polynomials over a prime-order field, vectors of such
//!   polynomials, and NTT polynomials / vectors.
//!
//! * Packing of polynomials into coefficients with a specified number of bits.
//!
//! * Utility functions such as truncating integers, flattening arrays of arrays, and unflattening
//!   arrays into arrays of arrays.
//!
//! While this is currently a module within the `ml_dsa` crate, the idea of pulling it out is that
//! it could be a separate crate on which both the `ml_dsa` crate and the `ml_kem` crate depend.

/// Utility functions such as truncating integers, flattening arrays of arrays, and unflattening
/// arrays into arrays of arrays.
pub mod util;

/// Linear algebra with degree-256 polynomials over a prime-order field, vectors of such
/// polynomials, and NTT polynomials / vectors
pub mod algebra;

/// Packing of polynomials into coefficients with a specified number of bits.
pub mod encode;
