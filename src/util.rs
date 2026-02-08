//! Utility functions for ML-DSA

use hybrid_array::{
    typenum::{U32, U64},
    Array,
};

/// 32-byte array type.
pub type B32 = Array<u8, U32>;

/// 64-byte array type.
pub type B64 = Array<u8, U64>;
