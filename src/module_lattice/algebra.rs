use super::util::Truncate;

use core::fmt::Debug;
use core::ops::{Add, Mul, Neg, Sub};
use hybrid_array::{typenum::U256, Array, ArraySize};
use num_traits::PrimInt;

#[cfg(feature = "zeroize")]
use zeroize::Zeroize;

/// Trait defining a finite field with prime order
pub trait Field: Copy + Default + Debug + PartialEq {
    /// The primitive integer type used to represent field elements
    type Int: PrimInt + Default + Debug + From<u8> + Into<u128> + Into<Self::Long> + Truncate<u128>;
    /// Twice the bit width of Int, used for intermediate calculations
    type Long: PrimInt + From<Self::Int>;
    /// Four times the bit width of Int, used for Barrett reduction
    type LongLong: PrimInt;

    /// The prime modulus defining the field
    const Q: Self::Int;
    /// The prime modulus as a Long type
    const QL: Self::Long;
    /// The prime modulus as a LongLong type
    const QLL: Self::LongLong;

    /// Shift amount for Barrett reduction
    const BARRETT_SHIFT: usize;
    /// Multiplier for Barrett reduction
    const BARRETT_MULTIPLIER: Self::LongLong;

    /// Reduce a value that is at most 2*Q-1 to the range [0, Q)
    fn small_reduce(x: Self::Int) -> Self::Int;
    /// Reduce a larger value using Barrett reduction
    ///
    /// Barrett reduction efficiently computes x mod q without expensive division
    /// It uses the precomputed multiplier to approximate division by q.
    /// Mathematical principle: ⌊x/q⌋ ≈ ⌊(x * ⌊2^k/q⌋) / 2^k⌋
    /// where k = BARRETT_SHIFT is chosen large enough for accuracy
    ///
    /// This replaces slow division with fast operations:
    /// 1. One multiplication: x * BARRETT_MULTIPLIER
    /// 2. One bit shift: >> BARRETT_SHIFT (equivalent to division by 2^k)
    /// 3. One subtraction: x - quotient * q
    ///
    /// The approximation may be off by at most 1, which small_reduce() handles
    fn barrett_reduce(x: Self::Long) -> Self::Int;
}

/// The `define_field` macro creates a zero-sized struct and an implementation of the Field trait
/// for that struct.  The caller must specify:
///
/// * `$field`: The name of the zero-sized struct to be created
/// * `$int`: The primitive integer type to be used to represent members of the field
/// * `$long`: The primitive integer type to be used to represent products of two field members.
///   This type should have roughly twice the bits of `$int`.
/// * `$longlong`: The primitive integer type to be used to represent products of three field
///   members. This type should have roughly four times the bits of `$int`.
/// * `$q`: The prime number that defines the field.
#[macro_export]
macro_rules! define_field {
    ($field:ident, $int:ty, $long:ty, $longlong:ty, $q:literal) => {
        #[doc = concat!(
                    "Finite field marker type for modulus ",
                    stringify!($q),
                    "."
                )]
        #[derive(Copy, Clone, Default, Debug, PartialEq)]
        pub struct $field;

        impl Field for $field {
            type Int = $int;
            type Long = $long;
            type LongLong = $longlong;

            const Q: Self::Int = $q;
            const QL: Self::Long = $q;
            const QLL: Self::LongLong = $q;

            #[allow(clippy::as_conversions)]
            const BARRETT_SHIFT: usize = 2 * (Self::Q.ilog2() + 1) as usize;
            #[allow(clippy::integer_division_remainder_used)]
            // Precomputed ⌊2^k/q⌋ where k = BARRETT_SHIFT
            // This approximates 1/q as a rational number for fast division
            const BARRETT_MULTIPLIER: Self::LongLong = (1 << Self::BARRETT_SHIFT) / Self::QLL;

            fn small_reduce(x: Self::Int) -> Self::Int {
                if x < Self::Q {
                    x
                } else {
                    x - Self::Q
                }
            }

            fn barrett_reduce(x: Self::Long) -> Self::Int {
                let x: Self::LongLong = x.into();
                let product = x * Self::BARRETT_MULTIPLIER;
                let quotient = product >> Self::BARRETT_SHIFT;
                let remainder = x - quotient * Self::QLL;
                Self::small_reduce(Truncate::truncate(remainder))
            }
        }
    };
}

/// An `Elem` is a member of the specified prime-order field.  Elements can be added,
/// subtracted, multiplied, and negated, and the overloaded operators will ensure both that the
/// integer values remain in the field, and that the reductions are done efficiently.  For
/// addition and subtraction, a simple conditional subtraction is used; for multiplication,
/// Barrett reduction.
#[derive(Copy, Clone, Default, Debug, PartialEq)]
pub struct Elem<F: Field>(pub F::Int);

impl<F: Field> Elem<F> {
    /// Create a new field element from an integer value
    pub const fn new(x: F::Int) -> Self {
        Self(x)
    }
}

#[cfg(feature = "zeroize")]
impl<F: Field> Zeroize for Elem<F>
where
    F::Int: Zeroize,
{
    fn zeroize(&mut self) {
        self.0.zeroize();
    }
}

impl<F: Field> Neg for Elem<F> {
    type Output = Elem<F>;

    fn neg(self) -> Elem<F> {
        Elem(F::small_reduce(F::Q - self.0))
    }
}

impl<F: Field> Add<Elem<F>> for Elem<F> {
    type Output = Elem<F>;

    fn add(self, rhs: Elem<F>) -> Elem<F> {
        Elem(F::small_reduce(self.0 + rhs.0))
    }
}

impl<F: Field> Sub<Elem<F>> for Elem<F> {
    type Output = Elem<F>;

    fn sub(self, rhs: Elem<F>) -> Elem<F> {
        Elem(F::small_reduce(self.0 + F::Q - rhs.0))
    }
}

impl<F: Field> Mul<Elem<F>> for Elem<F> {
    type Output = Elem<F>;

    fn mul(self, rhs: Elem<F>) -> Elem<F> {
        let lhs: F::Long = self.0.into();
        let rhs: F::Long = rhs.0.into();
        let prod = lhs * rhs;
        Elem(F::barrett_reduce(prod))
    }
}

/// A `Polynomial` is a member of the ring `R_q = Z_q[X] / (X^256)` of degree-256 polynomials
/// over the finite field with prime order `q`.  Polynomials can be added, subtracted, negated,
/// and multiplied by field elements.  We do not define multiplication of polynomials here.
#[derive(Clone, Default, Debug, PartialEq)]
pub struct Polynomial<F: Field>(pub Array<Elem<F>, U256>);

impl<F: Field> Polynomial<F> {
    /// Create a new polynomial from an array of field elements
    pub const fn new(x: Array<Elem<F>, U256>) -> Self {
        Self(x)
    }
}

#[cfg(feature = "zeroize")]
impl<F: Field> Zeroize for Polynomial<F>
where
    F::Int: Zeroize,
{
    fn zeroize(&mut self) {
        self.0.zeroize();
    }
}

impl<F: Field> Add<&Polynomial<F>> for &Polynomial<F> {
    type Output = Polynomial<F>;

    fn add(self, rhs: &Polynomial<F>) -> Polynomial<F> {
        Polynomial(
            self.0
                .iter()
                .zip(rhs.0.iter())
                .map(|(&x, &y)| x + y)
                .collect(),
        )
    }
}

impl<F: Field> Sub<&Polynomial<F>> for &Polynomial<F> {
    type Output = Polynomial<F>;

    fn sub(self, rhs: &Polynomial<F>) -> Polynomial<F> {
        Polynomial(
            self.0
                .iter()
                .zip(rhs.0.iter())
                .map(|(&x, &y)| x - y)
                .collect(),
        )
    }
}

impl<F: Field> Mul<&Polynomial<F>> for Elem<F> {
    type Output = Polynomial<F>;

    fn mul(self, rhs: &Polynomial<F>) -> Polynomial<F> {
        Polynomial(rhs.0.iter().map(|&x| self * x).collect())
    }
}

impl<F: Field> Neg for &Polynomial<F> {
    type Output = Polynomial<F>;

    fn neg(self) -> Polynomial<F> {
        Polynomial(self.0.iter().map(|&x| -x).collect())
    }
}

/// A `Vector` is a vector of polynomials from `R_q` of length `K`.  Vectors can be
/// added, subtracted, negated, and multiplied by field elements.
#[derive(Clone, Default, Debug, PartialEq)]
pub struct Vector<F: Field, K: ArraySize>(pub Array<Polynomial<F>, K>);

impl<F: Field, K: ArraySize> Vector<F, K> {
    /// Create a new vector from an array of polynomials
    pub const fn new(x: Array<Polynomial<F>, K>) -> Self {
        Self(x)
    }
}

#[cfg(feature = "zeroize")]
impl<F: Field, K: ArraySize> Zeroize for Vector<F, K>
where
    F::Int: Zeroize,
{
    fn zeroize(&mut self) {
        self.0.zeroize();
    }
}

impl<F: Field, K: ArraySize> Add<&Vector<F, K>> for &Vector<F, K> {
    type Output = Vector<F, K>;

    fn add(self, rhs: &Vector<F, K>) -> Vector<F, K> {
        Vector(
            self.0
                .iter()
                .zip(rhs.0.iter())
                .map(|(x, y)| x + y)
                .collect(),
        )
    }
}

impl<F: Field, K: ArraySize> Sub<&Vector<F, K>> for &Vector<F, K> {
    type Output = Vector<F, K>;

    fn sub(self, rhs: &Vector<F, K>) -> Vector<F, K> {
        Vector(
            self.0
                .iter()
                .zip(rhs.0.iter())
                .map(|(x, y)| x - y)
                .collect(),
        )
    }
}

impl<F: Field, K: ArraySize> Mul<&Vector<F, K>> for Elem<F> {
    type Output = Vector<F, K>;

    fn mul(self, rhs: &Vector<F, K>) -> Vector<F, K> {
        Vector(rhs.0.iter().map(|x| self * x).collect())
    }
}

impl<F: Field, K: ArraySize> Neg for &Vector<F, K> {
    type Output = Vector<F, K>;

    fn neg(self) -> Vector<F, K> {
        Vector(self.0.iter().map(|x| -x).collect())
    }
}

/// An `NttPolynomial` is a member of the NTT algebra `T_q = Z_q[X]^256` of 256-tuples of field
/// elements.  NTT polynomials can be added and
/// subtracted, negated, and multiplied by scalars.
/// We do not define multiplication of NTT polynomials here.  We also do not define the
/// mappings between normal polynomials and NTT polynomials (i.e., between `R_q` and `T_q`).
#[derive(Clone, Default, Debug, PartialEq)]
pub struct NttPolynomial<F: Field>(pub Array<Elem<F>, U256>);

impl<F: Field> NttPolynomial<F> {
    /// Create a new NTT polynomial from an array of field elements
    pub const fn new(x: Array<Elem<F>, U256>) -> Self {
        Self(x)
    }
}

#[cfg(feature = "zeroize")]
impl<F: Field> Zeroize for NttPolynomial<F>
where
    F::Int: Zeroize,
{
    fn zeroize(&mut self) {
        self.0.zeroize();
    }
}

impl<F: Field> Add<&NttPolynomial<F>> for &NttPolynomial<F> {
    type Output = NttPolynomial<F>;

    fn add(self, rhs: &NttPolynomial<F>) -> NttPolynomial<F> {
        NttPolynomial(
            self.0
                .iter()
                .zip(rhs.0.iter())
                .map(|(&x, &y)| x + y)
                .collect(),
        )
    }
}

impl<F: Field> Sub<&NttPolynomial<F>> for &NttPolynomial<F> {
    type Output = NttPolynomial<F>;

    fn sub(self, rhs: &NttPolynomial<F>) -> NttPolynomial<F> {
        NttPolynomial(
            self.0
                .iter()
                .zip(rhs.0.iter())
                .map(|(&x, &y)| x - y)
                .collect(),
        )
    }
}

impl<F: Field> Mul<&NttPolynomial<F>> for Elem<F> {
    type Output = NttPolynomial<F>;

    fn mul(self, rhs: &NttPolynomial<F>) -> NttPolynomial<F> {
        NttPolynomial(rhs.0.iter().map(|&x| self * x).collect())
    }
}

impl<F> Mul<&NttPolynomial<F>> for &NttPolynomial<F>
where
    F: Field + MultiplyNtt,
{
    type Output = NttPolynomial<F>;

    fn mul(self, rhs: &NttPolynomial<F>) -> NttPolynomial<F> {
        F::multiply_ntt(self, rhs)
    }
}

/// Perform multiplication in the NTT domain.
///
/// This trait is implemented by downstream crates (e.g., ML-DSA, ML-KEM) to define
/// how NTT polynomials are multiplied, since the multiplication rule depends on the
/// specific NTT variant used by each scheme.
pub trait MultiplyNtt: Field {
    /// Multiply two NTT polynomials.
    fn multiply_ntt(lhs: &NttPolynomial<Self>, rhs: &NttPolynomial<Self>) -> NttPolynomial<Self>;
}

impl<F: Field> Neg for &NttPolynomial<F> {
    type Output = NttPolynomial<F>;

    fn neg(self) -> NttPolynomial<F> {
        NttPolynomial(self.0.iter().map(|&x| -x).collect())
    }
}

/// An `NttVector` is a vector of polynomials from `T_q` of length `K`.  NTT vectors can be
/// added and subtracted.  If multiplication is defined for NTT polynomials, then NTT vectors
/// can be multiplied by NTT polynomials, and "multiplied" with each other to produce a dot
/// product.
#[derive(Clone, Default, Debug, PartialEq)]
pub struct NttVector<F: Field, K: ArraySize>(pub Array<NttPolynomial<F>, K>);

impl<F: Field, K: ArraySize> NttVector<F, K> {
    /// Create a new NTT vector from an array of NTT polynomials
    pub const fn new(x: Array<NttPolynomial<F>, K>) -> Self {
        Self(x)
    }
}

#[cfg(feature = "zeroize")]
impl<F: Field, K: ArraySize> Zeroize for NttVector<F, K>
where
    F::Int: Zeroize,
{
    fn zeroize(&mut self) {
        self.0.zeroize();
    }
}

impl<F: Field, K: ArraySize> Add<&NttVector<F, K>> for &NttVector<F, K> {
    type Output = NttVector<F, K>;

    fn add(self, rhs: &NttVector<F, K>) -> NttVector<F, K> {
        NttVector(
            self.0
                .iter()
                .zip(rhs.0.iter())
                .map(|(x, y)| x + y)
                .collect(),
        )
    }
}

impl<F: Field, K: ArraySize> Sub<&NttVector<F, K>> for &NttVector<F, K> {
    type Output = NttVector<F, K>;

    fn sub(self, rhs: &NttVector<F, K>) -> NttVector<F, K> {
        NttVector(
            self.0
                .iter()
                .zip(rhs.0.iter())
                .map(|(x, y)| x - y)
                .collect(),
        )
    }
}

impl<F: Field, K: ArraySize> Mul<&NttVector<F, K>> for &NttPolynomial<F>
where
    for<'a> &'a NttPolynomial<F>: Mul<&'a NttPolynomial<F>, Output = NttPolynomial<F>>,
{
    type Output = NttVector<F, K>;

    fn mul(self, rhs: &NttVector<F, K>) -> NttVector<F, K> {
        NttVector(rhs.0.iter().map(|x| self * x).collect())
    }
}

impl<F: Field, K: ArraySize> Mul<&NttVector<F, K>> for &NttVector<F, K>
where
    for<'a> &'a NttPolynomial<F>: Mul<&'a NttPolynomial<F>, Output = NttPolynomial<F>>,
{
    type Output = NttPolynomial<F>;

    fn mul(self, rhs: &NttVector<F, K>) -> NttPolynomial<F> {
        self.0
            .iter()
            .zip(rhs.0.iter())
            .map(|(x, y)| x * y)
            .fold(NttPolynomial::default(), |x, y| &x + &y)
    }
}

/// A K x L matrix of NTT-domain polynomials.  Each vector represents a row of the matrix, so that
/// multiplying on the right just requires iteration.  Multiplication on the right by vectors
/// is the only defined operation, and is only defined when multiplication of NTT polynomials
/// is defined.
#[derive(Clone, Default, Debug, PartialEq)]
pub struct NttMatrix<F: Field, K: ArraySize, L: ArraySize>(pub Array<NttVector<F, L>, K>);

impl<F: Field, K: ArraySize, L: ArraySize> NttMatrix<F, K, L> {
    /// Create a new NTT matrix from an array of NTT vectors
    pub const fn new(x: Array<NttVector<F, L>, K>) -> Self {
        Self(x)
    }
}

impl<F: Field, K: ArraySize, L: ArraySize> Mul<&NttVector<F, L>> for &NttMatrix<F, K, L>
where
    for<'a> &'a NttPolynomial<F>: Mul<&'a NttPolynomial<F>, Output = NttPolynomial<F>>,
{
    type Output = NttVector<F, K>;

    fn mul(self, rhs: &NttVector<F, L>) -> NttVector<F, K> {
        NttVector(self.0.iter().map(|x| x * rhs).collect())
    }
}

#[cfg(test)]
mod test {
    use super::*;
    // Test types available for future use
    // use hybrid_array::typenum::{U2, U3};

    // Define a simple test field for testing
    define_field!(TestField, u32, u64, u128, 17);

    #[test]
    fn field_arithmetic() {
        // Test basic field operations
        let a = Elem::<TestField>::new(5);
        let b = Elem::<TestField>::new(7);

        // Addition
        let sum = a + b;
        assert_eq!(sum.0, 12);

        // Subtraction
        let diff = a - b;
        assert_eq!(diff.0, 15); // 5 - 7 = -2 ≡ 15 (mod 17)

        // Multiplication
        let prod = a * b;
        assert_eq!(prod.0, 1); // 5 * 7 = 35 ≡ 1 (mod 17)

        // Negation
        let neg_a = -a;
        assert_eq!(neg_a.0, 12); // -5 ≡ 12 (mod 17)
    }

    #[test]
    fn field_reduction() {
        // Test small reduction
        assert_eq!(TestField::small_reduce(16), 16);
        assert_eq!(TestField::small_reduce(17), 0);
        assert_eq!(TestField::small_reduce(18), 1);

        // Test Barrett reduction
        assert_eq!(TestField::barrett_reduce(35), 1); // 35 ≡ 1 (mod 17)
        assert_eq!(TestField::barrett_reduce(34), 0); // 34 ≡ 0 (mod 17)
    }

    #[test]
    fn polynomial_arithmetic() {
        let mut p1_coeffs = Array::default();
        let mut p2_coeffs = Array::default();

        // Set some test coefficients
        p1_coeffs[0] = Elem::<TestField>::new(1);
        p1_coeffs[1] = Elem::<TestField>::new(2);
        p2_coeffs[0] = Elem::<TestField>::new(3);
        p2_coeffs[1] = Elem::<TestField>::new(4);

        let p1 = Polynomial::new(p1_coeffs);
        let p2 = Polynomial::new(p2_coeffs);

        // Addition
        let sum = &p1 + &p2;
        assert_eq!(sum.0[0].0, 4);
        assert_eq!(sum.0[1].0, 6);

        // Subtraction
        let diff = &p1 - &p2;
        assert_eq!(diff.0[0].0, 15); // 1 - 3 = -2 ≡ 15 (mod 17)
        assert_eq!(diff.0[1].0, 15); // 2 - 4 = -2 ≡ 15 (mod 17)
    }
}
